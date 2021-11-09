#!/usr/bin/env python

# Copyright (C) 2016  Shengwei Hou
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function
import os
import sys
import argparse
from Bio import SeqIO
import subprocess


rev_dict = {"A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N"}


class Fasta(object):

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq.upper()

    def get_gc_content(self):
        ret = 0
        if len(self.seq) == 0:
            return ret
        for c in self.seq:
            if c == "G" or c == "C":
                ret += 1
        return float(ret)/len(self.seq)

    def __str__(self):
        return ">"+self.header+"\n"+self.seq+"\n"


def parse_fasta(fasta_file):
    """
    :param fasta_file: input fasta file
    :return:           yield Fasta record as a generator
    """
    header = ""
    seq = []
    with open(fasta_file, "r") as ih:
        for line in ih:
            if line.startswith(">"):
                if header:
                    yield Fasta(header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield Fasta(header, "".join(seq))


def get_pretty_seq(seq, width):
    ret = ""
    for i, char in enumerate(seq):
        if i > 0 and i % width == 0:
            ret += "\n"
        ret += char
    return ret


def run_CD_HIT_2D(input_faa1, input_faa2, out_file, identity=1.0, word_size=5):

    cmd = "cd-hit-2d -i {first_faa} -i2 {second_faa} -o {out_file} -c {identity} -n {word_size}".format(
        first_faa=input_faa1,
        second_faa=input_faa2,
        out_file=out_file,
        identity=identity,
        word_size=word_size
    )

    proc = subprocess.Popen(cmd, shell=True)
    proc.wait()

    return out_file + ".clstr"


def parse_cd_hit_cluster_result(input_cd_hit_cluster_file):
    """
    :param input_cd_hit_cluster_file: cd-hit clstr file
    :return: a list of sequentially identical locus_tag list, the length could be 1/2/3/...
    """

    # get sequentially identical pairs
    pair_lists = [] # [[BSR30_1, HG001_1], [BSR30_2, HG001_2]]
    with open(input_cd_hit_cluster_file, "r") as ih:
        pair = []
        for line in ih:
            if line.startswith(">"):
                if pair:
                    pair_lists.append(pair)
                pair = []
            else:
                line = line.strip().split(">")[1].split("...")[0]
                pair.append(line)

    return pair_lists


def parse_locustag2coordinates(input_faa1, input_faa2):
    """
    :param input_faa1:  first faa file contains genomic coordinates in the headers
    :param input_faa2: second faa file contains genomic coordinates in the headers
    :return: a locustag2coordinates dict
    """
    locustag2coordinates = {}

    for input_faa in [input_faa1, input_faa2]:
        for fasta in parse_fasta(input_faa):
            header = fasta.header
            locus_tag = header.strip().split()[0]
            coordinate = header.strip().split("coordinate=")[-1]
            locustag2coordinates[locus_tag] = coordinate

    return locustag2coordinates


def refine_identical_pairs_by_coordinates(cd_hit_identical_pairs, locustag2coordinates, ignore_seqname=True):

    # get coordinate lists
    coordinate_lists = []
    pair_lists = cd_hit_identical_pairs
    for pairs in pair_lists:
        coordinate_list = []
        for locustag in pairs:
            print(locustag)
            coordinate = locustag2coordinates[locustag]
            if ignore_seqname:
                coordinate = coordinate.split(":")[-1]
            coordinate_list.append(coordinate)
        coordinate_lists.append(coordinate_list)

    # refine identical pairs with coordinates
    refined_identical_pairs = []
    for pair_list, coordinate_list in zip(pair_lists, coordinate_lists):
        if len(pair_list) == 1:
            refined_identical_pairs.append(pair_list)
        elif len(pair_list) == 2:
            start_1, end_1, strand_1 = coordinate_list[0].split("_")
            start_2, end_2, strand_2 = coordinate_list[1].split("_")
            if strand_1 == strand_2:
                if int(start_1) <= int(start_2) <= int(end_1) or int(start_1) <= int(end_2) <= int(end_1):
                    refined_identical_pairs.append(pair_list)
                else:
                    print(pair_list)
        # for more than 3 orthologs, compare each with the first one
        else:
            start_1, end_1, strand_1 = coordinate_list[0].split("_")
            for i, other_coordinate in enumerate(coordinate_list[1:]):
                start_i, end_i, strand_i = other_coordinate.split("_")
                if int(start_1) <= int(start_i) <= int(end_1) or int(start_1) <= int(end_i) <= int(end_1):
                    refined_identical_pairs.append([pair_list[0], pair_list[i+1]])
                    break

    return refined_identical_pairs



def pair_locus_tags(input_faa1, input_faa2, out_file):

    # get cd-hit cluster file
    cd_hit_result = "cd_hit_results.txt"
    cd_hit_cluster = run_CD_HIT_2D(input_faa1, input_faa2, out_file=cd_hit_result, identity=1.0, word_size=5)

    # cd_hit identical pairs
    cd_hit_identical_pairs = parse_cd_hit_cluster_result(cd_hit_cluster)

    # get locustag2coordinates
    locustag2coordinates = parse_locustag2coordinates(input_faa1, input_faa2)

    # get coordinate constraint identical pairs
    refined_identical_pairs = refine_identical_pairs_by_coordinates(cd_hit_identical_pairs, locustag2coordinates)

    # write out identical pairs
    with open(out_file, "w") as oh:
        for pair_list in refined_identical_pairs:
            if len(pair_list) == 1:
                oh.write(pair_list[0] + "\tNA\n")
            else:
                oh.write(pair_list[0] +"\t"+ pair_list[1] +"\n")


def main():

    # main parser
    parser = argparse.ArgumentParser(description="pair locus_tag names between two faa files of same genome \
                         using genomic coordinates and sequence identities")
    parser.add_argument("input_faa1", help="input the first faa file")
    parser.add_argument("input_faa2", help="input the second faa file")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    if args.prefix:
        prefix = args.prefix
        out_file = os.path.join(args.out_folder, prefix+".tsv")
    else:
        basename1 = os.path.basename(args.input_faa1)
        basename2 = os.path.basename(args.input_faa2)
        prefix = basename1 +"_"+ basename2
        out_file = os.path.join(args.out_folder, prefix+".tsv")

    if os.path.exists(out_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # pair
    pair_locus_tags(args.input_faa1, args.input_faa2, out_file)


if __name__ == "__main__":
    main()

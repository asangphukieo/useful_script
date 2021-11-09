#!/usr/bin/env python

# Copyright (C) 2018  Shengwei Hou
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


def ipr2goatools(input_ipr_tsv, out_file):
    """ slice interproscan tsv file to get protein to GO lookup table
    """
    gene2GOs = {} # {gene1:[GO1, GO2, ...]}
    with open(input_ipr_tsv, "r") as ih:
        for line in ih:
            if "GO:" not in line:
                continue
            line = line.strip().split()
            gene = line[0]
            GO_terms = []
            for _ in line[1:]:
                if "GO:" in _:
                    GO_terms.extend(_.split("|"))
            if gene not in gene2GOs:
                gene2GOs[gene] = GO_terms
            else:
                gene2GOs[gene].extend(GO_terms)

    with open(out_file, "w") as oh:
        for gene, GOs in gene2GOs.items():
            oh.write(gene +"\t"+ ";".join(set(GOs)) +"\n")


def main():


    # main parser
    parser = argparse.ArgumentParser(description="slice interproscan tsv file \
                         to get protein to GO lookup table")
    parser.add_argument("input_ipr_tsv", help="input interproscan file in tsv format")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    if args.prefix:
        prefix = args.prefix
        out_file = os.path.join(args.out_folder, prefix+"_ipr2goatools.tsv")
    else:
        basename = os.path.basename(args.input_ipr_tsv)
        prefix = basename.split(".")[0]
        out_file = os.path.join(args.out_folder, prefix+"_ipr2goatools.tsv")

    if os.path.exists(out_file):
        if args.force:
            sys.stdout.write("Warning: output file exists, will be overwriten!")
        else:
            sys.stderr.write("Error: output file detected, please backup it at first")
            sys.exit(0)

    # convert
    ipr2goatools(args.input_ipr_tsv, out_file)


if __name__ == "__main__":
    main()

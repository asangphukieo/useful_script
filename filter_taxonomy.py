#!/usr/bin/env python3

#modify from require creteJSON-1.0 

import sys,os
from Bio import SeqIO

import math
from ete3 import NCBITaxa
import sys
import os
import warnings
warnings.filterwarnings('ignore', '.*was translated into.*',)


# absolute path of this script
scriptPath = os.path.dirname(os.path.realpath(__file__))

# absolute path to local NCBI taxa DB file
NCBITaxaDbFile = scriptPath + "/taxa.sqlite"

def anlyse_input_file(in_file_new_line, in_file_glassgo):
    unique_ids = dict()
    counter = 0
    if in_file_new_line == "" and in_file_glassgo == "":
        print("Please specify your input!")
        exit()
    if in_file_new_line != "" and in_file_glassgo == "":
        handle = open(in_file_new_line, "r")
        for tax_id in handle:
            tax_id = tax_id.rstrip()
            if tax_id in unique_ids:
                unique_ids[tax_id] += 1
            else:
                unique_ids[tax_id] = 1
            counter += 1
        handle.close()
    if in_file_new_line == "" and in_file_glassgo != "":
        handle = open(in_file_glassgo, "r")
        for line in handle:
            if line.startswith(">"):
                tax_id_all = line.split("-taxID:")
                if len(tax_id_all) == 2:
                    tax_id = tax_id_all[1].split(";")[0]
                    if tax_id in unique_ids:
                        unique_ids[tax_id] += 1
                    else:
                        unique_ids[tax_id] = 1
                    counter += 1
        handle.close()
    if counter == 0:
        #sys.stdout.write("*** ERROR: No tax-ids were found! Please check your input file! ***\n")
        return '*** ERROR: No tax-ids were found! Please check your input file! ***'
        #exit()
    return unique_ids


def compute_taxid_paths(unique_tax_id_hash, ):
    #ncbi = NCBITaxa()

    ncbi = NCBITaxa(NCBITaxaDbFile)
    pathways = list()
    tax_name_ctr = dict()
    max_scalable_hits = 1000
    max_value = 40
    for tax_id in unique_tax_id_hash:
        # save mode; because the tax id can also be a not parsable string
        try:
            # get pathway (ete3 package) => "['root', 'bacteria', 'bac1']"
            global_scaling_val = unique_tax_id_hash[tax_id]
            lineage = ncbi.get_lineage(int(tax_id))

            # prepare output for CopraRNA
            path_output1 = ncbi.get_rank(lineage)
            path_output2 = lineage

            names = ncbi.get_taxid_translator(lineage)
            tmp_path = list()
            for tax_id2 in lineage:
                tax_name = str(tax_id2) + ":" + str(names[tax_id2])
                if tax_name in tax_name_ctr:
                    tax_name_ctr[tax_name][0] += global_scaling_val
                else:
                    tax_name_ctr[tax_name] = list()
                    tax_name_ctr[tax_name].append(global_scaling_val)
                    #tax_name_ctr[tax_name][0] += unique_tax_id_hash[tax_id]
                    tax_name_ctr[tax_name].append(0)
                    tax_name_ctr[tax_name].append(0)
                tmp_path.append(tax_name)
            # normalize node values
            for tax_name in tax_name_ctr:
                if (tax_name_ctr[tax_name][0]) <= max_scalable_hits:
                    tax_name_ctr[tax_name][1] = math.sqrt(float(tax_name_ctr[tax_name][0])) * 1.26
                    tax_name_ctr[tax_name][2] = "passed"
                else:
                    tax_name_ctr[tax_name][1] = max_value
                    tax_name_ctr[tax_name][2] = "failed"
            # append sub-pathway to pathways
            pathways.append(tmp_path)
        except ValueError:
            path_output1=0
            path_output2=0
            pass

    return pathways, tax_name_ctr, path_output1, path_output2

#####################################################
#ncbi = NCBITaxa(NCBITaxaDbFile)
#ncbi.update_taxonomy_database()
#exit()

input_file=sys.argv[1]
tax_id_filter=sys.argv[2] #cyano phylum = 1117
out_file=input_file #org_tax_filtering.fasta

count=1
count_pass=0
out_write=''
print(input_file)
for sRNA in SeqIO.parse(input_file, "fasta"):

	header = sRNA.description
	seq = sRNA.seq
	#seq_rec.description
	#print(header)
	#print('>'+header+'\n'+seq+'\n')
	w2file=open(str(count)+'_tax_filtering_fasta.temp','w')
	w2file.write('>'+str(header)+'\n'+str(seq)+'\n')
	w2file.close()

	unique_ids = anlyse_input_file("", str(count)+'_tax_filtering_fasta.temp')
	#print(unique_ids)
	if count == 1: #first sequence id is GLASSgo query file, skip it
		out_write+='>'+str(header)+'\n'+str(seq)+'\n'
		#w2file=open('org_tax_filtering.fasta','a')
		#w2file.write('>'+str(header)+'\n'+str(seq)+'\n')  
		count_pass+=1  	
	else:
		if 'ERROR: No tax-ids were found!' not in unique_ids:
			master_paths, tax_name_ctr, path_output_1, path_output_2 = compute_taxid_paths(unique_ids)

			if path_output_1 != 0:
				#print(path_output1)
				#dict_map=path_output1

				try:
					if int(tax_id_filter) in path_output_1:
						print(path_output_1)
						out_write+='>'+str(header)+'\n'+str(seq)+'\n'
						#w2file=open('org_tax_filtering.fasta','a')
						#w2file.write('>'+str(header)+'\n'+str(seq)+'\n')
						count_pass+=1					    	    			
				except KeyError:
					pass
	os.system('rm '+str(count)+'_tax_filtering_fasta.temp')

	count+=1

os.system('mv '+input_file+' '+input_file+'_x') #change the previous file name 

w2file=open(out_file,'w') #replace the previous file
w2file.write(out_write)
w2file.close()

print('From', count, 'sequences' , 'there are', count_pass ,'sequences passing taxonmy filter')

#usage: python filter_taxonomy.py Test_GG_1_modi.fasta 1117
#usage: python filter_taxonomy.py Test_GG_1_modi.fasta 1117

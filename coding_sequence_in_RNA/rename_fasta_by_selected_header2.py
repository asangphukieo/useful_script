import sys
from Bio import SeqIO
#require biopython

#input
#>BA000004_c2762887_1

#output
#>B_halodurans BA000004_c2762887_1

def select_fasta_from_name(name_list,fasta,outfile): #name_list,fasta,outfile
	dict_name={}
	for code in open(name_list):
		code=code.rstrip()
		#print(code)
		code=code.split('\t') #split tab to select the second column which is accession id	
		dict_name[code[1]]=code[0]
	out_file=open(outfile,'w')				
	for seq_rec in SeqIO.parse(fasta,"fasta"):
		header=str(seq_rec.description).split('_')
		modify_name=dict_name[header[0]].split(' ')[0]+'_'+dict_name[header[0]].split(' ')[1]
		header2=modify_name+'_'+str(seq_rec.description)

		out_file.write('>'+str(header2)+'\n'+str(seq_rec.seq)+'\n')
	out_file.close()


select_fasta_from_name(sys.argv[1],sys.argv[2],sys.argv[3])

#usage:
#python rename_fasta_by_selected_header2.py list_included_genomes.txt u_orfm_SR1P.fasta re_u_orfm_SR1P.fasta
#input list_included_genomes >> 2 columns separated by tab , which the second will be used as a key for fasta selecion

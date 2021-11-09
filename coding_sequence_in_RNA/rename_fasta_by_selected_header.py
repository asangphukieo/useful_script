import sys
from Bio import SeqIO
#require biopython

#input
#>lcl|ORF1_BA000004_c2762887:61:186 unnamed protein product, partial
#output
#>Bacillus halodurans C-125 BA000004_c2762887

def select_fasta_from_name(name_list,fasta,outfile): #name_list,fasta,outfile
	dict_name={}
	for code in open(name_list):
		code=code.rstrip()
		#print(code)
		code=code.split('\t') #split tab to select the seconf colund which is accession id	
		dict_name[code[1]]=code[0]
	out_file=open(outfile,'w')				
	for seq_rec in SeqIO.parse(fasta,"fasta"):
		header=str(seq_rec.description).split('_')
		header2=dict_name[header[1]]+' '+header[1]+'_'+header[2].split(':')[0]
		header2=header2.replace('(',' ')
		header2=header2.replace(')',' ')
		out_file.write('>'+str(header2)+'\n'+str(seq_rec.seq)+'\n')
	out_file.close()


select_fasta_from_name(sys.argv[1],sys.argv[2],sys.argv[3])

#usage:
#python select_fasta_by_header.py list_included_genomes.txt sr1_de-duplicated.fasta filtered_sr1_de-duplicated.fasta
#input list_included_genomes >> 2 columns separated by tab , which the second will be used as a key for fasta selecion

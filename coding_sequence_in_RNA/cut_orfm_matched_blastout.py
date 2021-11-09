from Bio import SeqIO
import sys
import subprocess
from subprocess import Popen, STDOUT, PIPE, call

# select possible ORF from blast output
# cut blast match from starting site that matched to SR1P till the end of ORF
# ------------------------ Blast matched , ORF 
#      --------------- 	   SR1P
#      ------------------- Example Output

#### example of input files ####
# orfm output
# >CP000764_c2753285_25_1_2
# DHSIRGKRKMGTIVCQDCESTIAYFEEEKTTVLYGKCGSHCECGHEKHVKA

# blast output
# qseqid stitle evalue pident qstart qend sstart send sseq
# tr|A0A0C3FE40|A0A0C3FE40_BACIU	CP000764_c2753285_25_1_2	6.51e-13	64.706	1	34	10	43
################################

def cut_orfm_matched_2_blast(orfm_fasta,blast_output,outfile): 

	#make list of orf
	allorf_matched=Popen('cut -f2,7 '+blast_output,shell=True,stdout=PIPE,stderr=STDOUT,encoding='utf8').stdout.read().split('\n') #,encoding='utf8'
	list_hits=[]
	for hit in allorf_matched:
		if hit != '':
			hit=hit.rstrip()
			hit=hit.split()
			list_hits.append(hit)

	#print(list_hits)

	out_file=open(outfile,'w')

	for seq_rec in SeqIO.parse(orfm_fasta,"fasta"):
		org_name=''
		#print(seq_rec.description)
		header=seq_rec.description.rstrip()
		
		for i in list_hits:
			if header == i[0] :
				#print('Found') 
				out_file.write('>'+str(header)+'\n'+str(seq_rec.seq)[int(i[1])-1:]+'\n')
				#print(full_header)
		
	out_file.close()
	print('Select fasta file ... Done!')
cut_orfm_matched_2_blast(sys.argv[1],sys.argv[2],sys.argv[3])
#python cut_orfm_matched_blastout.py orfm_re_renamed_filtered_sr1_de-duplicated.fasta u_SR1P.blastout orfm_SR1P.fasta

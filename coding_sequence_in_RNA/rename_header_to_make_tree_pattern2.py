
from Bio import SeqIO
import sys

#input format
#>BA000004.3:c2762887-2762701 Bacillus halodurans C-125 DNA, complete genome-p.c.VAL:57.74%-taxID:272558
#output 
#>B_halodurans_BA000004_c2762887

def rename_seq(fasta,outfile): 

	out_file=open(outfile,'w')

	for seq_rec in SeqIO.parse(fasta,"fasta"):
		org_name=''
		if ',' in seq_rec.description:
			header=seq_rec.description.split(',')
			header2=header[0].split()
			acc=header2[0].split('.')[0]		
		elif '-p.c.VAL' in seq_rec.description:
			header=seq_rec.description.split('-p.c.VAL')
			header2=header[0].split()
			acc=header2[0].split('.')[0]
		
		count=0
		locus=''
		for i in header2:
			if count ==1:
				org_name+=i.split()[0]
			elif count ==2:
				org_name+='_'+i
			elif count ==0:
				locus=i.split('-')[0].split(':')[1]
				
			count+=1
		org_name=''+acc+'_'+locus

		full_header=org_name
		full_header=full_header.replace('DNA','')
		full_header=full_header.replace('complete genome','')
		full_header=full_header.replace('genome','')
		full_header=full_header.replace('chromosome','')
		full_header=full_header.replace('  ',' ')
		out_file.write('>'+str(full_header)+'\n'+str(seq_rec.seq)+'\n')
		
	out_file.close()

rename_seq(sys.argv[1],sys.argv[2])
#python rename_header_to_make_tree.py ../comb_glassgo_renamed/pass_glassgo2copra.fasta ../comb_glassgo_renamed/rename_pass_glassgo2copra.fasta

from Bio import SeqIO
import sys

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
		for i in header2:
			if count >=1:
				org_name+=i+' '
			count+=1
		org_name+='\t'+acc+''

		full_header=org_name.replace('(',' ').replace(')','')
		full_header=full_header.replace('DNA','')
		full_header=full_header.replace('complete genome','')
		full_header=full_header.replace('genome','')
		full_header=full_header.replace('chromosome','')
		full_header=full_header.replace('  ',' ')
		out_file.write('>'+str(full_header)+'\n'+str(seq_rec.seq)+'\n')
		
	out_file.close()

rename_seq(sys.argv[1],sys.argv[2])
#python rename_header_to_make_tree.py ../comb_glassgo_renamed/pass_glassgo2copra.fasta ../comb_glassgo_renamed/rename_pass_glassgo2copra.fasta

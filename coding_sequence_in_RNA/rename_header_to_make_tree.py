from Bio import SeqIO
import sys

#>gi|1040087903|gb|CP016194.1|:c4124154-4123974 Bacillus thuringiensis serovar coreanensis strain ST7, complete genome-p.c.VAL:57.76%-taxID:180843

def rename_seq(fasta,outfile): 

	out_file=open(outfile,'w')

	for seq_rec in SeqIO.parse(fasta,"fasta"):
		org_name=''
		if ',' in seq_rec.description:
			header=seq_rec.description.split(',')
			#print(header)
			header2=header[0].split('|')[4]
			acc=header[0].split('|')[3].split('.')[0]	
			locus=header2.replace(':','').split()[0]
			locus=locus.replace('-','_')
		else:
			header=seq_rec.description.split('-p.c.VAL')
			#print(header)
			header2=header[0].split('|')[4]
			acc=header[0].split('|')[3].split('.')[0]	
			locus=header2.replace(':','').split()[0]
			locus=locus.replace('-','_')
		
		org_name=acc+'_'+locus

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

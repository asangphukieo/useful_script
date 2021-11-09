import sys,os
from Bio import pairwise2, Align, AlignIO, SeqIO
import subprocess

#to cluster similar genes with specific threshold and select representative sequence

def usearch(fasta_file,ident_score,output):
	print('Start usearch ...')
	cmd='../library/usearch11.0.667_i86linux32 -cluster_fast '+fasta_file+' -id '+str(ident_score)+' -centroids '+output+' -uc clusters_data.uc'
	ps = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	output = ps.stdout.read()
	print('usearch finish^^')


def cluster_gene(list_name, fasta_file, ident_score,output):
	dict_name={}
	for i in open(list_name):
		i=i.rstrip()
		dict_name[i]='T'

	selected_seq=open(output+'.temp','w')
	print('Collect selected sequence ...')

	for seq_rec in SeqIO.parse(fasta_file,"fasta"):
		if seq_rec.description in dict_name:
			selected_seq.write('>'+seq_rec.description+'\n'+str(seq_rec.seq)+'\n')
	
	selected_seq.close()	
	
	#use usearch to cluster
	usearch(output+'.temp',ident_score,output)
	os.system('rm '+output+'.temp')

#cluster_gene(list_name, fasta_file, ident_score,output)
cluster_gene(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])

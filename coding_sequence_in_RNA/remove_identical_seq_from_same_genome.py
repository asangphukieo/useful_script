from Bio import SeqIO
import sys
import subprocess
from subprocess import Popen, STDOUT, PIPE, call

#### example of input files ####
#>CP002293_c2818810_18_3_3
#MGTIVCQTCDATIAHFEDEKVTTLYGKCSKCDCSDTKEDEQ
################################

## remove identical sequence and nested sequence


def Sort(sub_li): 
  
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using first element of  
    # sublist lambda has been used 
    sub_li.sort(reverse = True, key = lambda x: x[0]) 
    return sub_li 

def remove_identical_seq(orfm_fasta,outfile): 
	allorf_matched=Popen('grep ">" '+orfm_fasta+' | cut -f1 -d"_" | sed s/">"/""/g| sort -u ',shell=True,stdout=PIPE,stderr=STDOUT,encoding='utf8').stdout.read().split('\n')

	list_hits=[]
	for seq_rec in SeqIO.parse(orfm_fasta,"fasta"):
		header=seq_rec.description.rstrip()
		list_hits.append([len(seq_rec),header,seq_rec.seq])
	#print(list_hits)
	sorted_seq_len=Sort(list_hits)
	#print(sorted_seq_len)

	print('Number of unique genomes ==>',len(allorf_matched)-1)
	count=1
	collected_seq=[]
	for i in allorf_matched:
		if i != '':
			#print(i)
			i=i.rstrip()
			selected_seq=[]
			count_seq=1
			for j in sorted_seq_len: #j is query j[0]= length, j[1]= seq name, j[2]= seq
				#print(j)
				
				if i in j[1]:
					head=j[1].split('_')[0]+'_'+j[1].split('_')[1]+'_'+str(count_seq)
					if len(selected_seq) >= 1: #compare j to existing sequences in seleted_list
						passport='Accept'
						for db in selected_seq: #db[0] = seq name, db[1]=seq
							if (j[2] in db[1]) or (j[2] == db[1]):
								passport='Unaccept'

						if passport=='Accept':
							selected_seq.append([head,j[2]])
							count+=1
							count_seq+=1
						else:
							count+=1
					else:
						#print('Add')
						selected_seq.append([head,j[2]]) 
						count+=1
						count_seq+=1
			
			#check uniqu e#if there is only one copy add *
			if len(selected_seq) == 1:
				
				new_head=selected_seq[0]
				
				selected_seq=[[new_head[0]+'*',new_head[1]]]
			for A in selected_seq:
				collected_seq.append(A)
	print('Total sequences ==>',count)
	out_file=open(outfile,'w')
	for seq_u in collected_seq:	
		out_file.write('>'+str(seq_u[0])+'\n'+str(seq_u[1])+'\n')
	out_file.close()

	print('Select fasta file ... Done!')
remove_identical_seq(sys.argv[1],sys.argv[2])
#python remove_identical_seq_from_same_genome.py orfm_SR1P.fasta u_orfm_SR1P.fasta

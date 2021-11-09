
#build_tree dependency
from Bio import SeqIO
import sys

#remove redundant sequences but keep paralogs 
#keep non-overlaped paralogs
#algorithm: maximize number of paralogs in the genome and keep only non-overlaped paralog. If more than one pattern found (>1 Hit per sRNA), keep longest pattern.

seq_dict={}

for seq_rec in SeqIO.parse(sys.argv[1],"fasta"):
	#print(seq_rec.id)
	pos=seq_rec.id.split(':')
	acc=pos[0]
	#print(len(seq_rec))

	if 'c' in pos[1]: #reverse transcipt
		col=pos[1].replace('c','')
		col=col.split('-')
		str_srna=col[1]
		stop_srna=col[0]
		#print(pos[1],str_srna,stop_srna)

	else:
		col=pos[1].split('-')
		str_srna=col[0]
		stop_srna=col[1]
		#print(pos[1],str_srna,stop_srna)

	if acc not in seq_dict:
		seq_dict[acc]=[[str_srna,stop_srna,seq_rec.description,seq_rec.seq,len(seq_rec)],]

		#print('New data add',acc,seq_dict[acc])
	
	else:
		#print('Redundant',acc,len(seq_dict[acc]))
		#print(seq_dict[acc])
		collect_i=[]
		#print(acc,seq_dict[acc])
		for i in seq_dict[acc]:
			#collect_i.append(i)

			collect_i.append(i) #add existing genes to db
		
		collect_i.append([str_srna,stop_srna,seq_rec.description,seq_rec.seq,len(seq_rec)])
		del seq_dict[acc]
		seq_dict[acc]=collect_i

		#print(seq_dict[acc])

		#print('Update Redundant',acc,len(seq_dict[acc]))
	#seq_dict[seq_rec.id]=seq_rec.seq
	#key.append(seq_rec.id)

collect_k=[]

for k in seq_dict:
	#for only 3 paralogs	
	if len(seq_dict[k]) >=2:
		original_list=seq_dict[k]
		max_len_dblist=0
		selected_dblist=''
		sum_seqlong=0
		for j in range(0,len(seq_dict[k])):
			db_list=[]
			db_list.append(seq_dict[k][j])
			for a in original_list:	
				#print('Y')
				str_srna=a[0]
				stop_srna=a[1]
				len_srna=int(a[4])
				status='accept'

				sum_seq=0
				for n in db_list:#included list
					#print(n,'f')
					db_start=n[0]
					db_stop=n[1]
					sum_seq+=int(n[4])

					if ( int(str_srna) <= int(db_start) < int(stop_srna))  or ( int(str_srna) < int(db_stop) <= int(stop_srna))    or    ( int(db_start) <= int(str_srna) < int(db_stop))  or ( int(db_start) < int(stop_srna) <= int(db_stop)):
						status='Unaccept'
						#db_list.append(A)
							#original_list[j]='REMOVE'
						#print(original_list)
						#else:#overlaped found
							#original_list[j]='REMOVE'
				if status== 'accept':
					db_list.append(a)
				

			if len(db_list) >= max_len_dblist and (sum_seq+len_srna) > sum_seqlong:
				max_len_dblist=len(db_list)
				selected_dblist=db_list
				sum_seqlong = sum_seq+len_srna
		#print(max_len_dblist,selected_dblist)

		del seq_dict[k]
		seq_dict[k]=selected_dblist



w2stat=open('gene_duplicate_stat.temp','w')
w2file=open(sys.argv[2],'w')
for k in seq_dict:
	for j in seq_dict[k]:
		w2file.write('>'+j[2]+'\n'+str(j[3])+'\n')
	w2stat.write(seq_dict[k][0][2].split(' ')[1]+' '+seq_dict[k][0][2].split(' ')[2]+'\t'+k+'\t'+str(len(seq_dict[k]))+'\n')
w2file.close()
w2stat.close()
print('Done!!')

#usage
#python remove_redundant_seq.py ../test_dup.fasta outtest.fasta


				

#!/usr/bin/python

import sys
import commands
import subprocess

#map ref_seq from fasta file to the Copra genome database
#and remove unmatched organisms
#python2 Remove_unmatch_organism.py test.fasta CopraRNA_available_organisms.txt out.fasta

openfasta= open(sys.argv[1],'r')
readfasta=openfasta.read()
openfasta.close()
readfasta=readfasta.split('>')
count_seq=0
count_found=0

print 'Open file >> ',sys.argv[1]
commands.getoutput('mv '+sys.argv[1]+' '+sys.argv[1]+'.previous')

def search_in_list(listRefSeq,query):
	for i in listRefSeq:
		if query in i :
			return i

#create list of RefSeq ID from Copra file
RefSeq=commands.getoutput('cut -f1 '+sys.argv[2])
RefSeq=RefSeq.split('\n')
RefSeqList=[]
for i in RefSeq:
	if ' ' in i:
		i=i.split()
		for j in i:
			RefSeqList.append(j.replace('/',''))
	else:
		RefSeqList.append(i.replace('/',''))

#write output file
write2file=open(sys.argv[1],'w')

if len(readfasta) >0:
	for i in readfasta:
		if i != '':
			line=i.split('\n')# line[0]=gene label, line[1] = nt sequence
			count_seq+=1

			match=search_in_list(RefSeqList,line[0])
			if str(match) != 'None':
				#print ('Sequence',line[0],'Found')
				count_found+=1
				write2file.write('>'+str(match)+'\n'+str(line[1])+'\n')


			print 'Query sequence from GLASSgo is',line[0]

	print 'Total number of fasta sequence is',count_seq
	print 'Total number of fasta sequence found in CopraRNA database is',count_found

else:
	print 'Wrong Fasta file!!!'

write2file.close()
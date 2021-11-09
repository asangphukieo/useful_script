import os
import commands
import sys

def get_uniprot_seq(list_of_uni_ID):
	for i in list_of_uni_ID:
		os.system("wget http://www.uniprot.org/uniprot/"+str(i)+".fasta")
	os.system("cat *.fasta > all_seq")
	os.system("rm *.fasta ")
	os.system("mv all_seq all_seq3.faa")

file_uniprot_id= sys.argv[1] #"uniq_all_secific_photo.out" #"uniq_all_secific_photo.out"
prots = commands.getoutput("cut -f1 "+file_uniprot_id+" ")
prots=prots.split('\n')
get_uniprot_seq(prots)

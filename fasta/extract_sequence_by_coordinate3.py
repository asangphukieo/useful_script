from Bio import SeqIO
import sys
import random as rd

#use conda of build_tree
#allow to get franking region

#coordinate input
#1000027,1000218,+,BSU_09259,yhzG
#1000364,1001716,+,BSU_09260,yhxA

## requirement 
#[1] Genome in FASTA
#[2] Start,Stop position,Strand,locus_tag,gene_name  (in 5 columns separated by tab) 
#[3] 5' UTR length to include
#[4] 3' UTR length to include
#[5] 'Random' or 'NoRandom' start and stop position (keep seq length identical to original coordinate) 
#[6] output file 
#[7] excluding longer sequences (CLIP peak region)
#[8] call anti-sense rna ['T','F']

#usage
#python extract_sequence_by_coordinate3.py FQ312003.fasta Bacillus_subtilis_168_uniq.locus 100 50 'NoRandom' test.fasta 100000000 'T'

inFile = open(sys.argv[1],'r') #input fasta file
coordinate=sys.argv[2]
UTR5=int(sys.argv[3]) 
UTR3=int(sys.argv[4])  #including franking region
random=sys.argv[5]
fw=open(sys.argv[6],'w') #output fasta file
length_filter=int(sys.argv[7]) #filter out very long binding site ,default 75 folling paper
anti=sys.argv[8]

for record in SeqIO.parse(inFile,'fasta'):
	count=1
	for line in open(coordinate): #input coordinate file form blastout

		line= line.rstrip()
		col=line.split()
		start=int(col[0])-1
		end=int(col[1])
		strandS=str(col[2])

		locus_tag=str(col[3])
		gene_name=str(col[4])

		if (end - start) <length_filter:

			if random == 'Random':
				org_len=end-start
				start=rd.randrange(1, len(record))
				end=start+(org_len)

			F_start=start-UTR5
			F_end=end+UTR3	

			fw.write(">" + str(locus_tag) +';'+str(gene_name)+';'+str(F_start)+':'+str(F_end)+';'+str(strandS) +"\n")
			#print(record.seq)
			if strandS=='+':
				fww = (str(record.seq[F_start:F_end])+ '\n')
			else:
				fww = (str(record.seq[F_start:F_end].reverse_complement()) + '\n')
			
			fw.write(fww)

			if anti == 'T':	
				fw.write(">" + str(locus_tag) +'_anti;'+str(gene_name)+';'+str(F_start)+':'+str(F_end)+';'+str(strandS) +"\n")
				#print(record.seq)
				if strandS=='+':
					fww = (str(record.seq[F_start:F_end].reverse_complement()) + '\n')
				else:
					fww = (str(record.seq[F_start:F_end]) + '\n')
				
				fw.write(fww)

			count+=1
fw.close()
inFile.close()


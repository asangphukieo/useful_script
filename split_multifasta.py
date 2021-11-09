import sys,os
from Bio import SeqIO

input_file=sys.argv[1]
output_single_sRNA_fasta=sys.argv[2]
#output_single_sRNA_fasta='./SRNA_INPUT/SPLIT_SRNA'

if os.path.exists(output_single_sRNA_fasta)== True:
    os.system('rm -r '+output_single_sRNA_fasta)
else:
    os.system('mkdir '+output_single_sRNA_fasta)


count=1
for sRNA in SeqIO.parse(input_file, "fasta"):
    header = sRNA.name
    seq = sRNA.seq
    #print(record.id)
    #print('>'+header+'\n'+seq+'\n')
    w2file=open(output_single_sRNA_fasta+'/'+str(count)+'.fasta','w')
    w2file.write('>'+str(header)+'\n'+str(seq)+'\n')
    w2file.close()
    count+=1

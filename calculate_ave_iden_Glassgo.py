import sys

score=[]
for i in open(sys.argv[1]):
	if '>' in i:
		i=i.rstrip()
		
		if 'p.c.VAL:' in i:
			ident=i.split('p.c.VAL:')[1].split('%')[0]
			score.append(float(ident))

sum_score=sum(score)
total_num=len(score)
ave=sum_score/total_num

#print(score)
#print(sum_score)
#print(total_num)
print('average p.ident =>',ave)

import re
import random
import sys
inf1=open(sys.argv[1])#original heterozygous file
total_reads=sys.argv[3] #This number is 186978547 in our sample

list_b=list()
dict_o=dict()
out=open(sys.argv[2],"w")#simulated by binomial sampling
new=1
for line in inf1:
	if new==1:
		new=0
	else:
		mobj=re.match("(.*?)\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t(.*?)\t.*?\t.*",line)
		try:
			tag=mobj.group(1)+'_'+mobj.group(2)
			reads=int(mobj.group(5))
			MAF=float(mobj.group(3))
			bi_type=mobj.group(4)
			if bi_type=='bi':
				each_n=int(reads*(1-MAF))
				for j in range(each_n):
					list_b.append({tag:'L'})
				for j in range(each_n):
					list_b.append({tag:'R'})
		except:
			pass
inf1.close()
mx_r=len(list_b)
print(mx_r)
for j in range(total_reads):
	try:
		num_r=random.randint(0,mx_r)
		tuple_x=list_b[num_r].popitem()
		key=tuple_x[0]
		typ=tuple_x[1]
		if key in dict_o:
			if typ in dict_o[key]:
				dict_o[key][typ]+=1
			else:
				dict_o[key][typ]=1
		else:
			dict_o.update({key:{typ:1}})
	except:
		print(num_r)
for key in dict_o:
	if dict_o[key]['L']>dict_o[key]['R']:
		MAF=dict_o[key]['R']/(dict_o[key]['L']+dict_o[key]['R'])
	else:
		MAF=dict_o[key]['L']/(dict_o[key]['L']+dict_o[key]['R'])
	print ("%s\t%3f"%(key,MAF),file=out)

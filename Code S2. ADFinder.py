import os
import time
import sys,getopt
import re
from Bio import SeqIO
from os.path import getsize
def r_m_l(astr):
	mobj_l=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",astr)
	pos_m_l=int(mobj_l.group(1))  
	cigar_l=mobj_l.group(2)
	mobj1=re.match("(\d+)M(\d+)S$",cigar_l)
	mobj2=re.match("(\d+)S(\d+)M$",cigar_l)
	if mobj1:
		pos_return=pos_m_l+int(mobj1.group(1))
		return(pos_return)
	if mobj2:
		return(pos_m_l)
def addword2dict(thedict, key_a, key_b, val):
	if key_a in thedict:
		if key_b in thedict[key_a]:
			thedict[key_a][key_b].append(val)
		else:
			list_t=list()
			list_t.append(val)
			thedict[key_a].update({key_b: list_t})
	else:
		list_t=list()
		list_t.append(val)
		thedict.update({key_a:{key_b: list_t}})
def addword2dict2(thedict, key_a, key_b, val):
	if key_a in thedict:
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})
def addnum2dict(thedict, key_a, key_b):
	if key_a in thedict:
		if key_b in thedict[key_a]:
			thedict[key_a][key_b]+=1
		else:
			thedict[key_a].update({key_b:1})
	else:
		thedict.update({key_a:{key_b:1}})
def find_overlap(l_seq,r_seq,l_pos,r_pos,l_m,r_s,dis,l_l,l_r):
	for i in range(l_l-10,l_m,-1):
		max_same=l_l-i
		same=0
		for j in range(0,max_same+1):
			if l_seq[i+j:i+j+1]==r_seq[j:j+1]:
				same+=1
				if same<=j-5:
					break
		if same/max_same>0.9:
			alter_seq=l_seq[l_m:]+r_seq[max_same-1:r_s+dis]
			return alter_seq
	return '0'
help_str='ASFinder.py -x <genome reference in fasta> [--s1 <paired-end reads 1> --s2 <paired-end reads 2>]|[-U <Unpaired reads>] -o <output folder>'
gzinf1=''
gzinf2=''
outfile1=''
extend_key=0
ref_file=''
sam_file=''
splice_file=''
fa_file=''
con_flag=''
try:
	opts,args=getopt.getopt(sys.argv[1:],"hx:p:o:U:",["extend","help","s1=","s2=","continue="])
except getopt.GetoptError:
	print (help_str)
	sys.exit(2)
for opt,value in opts:
	if opt in ("-h","--help"):
		print (help_str)
		sys.exit()
	if opt in ("--s1"):
		gzinf1=value
	if opt in ("--s2"):
		gzinf2=value
	if opt in ("-o"):
		outfolder=value
	if opt in ("-p"):
		pnum=value
	else:
		pnum="1"
	if opt in ("-x"):
		fa_file=value
	# if opt in ("--splice"):
		# splice_file=value
	if opt in ("-U"):
		ref_file=value
	if opt in ("--continue"):
		con_flag=value
	# if opt in ("-S"):
		# sam_file=value
	if opt in ("--extend"):
		extend_key=1
	else:
		extend_key=0
if not(((((gzinf1 and gzinf2) or (ref_file)) and fa_file) and outfolder) or (con_flag)):
	print (help_str)
	sys.exit()
else:
	try:
		os.system('mkdir '+ outfolder)
	except:
		pass

if con_flag:
	start_from_splice=0
	start_from_sam=0
	list_dir=os.listdir(con_flag)
	for i in range(len(list_dir)):
		l_path = os.path.join(con_flag,list_dir[i])
		if os.path.isfile(l_path):
			mobj=re.match(".*output\.splice",list_dir[i])
			if mobj:
				size_splice=os.path.getsize(l_path)
				if size_splice==0:
					start_from_sam=1
				else:
					for line in open(l_path):
						mobj_M=re.search("insertion_M",line)
						if mobj_M:
							start_from_splice=1
							break
	if start_from_splice==1:
		inf_sam=open(sam_file)
		dict_rec_m=dict()
		for line in inf_sam:
			mobj=re.search("(.*?)\t.*?\t(.*?)\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",line)
			if mobj:
				tag_r=mobj.group(1)
				tag_c=mobj.group(2)
				pos_m=int(mobj.group(3))
				cigar=mobj.group(4)
				seq_r=mobj.group(5)
				mobj5=re.match("(\d+)M$",cigar)
				if mobj5:#total coverage count
					n3=int(mobj5.group(1))
					for i in range(pos_m,pos_m+n3):
						addnum2dict(dict_rec_m,tag_c,i)
		del inf_sam
		inf1=open(splice_file)
		out_ins_file=outfolder + '/' + 'output.insertion.tab'
		out_del_file=outfolder + '/' + 'output.deletion.tab'
		out_ins=open(out_ins_file,"w")
		out_del=open(out_del_file,"w")
	if start_from_sam==1:
		final_outfile=sam_file
		# final_outfile=outfolder + '/' + 'reads_mapping.sort'
		# samtools_cmd='samtools view -F 4 -b -u -h ' + sam_file + '| samtools sort -O SAM - > ' + final_outfile
		# os.system(samtools_cmd)
		print("#2: Reading SAM file...")
		#muti-thread parameters--------
		# threads=[]
		# queue=queue.Queue()
		# lock=threading.Lock()
		# num=threading.Semaphore(10)
		#------------------------------
		inf1=open(final_outfile)
		out1_file=outfolder + '/' + 'output.splice'
		out2_file=outfolder + '/' + 'output.uncompleted.splice'
		out_ins_file=outfolder + '/' + 'output.insertion.tab'
		out_del_file=outfolder + '/' + 'output.deletion.tab'
		out_uncom_file=outfolder + '/' + 'output.uncompleted_insertion.tab'
		out1e=open(out1_file,"w")
		out2e=open(out2_file,"w")
		out_ins=open(out_ins_file,"w")
		out_del=open(out_del_file,"w")
		out_uncom=open(out_uncom_file,"w")
		dict_rec_l=dict()
		dict_rec_r=dict()
		dict_rec_d=dict()
		dict_rec_m=dict()
		###############################################3
		print("#2 Reading SAM file finished")
		print("data preparation finished")
		t_start = time.clock()
		print("#3: Alternative splicing search part 1")
		for line in inf1:
			mobj=re.search("(.*?)\t.*?\t(.*?)\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",line)
			if mobj:
				tag_r=mobj.group(1)
				tag_c=mobj.group(2)
				pos_m=int(mobj.group(3))
				cigar=mobj.group(4)
				seq_r=mobj.group(5)
				mobj1=re.match("(\d+)S\d+M$",cigar)
				mobj2=re.match("\d+M(\d+)S$",cigar)
				mobj3=re.match("(\d+)M(\d+)D\d+M$",cigar)
				mobj4=re.match("(\d+)M(\d+)I\d+M$",cigar)
				mobj5=re.match("(\d+)M$",cigar)
				if mobj3:#long ref
					m_size=int(mobj3.group(1))
					ins_size=int(mobj3.group(2))
					if ins_size>=15:
						alter_seq=dict_chr[tag_c][pos_m+m_size-1:pos_m+m_size-1+ins_size]
						start=pos_m+m_size
						end=pos_m+m_size+ins_size-1
						print ("deletion_S_S\t%s\t%s\t%d\t%s\t%s\t%s" % (tag_c,tag_r,pos_m,cigar,seq_r,alter_seq),file=out1e)
				if mobj4:#long read
					m_size=int(mobj4.group(1))
					ins_size=int(mobj4.group(2))
					if ins_size>=15:
						insert_point=pos_m+m_size-1
						alter_seq=seq_r[m_size:m_size+ins_size]
						print ("insertion_S_S\t%s\t%s\t%d\t%s\t%s\t%s" % (tag_c,tag_r,pos_m,cigar,seq_r,alter_seq),file=out1e)
				if mobj1:
					n1=int(mobj1.group(1))
					if n1>20:
						addword2dict(dict_rec_d,tag_c,tag_r,line)
				if mobj2:
					n2=int(mobj2.group(1))
					if n2>20:
						addword2dict(dict_rec_d,tag_c,tag_r,line)
				if mobj5:#total coverage count
					n3=int(mobj5.group(1))
					for i in range(pos_m,pos_m+n3):
						addnum2dict(dict_rec_m,tag_c,i)
	###################################################################
		del inf1
		print("Data structure initiation for large deletions")
		for key1 in dict_rec_d:
			list_keys=list(dict_rec_d[key1].keys())
			for key2 in list_keys:
				if len(dict_rec_d[key1][key2])==1:
					mobj=re.search(".*?\t.*?\t.*?\t.*?\t.*?\t(.*?)\t.*?\t.*?\t.*?\t.*?\t",dict_rec_d[key1][key2][0])
					if mobj:
						cigar=mobj.group(1)
						mobj1=re.match("(\d+)S(\d+)M$",cigar)
						mobj2=re.match("(\d+)M(\d+)S$",cigar)
						if mobj1:
							addword2dict2(dict_rec_r,key1,key2,dict_rec_d[key1][key2][0])
						if mobj2:
							addword2dict2(dict_rec_l,key1,key2,dict_rec_d[key1][key2][0])
					dict_rec_d[key1].pop(key2)
	###################################################################
		print("Searching for large deletions")
		for key1 in dict_rec_d:
			for key2 in dict_rec_d[key1]:
				left2_list_tag=list()
				right2_list_tag=list()
				left2_list_pos=list()
				right2_list_pos=list()
				left2_list_M=list()
				right2_list_M=list()
				left2_list_S=list()
				right2_list_S=list()
				left2_list_seq=list()
				right2_list_seq=list()
				left2_list_cigar=list()
				right2_list_cigar=list()
				#Formation of data structure (reads that have more than 2 matches)
				for vk in dict_rec_d[key1][key2]:
					mobj=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",vk)
					if mobj:
						pos_m=int(mobj.group(1))
						cigar=mobj.group(2)
						seq_r=mobj.group(3)
						mobj1=re.match("(\d+)S(\d+)M$",cigar)
						mobj2=re.match("(\d+)M(\d+)S$",cigar)
						if mobj1:
							SM_S=int(mobj1.group(1))
							SM_M=int(mobj1.group(2))
							right2_list_tag.append(key2)
							right2_list_pos.append(pos_m)
							right2_list_seq.append(seq_r)
							right2_list_M.append(SM_M)
							right2_list_S.append(SM_S)
							right2_list_cigar.append(cigar)
						if mobj2:
							MS_M=int(mobj2.group(1))
							MS_S=int(mobj2.group(2))
							left2_list_tag.append(key2)
							left2_list_pos.append(pos_m)
							left2_list_seq.append(seq_r)
							left2_list_M.append(MS_M)
							left2_list_S.append(MS_S)
							left2_list_cigar.append(cigar)
				#Formation finished
				for i in range(len(left2_list_tag)):
					for j in range(len(right2_list_tag)):
						dis_lr=left2_list_pos[i]+left2_list_M[i]-right2_list_pos[j]
						dis_2=left2_list_M[i]-right2_list_S[j]
						length_read_l=left2_list_M[i]+left2_list_S[i]
						if dis_lr>=0 and dis_lr<=10:#long read
							alter_seq=seq_r[left2_list_M[i]:length_read_l-right2_list_M[j]+dis_lr]
							print ("insertion_S_T\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d" % (key1,left2_list_tag[i],left2_list_pos[i],right2_list_pos[j],left2_list_cigar[i],right2_list_cigar[j],left2_list_seq[i],alter_seq,dis_lr),file=out1e)
						if dis_lr<0 and dis_2<=10 and dis_2>=0:#long ref
							alter_seq=dict_chr[key1][left2_list_pos[i]+left2_list_M[i]-1:right2_list_pos[j]-1+dis_2]
							print ("deletion_S_T\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d" % (key1,left2_list_tag[i],left2_list_pos[i],right2_list_pos[j],left2_list_cigar[i],right2_list_cigar[j],left2_list_seq[i],alter_seq,dis_2),file=out1e)
		del dict_rec_d
		###############################################################################
		print("Searching for large insertions")
		for key1 in dict_rec_l:
			if key1 in dict_rec_r:
				sorted_list_r=sorted(dict_rec_r[key1].keys(),key=lambda x:r_m_l(dict_rec_r[key1][x]))
				sorted_list_l=sorted(dict_rec_l[key1].keys(),key=lambda x:r_m_l(dict_rec_l[key1][x]))
				for key2_1 in sorted_list_l:
					mobj_l=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",dict_rec_l[key1][key2_1])
					pos_m_l=int(mobj_l.group(1))  
					cigar_l=mobj_l.group(2)
					seq_r_l=mobj_l.group(3)
					mobj2=re.match("(\d+)M(\d+)S$",cigar_l)
					if mobj2:
						MS_M=int(mobj2.group(1))
						MS_S=int(mobj2.group(2))
						for key2_2 in sorted_list_r:
							if key2_2 in dict_rec_r[key1]:
								mobj_r=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",dict_rec_r[key1][key2_2])
								pos_m_r=int(mobj_r.group(1)) 
								cigar_r=mobj_r.group(2)
								seq_r_r=mobj_r.group(3)
								mobj3=re.match("(\d+)S(\d+)M$",cigar_r)
								if mobj3:
									SM_M=int(mobj3.group(1))
									SM_S=int(mobj3.group(2))
									length_read_l=MS_M+MS_S
									length_read_r=SM_M+SM_S
									dis_lr1=pos_m_l+MS_M-pos_m_r
									if dis_lr1<(-100):
										break
									if dis_lr1>100:
										dict_rec_r[key1].pop(key2_2)
									if dis_lr1>=0 and dis_lr1<=10:
										alter_seq=find_overlap(seq_r_l,seq_r_r,pos_m_l,pos_m_r,MS_M,SM_S,dis_lr1,length_read_l,length_read_r)
										if not alter_seq=='0':
											print ("insertion_M\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d" % (key1,key2_1,key2_2,pos_m_l,pos_m_r,cigar_l,cigar_r,seq_r_l,seq_r_r,alter_seq,dis_lr1),file=out1e)
											dict_rec_r[key1].pop(key2_2)
											break
										else:
											if MS_M/length_read_l>=0.5 and SM_M/length_read_r>=0.5:
												print ("%s\t%d\t%d\t%s\t%s"%(key1,pos_m_l+MS_M,pos_m_r,key2_1,key2_2),file=out2e)
		print("Releasing RAM")
		del dict_rec_l
		del dict_rec_r
		t_elapsed = (time.clock() - t_start)
		print("Time used:",t_elapsed)
		out1e.close()
		out2e.close()
		inf1=open(out1_file)#completed reads
		inf2=open(out2_file)#uncompleted reads
		uncom_num_dict=dict()
		for line in inf2:
			mobj1=re.match("(.*?)\t(\d+)\t(\d+)\t.*?\t.*",line)
			if mobj1:
				oritag=mobj1.group(1)
				start1=int(mobj1.group(2))
				start2=int(mobj1.group(3))
				this_key=str(start1)+'_'+str(start2)
				addnum2dict(uncom_num_dict,oritag,this_key)
		num_uncom=0
		print("tag\tchromosome\tposition\tsplicing_depth\tnon-splicing_depth",file=out_uncom)
		for key1 in uncom_num_dict:
			mobj_key=re.match("(.*)_(\d+)_(\d+)_(\d+)",key1)
			ori_key=mobj_key.group(1)
			num_key=mobj_key.group(2)
			start_key=int(mobj_key.group(3))
			end_key=int(mobj_key.group(4))
			for key2 in uncom_num_dict[key1]:
				if uncom_num_dict[key1][key2]>=5:
					num_uncom+=1
					mobj1=re.match("(\d+)_(\d+)",key2)
					pos1=int(mobj1.group(1))
					sumi=0
					numi=0
					meani=0
					if key1 in dict_rec_m:
						for i in range(pos1-5,pos1+5):
							if i in dict_rec_m[key1]:
								sumi+=dict_rec_m[key1][i]
								numi+=1
							else:
								numi+=1
						meani=sumi/numi
					else:
						meani=0
					#print("I%d\t%s\t%s\t%s\t%d\t%d\t%3f" % (num_ins,key1,mobj1.group(1),ins_dict[key1][key2],len(ins_dict[key1][key2]),ins_num_dict[key1][key2],meani),file=out_ins)
					print("I_U%d\t%s\t%d\t%d\t%3f" % (num_uncom,ori_key,int(mobj1.group(1))+start_key-1,uncom_num_dict[key1][key2],meani),file=out_uncom)
else:
	print("Reading fasta file")
	dict_chr=dict()
	new_fa=outfolder+'/'+'fasta_index.fa'
	out_new_fa=open(new_fa,"w")
	print("Extension of reads file")
	if extend_key==0:
		for rec in SeqIO.parse(fa_file,"fasta"):
			len_this=len(rec.seq)
			if len_this>=20000:
				seg_n=len_this//10000
				seg_l=int(len_this/seg_n)
				for i in range(seg_n-1):
					new_id=rec.id+'_'+str(i)+'_'+str(i*seg_l+1)+'_'+str((i+1)*seg_l)
					new_seq=str(rec.seq)[i*seg_l:(i+1)*seg_l]
					last_start=(i+1)*seg_l+1
					print(">%s\n%s" % (new_id,new_seq),file=out_new_fa)
					dict_chr.update({new_id:new_seq})
				last_id=rec.id+'_'+str(seg_n-1)+'_'+str(last_start)+'_'+str(len_this)
				last_seq=str(rec.seq)[last_start-1:len_this]
				print(">%s\n%s" % (last_id,last_seq),file=out_new_fa)
				dict_chr.update({last_id:last_seq})
			else:
				last_id=rec.id+'_'+str(0)+'_'+str(1)+'_'+str(len_this)
				print(">%s\n%s" % (last_id,str(rec.seq)),file=out_new_fa)
				dict_chr.update({last_id:str(rec.seq)})
	if extend_key==1:
		for rec in SeqIO.parse(fa_file,"fasta"):
			len_this=len(rec.seq)
			last_id=rec.id+'_'+str(0)+'_'+str(1)+'_'+str(len_this)
			print(">%s\n%s" % (last_id,str(rec.seq)),file=out_new_fa)
			dict_chr.update({last_id:str(rec.seq)})
	print("Reading fasta file finished")
	print("Unwrap fastq file and reference file")
	if (gzinf1 and gzinf2) or ref_file:
		temfile=outfolder + '/' + 'filtered_reads.fq'
		out1=open(temfile,"w")
		if gzinf1 and gzinf2:
			list_inf_s1=gzinf1.strip().split(",")
			list_inf_s2=gzinf2.strip().split(",")
			if len(list_inf_s1)!=len(list_inf_s2):
				print("the paired-end reads file not equal")
			else:
				for i in range(len(list_inf_s1)):
					mobj1=re.match(".*\/(.*)",list_inf_s1[i])
					mobj2=re.match(".*\/(.*)",list_inf_s2[i])
					if mobj1:
						filename1=mobj1.group(1)
						filename2=mobj2.group(1)
					else:
						print ('Fail to read input file')
						sys.exit()
					mobj1_gz=re.match(".*\.gz",filename1)
					mobj2_gz=re.match(".*\.gz",filename2)
					if mobj1_gz:
						gzout1=outfolder + '/' + filename1+'.fq'
						os.system('gunzip -c ' + list_inf_s1[i] + '>' + gzout1)
						teminf1=open(gzout1)
						#os.system('rm '+gzout1)
					else:
						teminf1=open(list_inf_s1[i])
					if mobj2_gz:
						gzout2=outfolder + '/' + filename2+'.fq'
						os.system('gunzip -c ' + list_inf_s2[i] + '>' + gzout2)
						teminf2=open(gzout2)
						#os.system('rm '+gzout2)
					else:
						teminf2=open(list_inf_s2[i])
					for line in teminf1:
						mobj_a1=re.match("(@.*?) (.*)",line)
						mobj_a2=re.match("(\+.*?) (.*)",line)
						if mobj_a1:
							newline=mobj_a1.group(1)+"L"+" "+mobj_a1.group(2)
							print(newline,file=out1)
						elif mobj_a2:
							newline=mobj_a2.group(1)+"L"+" "+mobj_a2.group(2)
							print(newline,file=out1)
						else:
							print(line,file=out1,end="")
					for line in teminf2:
						mobj_a1=re.match("(@.*?) (.*)",line)
						mobj_a2=re.match("(\+.*?) (.*)",line)
						if mobj_a1:
							newline=mobj_a1.group(1)+"R"+" "+mobj_a1.group(2)
							print(newline,file=out1)
						elif mobj_a2:
							newline=mobj_a2.group(1)+"R"+" "+mobj_a2.group(2)
							print(newline,file=out1)
						else:
							print(line,file=out1,end="")
					teminf1.close()
					teminf2.close()
		if ref_file:
			list_inf_u=ref_file.strip().split(",")
			for i in range(len(list_inf_u)):
				mobj3=re.match(".*\/(.*)",list_inf_u[i])
				if mobj3:
					filename=mobj3.group(1)
				else:
					print ('Fail to read input file')
					sys.exit()
				mobj_gz=re.match(".*\.gz",filename)
				if mobj_gz:
					gzout=outfolder + '/' + filename +'.fq'
					os.system('gunzip -c ' + list_inf_u[i] + '>' + gzout)
					teminf=open(gzout)
					#os.system('rm '+gzout)
				else:
					teminf=open(list_inf_u[i])
				for line in teminf:
					mobj=re.match("(@.*?) (.*)",line)
					if mobj:
						newline=mobj.group(1) +" "+mobj.group(2)
						print(newline,file=out1)
					else:
						print(line,file=out1,end="")
				teminf.close()
		out1.close()
		outfile=outfolder + '/' + 'reads_mapping.sam'
		bowtie2_cmd1='bowtie2-build ' + new_fa + " " + outfolder + '/' + 'bowtie2_reference'
		os.system(bowtie2_cmd1)
		bowtie2_cmd2='bowtie2 -p ' + pnum + ' -x ' + outfolder + '/' + 'bowtie2_reference' + " -U " + temfile + " -S " + outfile + " --local -k 5 --ma 3"
		os.system(bowtie2_cmd2)
		final_outfile=outfile
		# final_outfile=outfile + '.sort'
		# samtools_cmd='samtools view -F 4 -b -u -h ' + outfile + '| samtools sort -O SAM - > ' + final_outfile
		# os.system(samtools_cmd)
	if sam_file and not(gzinf1 and gzinf2) and not(ref_file):
		final_outfile=sam_file
		# final_outfile=outfolder + '/' + 'reads_mapping.sort'
		# samtools_cmd='samtools view -F 4 -b -u -h ' + sam_file + '| samtools sort -O SAM - > ' + final_outfile
		# os.system(samtools_cmd)

	print("#2: Reading SAM file...")
	#muti-thread parameters--------
	# threads=[]
	# queue=queue.Queue()
	# lock=threading.Lock()
	# num=threading.Semaphore(10)
	#------------------------------
	inf1=open(final_outfile)
	out1_file=outfolder + '/' + 'output.splice'
	out2_file=outfolder + '/' + 'output.uncompleted.splice'
	out_ins_file=outfolder + '/' + 'output.insertion.tab'
	out_del_file=outfolder + '/' + 'output.deletion.tab'
	out_uncom_file=outfolder + '/' + 'output.uncompleted_insertion.tab'
	out1e=open(out1_file,"w")
	out2e=open(out2_file,"w")
	out_ins=open(out_ins_file,"w")
	out_del=open(out_del_file,"w")
	out_uncom=open(out_uncom_file,"w")
	dict_rec_l=dict()
	dict_rec_r=dict()
	dict_rec_d=dict()
	dict_rec_m=dict()

	print("#2 Reading SAM file finished")
	print("data preparation finished")
	t_start = time.clock()
	print("#3: Alternative splicing search part 1")
	for line in inf1:
		mobj=re.search("(.*?)\t.*?\t(.*?)\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",line)
		if mobj:
			tag_r=mobj.group(1)
			tag_c=mobj.group(2)
			pos_m=int(mobj.group(3))
			cigar=mobj.group(4)
			seq_r=mobj.group(5)
			mobj1=re.match("(\d+)S\d+M$",cigar)
			mobj2=re.match("\d+M(\d+)S$",cigar)
			mobj3=re.match("(\d+)M(\d+)D\d+M$",cigar)
			mobj4=re.match("(\d+)M(\d+)I\d+M$",cigar)
			mobj5=re.match("(\d+)M$",cigar)
			if mobj3:#long ref
				m_size=int(mobj3.group(1))
				ins_size=int(mobj3.group(2))
				if ins_size>=15:
					alter_seq=dict_chr[tag_c][pos_m+m_size-1:pos_m+m_size-1+ins_size]
					start=pos_m+m_size
					end=pos_m+m_size+ins_size-1
					print ("deletion_S_S\t%s\t%s\t%d\t%s\t%s\t%s" % (tag_c,tag_r,pos_m,cigar,seq_r,alter_seq),file=out1e)
			if mobj4:#long read
				m_size=int(mobj4.group(1))
				ins_size=int(mobj4.group(2))
				if ins_size>=15:
					insert_point=pos_m+m_size-1
					alter_seq=seq_r[m_size:m_size+ins_size]
					print ("insertion_S_S\t%s\t%s\t%d\t%s\t%s\t%s" % (tag_c,tag_r,pos_m,cigar,seq_r,alter_seq),file=out1e)
			if mobj1:
				n1=int(mobj1.group(1))
				if n1>20:
					addword2dict(dict_rec_d,tag_c,tag_r,line)
			if mobj2:
				n2=int(mobj2.group(1))
				if n2>20:
					addword2dict(dict_rec_d,tag_c,tag_r,line)
			if mobj5:#total coverage count
				n3=int(mobj5.group(1))
				for i in range(pos_m,pos_m+n3):
					addnum2dict(dict_rec_m,tag_c,i)
###################################################################
	del inf1
	print("Data structure initiation for large deletions")
	for key1 in dict_rec_d:
		list_keys=list(dict_rec_d[key1].keys())
		for key2 in list_keys:
			if len(dict_rec_d[key1][key2])==1:
				mobj=re.search(".*?\t.*?\t.*?\t.*?\t.*?\t(.*?)\t.*?\t.*?\t.*?\t.*?\t",dict_rec_d[key1][key2][0])
				if mobj:
					cigar=mobj.group(1)
					mobj1=re.match("(\d+)S(\d+)M$",cigar)
					mobj2=re.match("(\d+)M(\d+)S$",cigar)
					if mobj1:
						addword2dict2(dict_rec_r,key1,key2,dict_rec_d[key1][key2][0])
					if mobj2:
						addword2dict2(dict_rec_l,key1,key2,dict_rec_d[key1][key2][0])
				dict_rec_d[key1].pop(key2)
###################################################################
	print("Searching for large deletions")
	for key1 in dict_rec_d:
		for key2 in dict_rec_d[key1]:
			left2_list_tag=list()
			right2_list_tag=list()
			left2_list_pos=list()
			right2_list_pos=list()
			left2_list_M=list()
			right2_list_M=list()
			left2_list_S=list()
			right2_list_S=list()
			left2_list_seq=list()
			right2_list_seq=list()
			left2_list_cigar=list()
			right2_list_cigar=list()
			#Formation of data structure (reads that have more than 2 matches)
			for vk in dict_rec_d[key1][key2]:
				mobj=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",vk)
				if mobj:
					pos_m=int(mobj.group(1))
					cigar=mobj.group(2)
					seq_r=mobj.group(3)
					mobj1=re.match("(\d+)S(\d+)M$",cigar)
					mobj2=re.match("(\d+)M(\d+)S$",cigar)
					if mobj1:
						SM_S=int(mobj1.group(1))
						SM_M=int(mobj1.group(2))
						right2_list_tag.append(key2)
						right2_list_pos.append(pos_m)
						right2_list_seq.append(seq_r)
						right2_list_M.append(SM_M)
						right2_list_S.append(SM_S)
						right2_list_cigar.append(cigar)
					if mobj2:
						MS_M=int(mobj2.group(1))
						MS_S=int(mobj2.group(2))
						left2_list_tag.append(key2)
						left2_list_pos.append(pos_m)
						left2_list_seq.append(seq_r)
						left2_list_M.append(MS_M)
						left2_list_S.append(MS_S)
						left2_list_cigar.append(cigar)
			#Formation finished
			for i in range(len(left2_list_tag)):
				for j in range(len(right2_list_tag)):
					dis_lr=left2_list_pos[i]+left2_list_M[i]-right2_list_pos[j]
					dis_2=left2_list_M[i]-right2_list_S[j]
					length_read_l=left2_list_M[i]+left2_list_S[i]
					if dis_lr>=0 and dis_lr<=10:#long read
						alter_seq=seq_r[left2_list_M[i]:length_read_l-right2_list_M[j]+dis_lr]
						print ("insertion_S_T\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d" % (key1,left2_list_tag[i],left2_list_pos[i],right2_list_pos[j],left2_list_cigar[i],right2_list_cigar[j],left2_list_seq[i],alter_seq,dis_lr),file=out1e)
					if dis_lr<0 and dis_2<=10 and dis_2>=0:#long ref
						alter_seq=dict_chr[key1][left2_list_pos[i]+left2_list_M[i]-1:right2_list_pos[j]-1+dis_2]
						print ("deletion_S_T\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d" % (key1,left2_list_tag[i],left2_list_pos[i],right2_list_pos[j],left2_list_cigar[i],right2_list_cigar[j],left2_list_seq[i],alter_seq,dis_2),file=out1e)
	del dict_rec_d
	###############################################################################
	print("Searching for large insertions")
	for key1 in dict_rec_l:
		if key1 in dict_rec_r:
			sorted_list_r=sorted(dict_rec_r[key1].keys(),key=lambda x:r_m_l(dict_rec_r[key1][x]))
			sorted_list_l=sorted(dict_rec_l[key1].keys(),key=lambda x:r_m_l(dict_rec_l[key1][x]))
			for key2_1 in sorted_list_l:
				mobj_l=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",dict_rec_l[key1][key2_1])
				pos_m_l=int(mobj_l.group(1))  
				cigar_l=mobj_l.group(2)
				seq_r_l=mobj_l.group(3)
				mobj2=re.match("(\d+)M(\d+)S$",cigar_l)
				if mobj2:
					MS_M=int(mobj2.group(1))
					MS_S=int(mobj2.group(2))
					for key2_2 in sorted_list_r:
						if key2_2 in dict_rec_r[key1]:
							mobj_r=re.search(".*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t.*?\t(.*?)\t",dict_rec_r[key1][key2_2])
							pos_m_r=int(mobj_r.group(1)) 
							cigar_r=mobj_r.group(2)
							seq_r_r=mobj_r.group(3)
							mobj3=re.match("(\d+)S(\d+)M$",cigar_r)
							if mobj3:
								SM_M=int(mobj3.group(1))
								SM_S=int(mobj3.group(2))
								length_read_l=MS_M+MS_S
								length_read_r=SM_M+SM_S
								dis_lr1=pos_m_l+MS_M-pos_m_r
								if dis_lr1<(-100):
									break
								if dis_lr1>100:
									dict_rec_r[key1].pop(key2_2)
								if dis_lr1>=0 and dis_lr1<=10:
									alter_seq=find_overlap(seq_r_l,seq_r_r,pos_m_l,pos_m_r,MS_M,SM_S,dis_lr1,length_read_l,length_read_r)
									if not alter_seq=='0':
										print ("insertion_M\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d" % (key1,key2_1,key2_2,pos_m_l,pos_m_r,cigar_l,cigar_r,seq_r_l,seq_r_r,alter_seq,dis_lr1),file=out1e)
										dict_rec_r[key1].pop(key2_2)
										break
									else:
										if MS_M/length_read_l>=0.5 and SM_M/length_read_r>=0.5:
											print ("%s\t%d\t%d\t%s\t%s"%(key1,pos_m_l+MS_M,pos_m_r,key2_1,key2_2),file=out2e)
	print("Releasing RAM")
	del dict_rec_l
	del dict_rec_r
	t_elapsed = (time.clock() - t_start)
	print("Time used:",t_elapsed)
	out1e.close()
	out2e.close()
	inf1=open(out1_file)#completed reads
	inf2=open(out2_file)#uncompleted reads
	uncom_num_dict=dict()
	for line in inf2:
		mobj1=re.match("(.*?)\t(\d+)\t(\d+)\t.*?\t.*",line)
		if mobj1:
			oritag=mobj1.group(1)
			start1=int(mobj1.group(2))
			start2=int(mobj1.group(3))
			this_key=str(start1)+'_'+str(start2)
			addnum2dict(uncom_num_dict,oritag,this_key)
	num_uncom=0
	print("tag\tchromosome\tposition\tsplicing_depth\tnon-splicing_depth",file=out_uncom)
	for key1 in uncom_num_dict:
		mobj_key=re.match("(.*)_(\d+)_(\d+)_(\d+)",key1)
		ori_key=mobj_key.group(1)
		num_key=mobj_key.group(2)
		start_key=int(mobj_key.group(3))
		end_key=int(mobj_key.group(4))
		for key2 in uncom_num_dict[key1]:
			if uncom_num_dict[key1][key2]>=5:
				num_uncom+=1
				mobj1=re.match("(\d+)_(\d+)",key2)
				pos1=int(mobj1.group(1))
				sumi=0
				numi=0
				meani=0
				if key1 in dict_rec_m:
					for i in range(pos1-5,pos1+5):
						if i in dict_rec_m[key1]:
							sumi+=dict_rec_m[key1][i]
							numi+=1
						else:
							numi+=1
					meani=sumi/numi
				else:
					meani=0
				#print("I%d\t%s\t%s\t%s\t%d\t%d\t%3f" % (num_ins,key1,mobj1.group(1),ins_dict[key1][key2],len(ins_dict[key1][key2]),ins_num_dict[key1][key2],meani),file=out_ins)
				print("I_U%d\t%s\t%d\t%d\t%3f" % (num_uncom,ori_key,int(mobj1.group(1))+start_key-1,uncom_num_dict[key1][key2],meani),file=out_uncom)


#############public region###############
ins_dict=dict()
ins_num_dict=dict()
del_dict=dict()
del_num_dict=dict()
print("Building hash-like data structure of output data")
#build hash-like data structure of output data
for line in inf1:
	mobj1=re.search("insertion_S_S",line)
	mobj2=re.search("deletion_S_S",line)
	mobj3=re.search("insertion_S_T",line)
	mobj4=re.search("deletion_S_T",line)
	mobj5=re.search("insertion_M",line)
	if mobj1:
		m1=re.match(".*?\t(.*?)\t.*?\t(.*?)\t(\d+)M(\d+)I\d+M\t.*?\t(.*)",line)
		tag_c=m1.group(1)
		pos_m=m1.group(2)
		m_size=m1.group(3)
		ins_size=m1.group(4)
		seq=m1.group(5)
		insert_point=int(pos_m)+int(m_size)-1
		this_key=str(insert_point)+'_'+ins_size
		addword2dict2(ins_dict,tag_c,this_key,seq)
		addnum2dict(ins_num_dict,tag_c,this_key)
	elif mobj2:
		m1=re.match(".*?\t(.*?)\t.*?\t(.*?)\t(\d+)M(\d+)D\d+M\t.*?\t(.*)",line)
		tag_c=m1.group(1)
		pos_m=m1.group(2)
		m_size=m1.group(3)
		ins_size=m1.group(4)
		seq=m1.group(5)
		start=int(pos_m)+int(m_size)
		end=int(pos_m)+int(m_size)+int(ins_size)-1
		this_key=str(start)+'_'+str(end)
		addword2dict2(del_dict,tag_c,this_key,seq)
		addnum2dict(del_num_dict,tag_c,this_key)
	elif mobj3:
		m1=re.match(".*?\t(.*?)\t.*?\t(.*?)\t(.*?)\t(\d+)M\d+S\t\d+S(\d+)M\t.*?\t(.*?)\t.*",line)
		tag_c=m1.group(1)
		l_pos=m1.group(2)
		r_pos=m1.group(3)
		lm=m1.group(4)
		rm=m1.group(5)
		seq=m1.group(6)
		insert_point=int(l_pos)+int(lm)-1
		ins_size=len(seq)
		this_key=str(insert_point)+'_'+str(ins_size)
		addword2dict2(ins_dict,tag_c,this_key,seq)
		addnum2dict(ins_num_dict,tag_c,this_key)
	elif mobj4:
		m1=re.match(".*?\t(.*?)\t.*?\t(.*?)\t(.*?)\t(\d+)M\d+S\t\d+S(\d+)M\t.*?\t(.*?)\t(.*)",line)
		tag_c=m1.group(1)
		l_pos=m1.group(2)
		r_pos=m1.group(3)
		lm=m1.group(4)
		rm=m1.group(5)
		seq=m1.group(6)
		start=int(l_pos)+int(lm)
		ins_size=len(seq)
		end=start+ins_size-1
		this_key=str(start)+'_'+str(end)
		addword2dict2(del_dict,tag_c,this_key,seq)
		addnum2dict(del_num_dict,tag_c,this_key)
	elif mobj5:
		m1=re.match(".*?\t(.*?)\t.*?\t.*?\t(.*?)\t(.*?)\t(\d+)M\d+S\t\d+S(\d+)M\t.*?\t.*?\t(.*?)\t.*",line)
		tag_c=m1.group(1)
		l_pos=m1.group(2)
		r_pos=m1.group(3)
		lm=m1.group(4)
		rm=m1.group(5)
		seq=m1.group(6)
		insert_point=int(l_pos)+int(lm)-1
		ins_size=len(seq)
		this_key=str(insert_point)+'_'+str(ins_size)
		addword2dict2(ins_dict,tag_c,this_key,seq)
		addnum2dict(ins_num_dict,tag_c,this_key)

num_ins=0
num_del=0
print("tag\tchromosome\tposition\tsequence\tlength\tsplicing_depth\tnon-splicing_depth",file=out_ins)
for key1 in ins_dict:
	mobj_key=re.match("(.*)_(\d+)_(\d+)_(\d+)",key1)
	ori_key=mobj_key.group(1)
	num_key=mobj_key.group(2)
	start_key=int(mobj_key.group(3))
	end_key=int(mobj_key.group(4))
	for key2 in ins_dict[key1]:
		if ins_dict[key1][key2]!='' and ins_num_dict[key1][key2]>=2:
			num_ins+=1
			mobj1=re.match("(\d+)_(\d+)",key2)
			pos1=int(mobj1.group(1))
			sumi=0
			numi=0
			meani=0
			if key1 in dict_rec_m:
				for i in range(pos1-5,pos1+5):
					if i in dict_rec_m[key1]:
						sumi+=dict_rec_m[key1][i]
						numi+=1
					else:
						numi+=1
				meani=sumi/numi
			else:
				meani=0
			#print("I%d\t%s\t%s\t%s\t%d\t%d\t%3f" % (num_ins,key1,mobj1.group(1),ins_dict[key1][key2],len(ins_dict[key1][key2]),ins_num_dict[key1][key2],meani),file=out_ins)
			print("I%d\t%s\t%d\t%s\t%d\t%d\t%3f" % (num_ins,ori_key,int(mobj1.group(1))+start_key-1,ins_dict[key1][key2],len(ins_dict[key1][key2]),ins_num_dict[key1][key2],meani),file=out_ins)
print("tag\tchromosome\tstart\tend\tsequence\tlength\tsplicing_depth\tnon-splicing_depth",file=out_del)
for key1 in del_dict:
	mobj_key=re.match("(.*)_(\d+)_(\d+)_(\d+)",key1)
	ori_key=mobj_key.group(1)
	num_key=mobj_key.group(2)
	start_key=int(mobj_key.group(3))
	end_key=int(mobj_key.group(4))
	for key2 in del_dict[key1]:
		if del_dict[key1][key2]!='' and del_num_dict[key1][key2]>=2:
			num_del+=1
			mobj2=re.match("(\d+)_(\d+)",key2)
			pos1=int(mobj2.group(1))
			sumi=0
			numi=0
			meani=0
			if key1 in dict_rec_m:
				for i in range(pos1-5,pos1+5):
					if i in dict_rec_m[key1]:
						sumi+=dict_rec_m[key1][i]
						numi+=1
					else:
						numi+=1
				meani=sumi/numi
			else:
				meani=0
			#print("D%d\t%s\t%s\t%s\t%s\t%d\t%d\t%3f" % (num_del,key1,mobj2.group(1),mobj2.group(2),del_dict[key1][key2],len(del_dict[key1][key2]),del_num_dict[key1][key2],meani),file=out_del)
			print("D%d\t%s\t%d\t%d\t%s\t%d\t%d\t%3f" % (num_del,ori_key,int(mobj2.group(1))+start_key-1,int(mobj2.group(2))+start_key-1,del_dict[key1][key2],len(del_dict[key1][key2]),del_num_dict[key1][key2],meani),file=out_del)


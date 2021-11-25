import os
import subprocess
import sys
from subprocess import call
from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable
#table = CodonTable.unambiguous_dna_by_name["Standard"]

args = sys.argv
#args = ['', 'aRNH_and_approximates_Afil_genome.fa', 'Afil.step3.rpsbproc.result','Afil']

db = dict()
for record in SeqIO.parse(args[1],'fasta'):
	db[record.id] = record

blast = open(args[2],'r').readlines()

d = {'RNase_H_like':'gRH', 'Retrotrans_gag':'GAG', 'RT_like':'gRT', 'rve':'INT', 'RNase_HI_like':'aRH', 'RNase_HI_RT_Ty3':'gRH', 'RT_LTR':'gRT','pepsin_retropepsin_like':'PRo','retropepsin_like':'PRo'}
d_length = {'gRH' : 121, 'GAG' : 97, 'gRT' : 177, 'INT' : 116, 'aRH' : 128, 'PRo' : 92}


with open('structure_information','w') as si:
	toWrite = 'ID\t\tStructure_a\t\tStructure_b\t\tStructure_c\n'
	si.write(toWrite) 
    
    
def inf_domain(s):
	s = s.split()
	coords = [int(s[4]),int(s[5])]
	name = d[s[9]]
	frame = int(s[1].split('[')[1].split(']')[0])
	'''
	if (0.85 * d_length[name]) <= float(coords[1]-coords[0]):
		condition = True
	else:
		condition = False
	'''
	if  (float(s[6])> 0.001) or (frame < 0) : #(coords[0] > coords[1])  or
		return []
	else:
		return [coords[0],coords[1],name,frame] #,condition]

def mine_domain(l,mass):
	'''
	if mass[2] == (1 or 2):
		fix = -1
	else:
		fix = 2
	'''
	translation = str(Seq.Seq(str(l[0].seq)[(mass[0]-1):(mass[1]-1)],generic_dna).translate()) #(table)
	text = '>' + str(l[0].id.split('_domain')[0]) + '_\"' + mass[2] + '\"-domain' + '|' + str(mass[0]) + '-' + str(mass[1]) + '\n' + translation + '\n'
	return text

qwerty = dict()
def sorting_elements(l,mass):
	global qwerty

	def find_tether_domain(lenght,mass,start,direction):
		if direction < 0:
			ma = start
			mi = start + 2500*direction
			if mi < 0:
				mi = 0 
		else:
			mi = start
			ma = start + 2500*direction
			if ma > lenght:
				ma = lenght
		a = []
		for i in mass:
			if (i[0] in range(mi,ma)) and (direction > 0):
				a.append(i)
			elif (i[1] in range(mi,ma)) and (direction < 0):
				a.append(i)
		if len(a) == 1:
			return a[0]
		else:
			mi = lenght
			ma = [] # now it is an element from mass
			if direction > 0:
				for i in a:
					if (i[1] - start < mi):
						mi = i[1] - start
						ma = i
			else:
				for i in a:
					if (start - i[0] < mi):
						mi = start - i[0]
						ma = i
			return ma

	def try_hard(lenght,mass,start,direction):
    		c = find_tether_domain(lenght,mass,start,direction)
    		if c == []:
        		return [['','','',''],'']
    		le = 1000000
    		if direction < 0:
        		le = start - c[1]
    		else:
        		le = c[0] - start
    		if le < 1300:
        		b = '.'
    		#elif le < 1600:
        	#	b = ';'
    		else:
        		b = ','
    		return(c,b)

	global is_good

	#compound = ''
	aRH = []
	top_gRT = []
	for i in mass:
		if 'aRH' in i:
			aRH.append(i)
	#print(aRH,len(aRH))
	if len(aRH) == 0:
		print('[ERROR] We`ve got loose an aRH domain!')
		return
	elif len(aRH) == 1:
		aRH = aRH[0] ### replacing aRH = [[c1,c2,name,fr]] to aRH = [c1,c2,name,fr]
	else:
		delta = 20000
		top_aRH = []
		for i in aRH:
			c = i[0] - (l[1]/2)
			if ((c > 0) and (c < delta)):
				delta = i[0] - (l[1]/2)
				top_aRH = i
			elif ((c < 0) and ((-c) < delta)):
				delta = (l[1]/2)-i[0]
				top_aRH = i
			else:
				continue
		aRH = top_aRH
	#is_good += ((aRH[1] - aRH[0]+1)/3)
	#print('aRH has ' + str((aRH[1] - aRH[0]+1)/3) + 'length!')

	top_len = 20000
	for i in mass:
		if 'gRT' in i:
			#print(i)
			c = aRH[0] - i[1]
			lenght = max(c,(-c))
			#print(c)
			if (lenght < 4500) and (lenght < top_len):
				#print(aRH,c, i,end='\t')
				#print('\n')
				top_len = lenght
				top_gRT = i
				#print(lenght,top_gRT)

	if top_gRT != []:
		a = aRH[0]
		left = aRH #left is a startpoint here
	else:
		print('\nThe element has no RT domain!')
		return

	structure = ''
	old_structure = ''
	# ADDING LEFT DOMAINS
	PRo_edge = '*'
	p = [left[0]-1,'','default']
	last_domain = ''
	domains_from_left = []
	while p != ['','','','']:
		new_p = try_hard(l[1],mass,p[0],-1)
		if len(p) == 3:
			structure=new_p[0][2]+'!'+structure
			last_domain = new_p[0][2]
			old_structure = new_p[0][2]+'.'+old_structure
			a = left[0]
		else:
			if new_p[0][2] == last_domain:
				structure=new_p[0][2]+new_p[1]+structure
				old_structure = '   .'+old_structure
			else:
				structure=new_p[0][2]+new_p[1]+structure
				old_structure = new_p[0][2]+'.'+old_structure
				last_domain = new_p[0][2]

			if new_p[0][2] == 'GAG':
				PRo_edge = a
			else:
				PRo_edge = '*'
			a = new_p[0][0]			
            
		p = new_p[0]
		if p[0] != '':
			short_len = (p[1] - p[0] +1) /3
			if short_len <= d_length[p[2]]:
				is_good += short_len
			else:
				is_good += (2 * d_length[p[2]] - short_len)
			domains_from_left.append(p)

	structure = structure[:-1]
	old_structure = old_structure[:-1]
	#edge0 = a 
	# ADDING RIGHT DOMAINS
	p = ['',left[0]-1,'default']
	domains_from_right = []
	last_domain = structure[-3:]
	#print(p)
	while (p != ['','','','']):
		if len(p) == 3:
			p = try_hard(15000,mass,p[1],1)
			structure=structure+'!'+p[0][2]
			last_domain = p[0][2]
			old_structure = old_structure +'.'+p[0][2]
			#a = left[1]
		else:
			p = try_hard(15000,mass,p[1],1)
			#print(p,old_structure)
			if last_domain == p[0][2]:
				structure=structure+p[1]+p[0][2]
				old_structure = old_structure + '.   '
			else:
				structure=structure+p[1]+p[0][2]
				old_structure = old_structure +'.'+p[0][2]
				last_domain = p[0][2]
			#a = p[0][1]
		p = p[0]
		if p[0] != '':
			short_len = (p[1] - p[0] +1) /3
			if short_len <= d_length[p[2]]:
				is_good += short_len
			else:
				is_good += (2 * d_length[p[2]] - short_len)
			domains_from_right.append(p)
	old_structure = old_structure[1:-1]
	#edge1 = a
	z2 = list()
	for i in structure:
		if i in ';.,!':
			z2.append(i)

	for i in range(len(z2)):
		if z2[i] == '!':
			center = i
	right = 0
	for i in range(center+1,len(z2)-1):
		if z2[i] == ',':
			if z2[i+1] in '.;':
				right = i
				break
			else:
				if old_structure.split('.')[i+1] == 'aRH':
					right = i+1
				else:
					right = i
				break
	if right == 0:
		right = len(z2)+1

	left = len(z2) # left became a len of z2 - massive of '!.,;' symbols
	i = center
	while i != 0:
		if z2[i] == ',':
			if z2[i-1] in '.;':
				left = i # now left should became a startpoint of final structure in z2 coordinates
				break
			else:
				if old_structure.split('.')[i-1] == 'aRH':
					left = i-1
				else:
					left = i # now left should became a startpoint of final structure in z2 coordinates
				break
		i -= 1

	# make new realistic edges great again
	
	if domains_from_left == []:
		edge0 = aRH[0]
	elif (i!=0):
		edge0 = domains_from_left[center-i-1][0]
	elif i == 0:
		edge0 = domains_from_left[center-i][0]
	if domains_from_right == []:
		edge1 = aRH[1]
	elif (right != len(z2)+1):
		edge1 = domains_from_right[right-center-1][1]
	elif right == len(z2)+1:
		edge1 = domains_from_right[right-center-2][1]
	#print(edge0,edge1)
    
	if (i == 0) and (left == len(z2)):
		left = -1 # left in z2 coordinates in the case when element startswith first found domain

	if left == -1:
		left = 0	 # left became a startpoint in old_structure coordinates
	else:
		left = 4+left*4  # left became a startpoint in old_structure coordinates

	if right > len(z2):
		right = len(structure)	# right became a startpoint in old_structure coordinates
	else:
		right = 4 + right*4-1   # right became a startpoint in old_structure coordinates

	end_structure = old_structure[left:right]
	end_structure = end_structure.replace(' ','')
	while end_structure.find('..') != -1:
		end_structure = end_structure.replace('..','.')
	if end_structure.startswith('.'):
		end_structure = end_structure[1:]
	if end_structure.endswith('.'):
		end_structure = end_structure[:-1]

	with open('structure_information','a') as si:
		toWrite = l[0].id + '\t\t' + structure + '\t\t' +  old_structure + '\t\t' + end_structure + '\n'
		si.write(toWrite)

	# NEW TEST OF FRAME
	one_frame = str(Seq.Seq(str(l[0].seq)[(edge0-1):(edge1)],generic_dna).translate())
	if PRo_edge != '*':
		one_frame_wo_GAG = str(Seq.Seq(str(l[0].seq)[(PRo_edge-1):(edge1)],generic_dna).translate())
	#print(one_frame, len(one_frame), (edge1-edge0+1)/3.0)
	if ('*' not in one_frame) and (len(one_frame)>200):
		is_good += 500.0
		#print('no Stop-codons!')
	elif (PRo_edge != '*') and ('*' not in one_frame_wo_GAG) and (len(one_frame)>200):
		is_good += 500.0
		#print('no Stop-codons!')
	

	# BLOCK OF FINDING LTRs
	with open('find_LTR.fa','w') as blastf:
		toWrite = '>' + str(l[0].id) + '\n' + str(l[0].seq)
		blastf.write(toWrite)
	#if '-trf' in args:
	#	call('/home/jack/Desktop/Mike/scripts/last_element_mining_algorithm/trf409_2.linux64 find_LTR.fa 2 5 7 80 20 30 500 -h -m', shell = 'True')
	#	call('makeblastdb -in find_LTR.fa.2.5.7.80.20.30.500.mask -dbtype nucl -out findLTRs_DB', shell = 'True', stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	#	call('blastn -query find_LTR.fa.2.5.7.80.20.30.500.mask -db findLTRs_DB -out blastn_LTR.out -outfmt 6 -evalue 1e-3 -num_threads 4 ', shell = 'True') #-max_hsps 5
	#else:
	call('makeblastdb -in find_LTR.fa -dbtype nucl -out findLTRs_DB', shell = 'True', stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	call('blastn -query find_LTR.fa -db findLTRs_DB -out blastn_LTR.out -outfmt 6 -evalue 1e-3 -num_threads 4 ', shell = 'True') #-max_hsps 5
	
	with open('blastn_LTR.out','r') as blastf:
		r = blastf.readlines()
	#print(r)
	LTR_char = []
	new_edge0,new_edge1 = -1, -1
	delta_l, delta_r  = 20000, 20000
	for i in r:
		for j in r:
			if i != j:
				if float(i.split('\t')[11].split('\n')[0]) == float(j.split('\t')[11].split('\n')[0]) and (int(i.split('\t')[6]) <int(i.split('\t')[7])) and (int(i.split('\t')[8])< int(i.split('\t')[9])) and (int(j.split('\t')[3]) > 100):
					mi = min(int(j.split('\t')[6]),int(j.split('\t')[8]))
					ma = max(int(j.split('\t')[7]),int(j.split('\t')[9]))
					if (mi < edge0) and (ma > edge1) and (edge0 - mi < delta_l) and (ma - edge1 < delta_r): # and ('tg' in str(l[0].seq)[mi:mi+10].lower()) and ('ca' in str(l[0].seq)[ma-2:ma].lower()): 
						delta_r = ma - edge1
						delta_l = edge0 - mi
						new_edge0 = mi
						new_edge1 = ma
						LTR_char = [float(j.split('\t')[2]),int(j.split('\t')[3])] # %identity, legnth of LTRs

	#global non_LTR
	#non_LTR = '>' + str(l[0].id.split('ID')[0]) + 'ID' + str(l[0].id.split('ID')[1].split('_')[0]) + '_domain-structure:' + structure + '\n' + str(Seq.Seq(str(l[0].seq)[edge0:edge1],generic_dna).translate()) +'\n'

	if (new_edge1 != -1) and (new_edge0 != -1):
		edge0 = new_edge0
		edge1 = new_edge1


	if LTR_char == []:
		LTR_char = ''
	else:
		is_good += float(LTR_char[0])
		#print(str(LTR_char[0]) + ' is an LTRs identity!')
		LTR_char = '|LTR%id' + str(LTR_char[0]) + '-len' + str(LTR_char[1]) 

	'''
	m1 = [is_GAG,is_PR,is_gRH,is_INT]
	m2 = [test_condition,test_frame]

	for i in m1:
		if i == True:
			is_good += 1
	for i in m2:
		if i == True:
			is_good += 2
	#print(m1,m2,is_good)
	'''

	p = '>' + str(l[0].id.split('ID')[0]) + 'ID' + str(l[0].id.split('ID')[1].split('_')[0]) + '_domain-structure:' + end_structure + LTR_char + '|Score=' + str(is_good) + '\n' + str(l[0].seq)[edge0:edge1]+ '\n'
	a = []
	a.append((edge1-edge0))
	a.append(edge0)
	a.append(edge1)
	a.append(structure)
	qwerty[l[0].id.split('_domain')[0]] = a
	
	#print('TEST OF CONDITION IS',test_condition,',TEST OF FRAME IS',test_frame, ',GAG', is_GAG, ',PR', is_PR, ',gRH', is_gRH, ',INT', is_INT,end = '\t')
	#print()
	#print(p.split('\n')[0])

	return p
	#9/0



with open('%s_little_damaged_elements' % args[3],'w') as damaged_elems:
	with open('%s_bad_elements' % args[3],'w') as bad_elems:
		with open('%s_best_elements' % args[3],'w') as best_elems:
			with open('prot_domains.fa','w') as f:
				with open('%s_elements' % args[3],'w') as all_elements:
					table = []
					str_counter = 0
					for string in blast:
						str_counter += 1
						if (string[0] in '0123456789Q') or (str_counter == len(blast)):
							if string[0] == 'Q':
								# func of sorting previous set of domains and consturcting element
								if len(table) !=0:
									is_good = 0.0
									#non_LTR = ''
									good_element = sorting_elements(element,table)
									if good_element != None:
										all_elements.write(good_element)
										if is_good >= 1150:
											best_elems.write(good_element)
										elif is_good >=900:
											damaged_elems.write(good_element)
										else:
											bad_elems.write(good_element)
									#print(is_good)
									#print('\n\n')
									#9/0
									'''
									if good_element != None:
										for to_mining in table:
											#print(qwerty)
											#print(qwerty[element[0].id])
											#print(qwerty[element[0].id][1],qwerty[element[0].id][2],to_mining[0],to_mining[1])
											ss = qwerty[element[0].id][1]
											ee = qwerty[element[0].id][2]
											#9/0
											if (to_mining[0] in range(ss,ee)) and (to_mining[1] in range(ss,ee)):
												f.write(mine_domain(element, to_mining))
									'''
								string = string.split()
								el_num = string[1].split('_')[1]
								element = [db[string[4]],int(string[3])]
								#print(string[4])
								table = []
								f.write('\n')
							elif (string[0] in '0123456789'):
								# func of scrining domains, mining domains, writing domains
								s = inf_domain(string)
								#print(s)
								if s == []:
									continue
								f.write(mine_domain(element, s))
								#print(m)
								table.append(s)
								#print(table)
							else:
								# func of sorting previous (and last) set of domains and consturcting element
								is_good = 0.0
								good_element = sorting_elements(element,table)
								if good_element != None:
									all_elements.write(good_element)
									if is_good >= 1150:
										best_elems.write(good_element)
									elif is_good >=900:
										damaged_elems.write(good_element)
									else:
										bad_elems.write(good_element)
								#print(is_good)
with open('qwerty_table','w') as qw:
	for i in qwerty:
		qw.write(str(i) + ' : ' + str(qwerty[i]) + '\n')
'''
min_qwerty = 20000
for i in qwerty:
	#print(qwerty[i])
	if qwerty[i][0] < min_qwerty:
		min_qwerty = qwerty[i][0]
print(min_qwerty)
'''

call('rm log.txt', shell = 'True')

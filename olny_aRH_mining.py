import sys
from subprocess import call
from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable
#table = CodonTable.unambiguous_dna_by_name["Standard"]

args = sys.argv
#args = ['', 'TAT-ATHILA-CHROMOVIRIDAE.fa_clusters', 'gyDBfrom.step3.rpsbproc.result','gyDBfrom', 'project_name']

db = dict()
for record in SeqIO.parse(args[1],'fasta'):
	db[record.id] = record

blast = open(args[2],'r').readlines()

d = {'RNase_H_like':'gRH', 'Retrotrans_gag':'GAG', 'RT_like':'gRT', 'rve':'INT', 'RNase_HI_like':'aRH', 'RNase_HI_RT_Ty3':'gRH', 'RT_LTR':'gRT','pepsin_retropepsin_like':'PRo','retropepsin_like':'PRo'}
d_length = {'gRH' : 121, 'GAG' : 97, 'gRT' : 177, 'INT' : 116, 'aRH' : 128, 'PRo' : 92}

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
	#translation = translation
	#print(translation)
	text = '>' + str(l[0].id) + '_\"' + mass[2] + '\"-domain' + '|' + str(mass[0]) + '-' + str(mass[1]) + '\n' + translation + '\n'
	return text

noRT = True
no_all = True
toWrite = ''
with open('%s_aRH_without_others.fa' %args[3],'w') as f:
	table = []
	str_counter = 0
	for string in blast:
		str_counter += 1
		if (string[0] in '0123456789Q'):
			if string[0] == 'Q':
				string = string.split()
				el_num = string[1].split('_')[1]
				element = [db[string[4]],int(string[3])]
				#print(string[4])
				#table = []
				if (noRT == True) or (no_all == True):
					f.write(toWrite)
				noRT = True
				no_all = True
				toWrite = ''
			else: #if (string[0] in '0123456789'):
				# func of scrining domains, mining domains, writing domains
				s = inf_domain(string)
				#print(s)
				if s == []:
					continue
				if s[2] == 'aRH':
					toWrite = mine_domain(element, s)
				elif s[2] == 'gRT':
					noRT = False
					no_all = False
				else:
					no_all = False


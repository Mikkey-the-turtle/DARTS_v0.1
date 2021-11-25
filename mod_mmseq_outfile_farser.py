import sys
from Bio import Seq
from Bio import SeqIO

args = sys.argv
#args = ['', 'clusterRes_cluster.tsv', 'gRT_Afil_elements.faa']

mmseq_outfile_str = open(args[1],'r').readlines()

db = dict()
for rec in SeqIO.parse(args[2],'fasta'):
	db[rec.id] = rec.seq 

def add_elem(l):
	return [l.split('=')[1],l]

def choose_delegate(t):
	global tt
	true_table = []
	global db
	if len(t) < 3:
		return ''
	t.sort(reverse=True)
	max = 0.0
	delegate = ''
	counter = 0
	for j in t:
		if (float(j[0]) > max) and (float(j[0]) < 1300) and (len(j[1].split(':')[1].split('|')[0]) < 24) and (j[1].count('gRT') == 1) and (j[1].count('INT') == 1) and (j[1].count('gRH') == 1) and (j[1].count('GAG') < 2) and (j[1].count('PRo') < 2):
			max = float(j[0])
			delegate = j[1]
	for j in t:
		if ((float(j[0]) < 300) and (':gRT|' in j[1])):
			continue
		elif (len(j[1].split(':')[1].split('|')[0]) > 24) or (j[1].count('gRT') >= 2) or (j[1].count('INT') >= 2) or (j[1].count('gRH') >= 2) or (j[1].count('GAG') >= 2) or (j[1].count('PRo') >= 2) or (j[1].count('gRH') == 0) or (j[1].count('INT') == 0):
			continue
		elif ((delegate.find('aRH') != -1) and (j[1].find('aRH') != -1)) or ((delegate.find('aRH') == -1) and (j[1].find('aRH') == -1)):
			counter += 1
			true_table.append(j[1])
	if (delegate == '') or (counter < 3):
		return ''
	else:
		for i in true_table:
			wrt = delegate + '\t' + i + '\n'
			tt.write(wrt)
		return ('>[' + str(counter) + ']' + delegate + '\n' + str(db[delegate]) + '\n')

with open('Natural_table','w') as tt:
	with open('mmseq_delegates_%s' %args[2],'w') as f:
		name, clust = '', ''
		table = []
		for i in mmseq_outfile_str:
			#print(i)
			if name != i.split('\t')[0]:
				if name == '':
					name = i.split('\t')[0]
					element = i.split('\t')[1].split('\n')[0]
					table.append(add_elem(element))
					continue
				toWrite = choose_delegate(table)
				f.write(toWrite)
				#print(toWrite)
				table = []
				name = i.split('\t')[0]
				element = i.split('\t')[1].split('\n')[0]
				#print(element)
				table.append(add_elem(element))
			else:
				element = i.split('\t')[1].split('\n')[0]
				#print(element)
				table.append(add_elem(element))

		# When no one string left
		toWrite = choose_delegate(table)
		#print(toWrite)
		f.write(toWrite)




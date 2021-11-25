from Bio import SeqIO
from Bio import Seq
import sys

args = sys.argv
#args = ['', 'prot_domains.fa', 'Afil_elements']

# gene set
d = ['gRH', 'GAG', 'gRT', 'INT', 'aRH', 'PRo']
0
# parse all prots from file of parsed rpsbproc domains-aa-seq
a = dict()
b = open('qwerty_table','r').readlines()
for i in b:
	a[i.split(' : ')[0]] = [int(i.split(',')[1]), int(i.split(', ')[2])]


domains = dict()
for record in SeqIO.parse(args[1],'fasta'):
	temp_id = str(record.id).split('ID')[0] + 'ID' + str(record.id).split('ID')[1].split('_')[0]
	if (temp_id not in a) or (int(record.id.split('|')[1].split('-')[0]) not in range(a[temp_id][0],a[temp_id][1])):
		continue
	elif ((temp_id + '|' + str(record.id).split('\"')[1]) not in domains):
		#print(temp_id+ '|' + str(record.id).split('\"')[1])
		domains[temp_id + '|' + str(record.id).split('\"')[1]] = record
	elif (len(record.seq) > len(domains[temp_id+ '|' + str(record.id).split('\"')[1]].seq)): # and (int(record.id.split('|')[1].split('-')[0]) in range(a[temp_id][0],a[temp_id][1])):
		# Need to concatinate ?
		#print(temp_id+ '|' + str(record.id).split('\"')[1])
		#print('New seq', str(record.id), str(record.seq), len(record.seq), int(record.id.split('|')[1].split('-')[0]), int(record.id.split('-')[2]), '\nExisting seq', str(domains[temp_id+ '|' + str(record.id).split('\"')[1]].id),str(domains[temp_id+ '|' + str(record.id).split('\"')[1]].seq), len(domains[temp_id+ '|' + str(record.id).split('\"')[1]].seq), 'elements coordinates', a[temp_id][0],a[temp_id][1])
		domains[temp_id + '|' + str(record.id).split('\"')[1]] = record

#print(a,len(a))

# parse file of all finded elements (needs only ids)
elements = []
for which_needs in SeqIO.parse(args[2],'fasta'):
	elements.append(str(which_needs.id))
	#elements.append(str(which_needs.id).split('_domain')[0])

#print(domains)
#print(elements)

def one_gene(de,dd,gene):
	short = dict()
	for i in de:
		short[(str(i).split('_domain')[0])] = i
		#print(short)
		#9/0
	with open('%s_%s.faa' % (gene,args[2]),'w') as f:
		for elem in short:
			if (elem+ '|' + gene) not in dd:
				print('There is no gene '+ gene + ' in ' + elem + ' !')
			else:
				#print(elem, short[elem],(dd[(elem+ '|' + gene)].seq),dd[(elem+ '|' + gene)].id,'\n')
				toWrite = '>' + str(short[elem]) +  '\n' + str(dd[(elem+ '|' + gene)].seq) + '\n'
				#print(toWrite)
				f.write(toWrite)
	#print(short)
for g in d:
	#print(g)
	one_gene(elements,domains,g)

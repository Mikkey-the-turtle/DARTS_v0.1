import sys
from subprocess import call
from Bio import Seq
from Bio import SeqIO

args = sys.argv
#args = ['', 'Afil', 'gRT_%s_elements_after_cdhit80_new_delegates.faa' OR 'mmseq_delegates_%s']


#gRT_Afil_elements_after_cdhit80_new_delegates.faa
db = dict()
for rec in SeqIO.parse(args[2], 'fasta'):
	db[rec.id] = rec

domains = ['gRH', 'GAG', 'INT', 'aRH', 'PRo']

def find_delegates_domains_one_type(base,dom):

	dombase = dict()
	for record in SeqIO.parse('%s_%s_elements.faa' % (dom,args[1]),'fasta'): #INT_Afil_elements.faa
		dombase[record.id] = record
	#print(dombase)

	with open('%s-domain_%s_elements_delegates.faa' % (dom,args[1]),'w') as f:
		for i in base:
			#print(i)
			#print('>' + str(dombase[i].id) + '\n' + str(dombase[i].seq) + '\n')
			#print(dombase[i].seq,len(dombase[i].seq))
			#print(dombase[i])
			if i in dombase:
				toWrite = '>' + str(dombase[i].id) + '\n' + str(dombase[i].seq) + '\n'
				f.write(toWrite)

for j in domains:
	print(j)
	find_delegates_domains_one_type(db,j)
	print()

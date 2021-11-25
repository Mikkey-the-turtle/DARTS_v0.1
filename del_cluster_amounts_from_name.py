import sys
from Bio import Seq
from Bio import SeqIO


args = sys.argv
#args = ['', 'fasta-file_without_sizes', 'fasta-file_with_sizes']

d = dict() #exististed gRT elements
for record in SeqIO.parse(args[1],'fasta'):
	d[record.id] = record
	#print(record)
	

#domains = ['aRH','gRH', 'GAG', 'INT', 'PRo']
#print(d)

a = dict()
for record in SeqIO.parse(args[2],'fasta'):
	if ']' in record.id:
		a[record.id.split(']')[1]] = record
	else:
		a[record.id] = record
with open('rewrited_%s' % (args[2]),'w') as f:
	for i in a:
		#print(str(d[i]))
		toWrite = '>' + str(d[i].id) + '\n' + str(a[i].seq) + '\n'
		f.write(toWrite)

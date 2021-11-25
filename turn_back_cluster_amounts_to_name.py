import sys
from Bio import SeqIO
from Bio import Seq

args = sys.argv
#args = ['', 'fasta-file_with_sizes', 'fasta-file_without_sizes']

d = dict() #exististed gRT elements
for record in SeqIO.parse(args[1],'fasta'):
	if ']' in record.id:
		d[record.id.split(']')[1]] = record
	else:
		d[record.id] = record
	#print(record)
	

#domains = ['aRH','gRH', 'GAG', 'INT', 'PRo']
#print(d)

a = dict()
for record in SeqIO.parse(args[2],'fasta'):
	a[record.id] = record
with open('rewrited_%s' % (args[2]),'w') as f:
	for i in a:
		#print(str(d[i]))
		toWrite = '>' + str(d[i].id) + '\n' + str(a[i].seq) + '\n'
		f.write(toWrite)

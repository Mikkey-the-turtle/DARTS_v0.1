import sys
from subprocess import call

args = sys.argv
#args = ['', 'file_to_split']

r = open(args[1],'r').read().split('>')

d = []
for i in r:
	d.append(len(i))

big_chromosome = max(d)

weight = 0
for i in d:
	weight += i

if weight > 5000000000:
	print('1') # need 'genome_splitter_script!'
elif big_chromosome > 100000000:
	print('2') # need 'sequence_splitter_script!'
else:
	print('3') # Use it w/o modifying!


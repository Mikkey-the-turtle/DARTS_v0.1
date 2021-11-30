import sys

args = sys.argv
#args = ['', 'input_file', 'output_file']

file = open(args[1],'r').read()

with open(args[2],'w') as write_f:
	for i in file:
		if i != ' ':
			write_f.write(i)
		else:
			write_f.write('_')
	
		

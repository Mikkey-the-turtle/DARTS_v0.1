import sys

args = sys.argv
#args = ['', 'Afil', 'gRT']

file = open('%s_%s_elements.faa' % (args[2],args[1]), 'r').read()

with open ('%s_%s_elements.faa' % (args[2],args[1]), 'w') as f:
	for i in file:
		if i == '*':
			f.write('X')
		else:
			f.write(i)

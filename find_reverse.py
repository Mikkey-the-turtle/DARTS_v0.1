import sys

args = sys.argv

b = open(args[1],'r').readlines()

for i in b:
	j = i.split('\t')
	if int(j[1]) > int(j[2]):
		print(i)

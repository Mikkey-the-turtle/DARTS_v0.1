import sys
from Bio import Seq
from Bio import SeqIO
#table = CodonTable.unambiguous_dna_by_name["Standard"]

args = sys.argv
#args = ['', 'Azolla_filiculoides.genome_v1.2.fa', 'Afil.step1.rpsbproc.result','Afil']

db = dict()
for record in SeqIO.parse(args[1],'fasta'):
	db[record.id] = record

blast = open(args[2],'r').readlines()

with open('Coords_1step','w') as bam:
		toBam = ''
		bam.write(toBam)

def new_elem(i):
	global db, args, l, is_left, is_right
	name = i[0]
	sc = db[i[0]]
	x = len(sc.seq)
	#s = 0
	#e = x
	#'''
	if is_left == True:
		s = i[1]
	elif i[1] - 7500 < 0:
		s = 0
	else:
		s = i[1] - 7500
	if is_right == True:
		e = i[2]
	elif i[2] + 7500 > x:
		e = x
	else:
		e = i[2] + 7500
	#'''
	if i[5] > 0:
		sc = '>' + args[3] + '_ID' + str(i[3]) +  '\n' + str(sc.seq[s:e]) + '\n'
		#print(sc)
		#translation = Seq.Seq(sc.split('\n')[1]).translate()
		#print(translation)
		#translation = str(Seq.Seq(str(l[0].seq)[(mass[0]-1):(mass[1]-1)],generic_dna).translate()) #(table)
	else:
		#print('Reverse!')
		sc = '>' + args[3] + '_ID' + str(i[3]) +  '\n' + str(Seq.Seq(str(sc.seq)[s+1:e+1]).reverse_complement()) + '\n'
		#print(sc.split('\n')[1])
		#print(Seq.Seq(sc.split('\n')[1]).translate())
		#9/0
	with open('Coords_1step','a') as bam:
		toBam = args[3] + '_ID' + str(i[3]) + '\t' + str(s) + '\t' + str(e) + '\t' + str(i[5]) + '\t' + str(i[0]) + '\n'
		bam.write(toBam)
	return sc

l = list()
ID = 1
for string in blast:
	if string[0] in '1234567890Q':
		string = string.split()
		#print(string,string[0])
		if string[0] == 'QUERY':
			scaff = string[4]
			#print(scaff)
		elif (string[0][0] in '1234567890') and (string[2] == 'Specific') and (float(string[6]) < 1e-3):
			start = int(string[4]) - 1
			end = int(string[5]) - 1
			name = string[9] + '_' + string[8]
			frame = int(string[1].split('[')[1].split(']')[0])
			l.append((scaff,start,end,ID,name,frame))
			ID += 1

l = sorted(l, key=lambda x: (x[0],x[1]))
exclude = []
with open('RT_and_approximates_%s_genome.fa' % args[3],'w') as f:
	for i in l:
		#print('\n')
		#rint(i,exclude)
		#if len(exclude) > 1:
		#	break
		if i in exclude:
			#print('lol!')
			continue
		else: 
			exclude.append(i)
		#if i in exclude:
		#	print('FUCK U, Misha!')
		mass = [i]
		for j in l:
			if (i[0] == j[0]) and ((-600 < i[1]-j[2] <600) or (-600 < i[2]-j[1] <600)) and (i[1] != j[1]) and (i[2] != j[2]):
				#print(i,j, exclude)
				mass.append(j)
				exclude.append(j)
		if len(mass) != 1:
			mi = min(x[1] for x in mass)
			ma = max(x[2] for x in mass)
			i = (i[0],mi,ma,i[3],i[4],i[5])
		mass = [i]
		for j in l:
			if (i[0] == j[0]) and ((-15000 < i[1]-j[2] <15000) or (-15000 < i[2]-j[1] <15000)) and (((i[5] > 0) and (j[5] >0)) or ((i[5] < 0) and (j[5] < 0))):
				mass.append(j)
		if len(mass) == 1:
			is_left = False
			is_right = False
			#toWrite = new_elem(i)
			#f.write(toWrite)
		else:
			mi = 0
			ma = len(db[i[0]].seq)
			for j in mass:
				if (j[1] < i[1]):
					temp = (i[1]+j[2])//2
					if temp > mi:
						mi = temp
				elif (j[1] > i[1]):
					temp = (i[2]+j[1])//2
					if temp < ma:
						ma = temp
			with open('lol','a') as lol:
				lol.write(str(i[0]) + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + str(i[3]) + '\n')
			if mi == 0:
				mi = i[1]
				is_left = False
			else:
				is_left = True
			if ma == len(db[i[0]].seq):
				ma = i[2]
				is_right = False
			else:
				is_right = True
			with open('lol','a') as lol:
				lol.write(str(i[0]) + '\t' + str(mi) + '\t' + str(ma) + '\t' + str(is_left) + '\t' + str(is_right) + '\n\n')
			i = (i[0],mi,ma,i[3],i[4],i[5])
		toWrite = new_elem(i)
		f.write(toWrite)
			#9/0




#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

args = sys.argv
#args = ['', 'genome_file', 'project_name']


db = dict()
for i in SeqIO.parse(args[1],'fasta'):
    db[i.id] = i.seq

counter = 0
i = 0
thr = 100000000
toWrite = ''

#z = open('blet','w')

for el in db:
    if len(db[el]) < thr - counter:
        toWrite += '>' + str(el) + '\n' + str(db[el]) + '\n'
        counter += len(db[el])
        #print(len(db[el]), counter)
    elif (len(db[el]) > thr - counter):
        c = len(db[el])
        #print(c)
        if c < (thr / 10):
            toWrite += '>' + str(el) + '\n' + str(db[el]) + '\n'
            with open('%s_%d.fa' %(args[2],i), 'w') as f:
                f.write(toWrite)
            counter = 0
            i += 1
            toWrite = ''
        else:
            while c > thr - counter:
                #print(c,thr,counter)
                toWrite += '>' + str(el) + '\n' + str(db[el][(len(db[el])-c):(len(db[el])-c+thr-counter)]) 
                with open('%s_%d.fa' %(args[2],i), 'w') as f:
                    f.write(toWrite)
                i += 1
                toWrite = ''
                c -= (thr - counter)
                counter = 0
            toWrite = '>' + str(el) + '\n' + str(db[el][(len(db[el])-c):]) + '\n'
            counter = c
        #print(c)
#    else:
#        print('lol, oshibo4ka!')
    #toWrite2 = toWrite + '\t' + str(i) + '\t\n'
    #z.write(toWrite2)
    #print(thr,counter, len(db[el]), el)

with open('%s_%d.fa' %(args[2],i), 'w') as f:
    f.write(toWrite)

#z.close()

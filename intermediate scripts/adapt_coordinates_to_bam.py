#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO

args = sys.argv
#args = ['', 'smansoni.fa_elements', 'Coords_1step.bam','qwerty_table']

db = dict()
for record in SeqIO.parse(args[1],'fasta'):
	db[record.id.split('_domain-structure:')[0]] = record.id

data1 = open(args[2],'r').read().split('\n')    
#data1 = ['Gbil_Ginkgo_biloba.scaf.fa_0.fa_ID10\t0\t12600\t+2\tcontig1', 'Gbil_Ginkgo_biloba.scaf.fa_0.fa_ID12\t45679\t61679\t+3' , 'Gbil_Ginkgo_biloba.scaf.fa_0.fa_ID17\t1012817\t1002619\t-1']

data2 = dict()
with open(args[3],'r') as qwerty_table:
    r = qwerty_table.readlines()
for i in range(len(r)):
    i = r[i]
    i = i.split(' : [')
    data2[i[0]] = [int(i[1].split(', ')[0]), int(i[1].split(', ')[1]), int(i[1].split(', ')[2]), str(i[1].split('\'')[1])]
#data2 = {'Gbil_Ginkgo_biloba.scaf.fa_0.fa_ID10' : [11513, 428, 11941, 'GAG.PRo,gRT.gRH!aRH.INT'], 'Gbil_Ginkgo_biloba.scaf.fa_0.fa_ID12' : [14425, 1105, 15530, 'GAG.PRo!aRH.gRT.gRH.INT'], 'Gbil_Ginkgo_biloba.scaf.fa_0.fa_ID17' : [10198, 2618, 12816, 'GAG.PRo!aRH.gRT.gRH.INT']}

#args = ['','smansoni.fa_elements','','']
with open('coordinates_table_of_%s' %args[1],'w') as f:
    for i in data1:
        #print('\n', i.split('\t')[0])

        parts = i.split('\t')
        if parts == ['']:
            break
        if i.split('\t')[0] not in db:
            continue
        elif '-len' in db[i.split('\t')[0]]:
            ltr_len = int(db[i.split('\t')[0]].split('-len')[1].split('|')[0])
        else:
            ltr_len = ''
            
        if (int(parts[-2]) > 0):
            #print(parts[1],parts[2],data2[parts[0]][1],data2[parts[0]][2])
            #print('start coord is : ' + str(int(parts[1])+data2[parts[0]][1]) + '\nend coord is : ' + str(int(parts[1])+data2[parts[0]][2]))
            contig_start = int(parts[1])+data2[parts[0]][1]
            contig_end = int(parts[1])+data2[parts[0]][2]
            if ltr_len != '':
                #print('Ltr_1',contig_start,contig_start+ltr_len,'Ltr_2',contig_end,contig_end-ltr_len)
# old           #toWrite = db[i.split('\t')[0]] + '\t' + parts[4] + '\t' + str(contig_start) + '\t' + str(contig_end) + '\t' + str(ltr_len) + '\t' + str(contig_start) + '-' + str(contig_start+ltr_len) + '\t' + str(contig_end) + '-' + str(contig_end-ltr_len) + '\n'
# old           #f.write(toWrite)
                toWrite = parts[4] + '\t' + str(contig_start) + '\t' + str(contig_start+ltr_len) + '\tLTR1\t' + db[i.split('\t')[0]] + '\t+\n'
                f.write(toWrite)
                toWrite = parts[4] + '\t' + str(contig_start+ltr_len) + '\t' + str(contig_end-ltr_len) + '\tinternal\t' + db[i.split('\t')[0]] + '\t+\n'
                f.write(toWrite)
                toWrite = parts[4] + '\t' + str(contig_end-ltr_len) + '\t' + str(contig_end) + '\tLTR2\t' + db[i.split('\t')[0]] + '\t+\n'
                f.write(toWrite)
            else:
                toWrite =  parts[4] + '\t' + str(contig_start) + '\t' + str(contig_end) + '\tinternal\t' + db[i.split('\t')[0]] + '\t+\n'
                f.write(toWrite)
        else:
            #print(parts[1],parts[2],data2[parts[0]][1],data2[parts[0]][2])
            #print('start coord is : ' + str(int(parts[2])-data2[parts[0]][1]) + '\nend coord is : ' + str(int(parts[2])-data2[parts[0]][2]))
            contig_start = int(parts[2])-data2[parts[0]][1]
            contig_end = int(parts[2])-data2[parts[0]][2]
            if ltr_len != '':
                #print('Ltr_1',contig_start,contig_start-ltr_len,'Ltr_2',contig_end,contig_end+ltr_len)
# old           #toWrite = db[i.split('\t')[0]] + '\t' + parts[4] + '\t' + str(contig_start) + '\t' + str(contig_end) + '\t' + str(ltr_len) + '\t' + str(contig_start) + '-' + str(contig_start-ltr_len) + '\t' + str(contig_end) + '-' + str(contig_end+ltr_len) + '\n'
# old           #f.write(toWrite)
                toWrite = parts[4] + '\t' + str(contig_end) + '\t' + str(contig_end+ltr_len) + '\tLTR2\t' + db[i.split('\t')[0]] + '\t-\n'
                f.write(toWrite)
                toWrite = parts[4] + '\t' + str(contig_end+ltr_len) + '\t' + str(contig_start-ltr_len) + '\tinternal\t' + db[i.split('\t')[0]] + '\t-\n'
                f.write(toWrite)
                toWrite = parts[4] + '\t' + str(contig_start-ltr_len) + '\t' + str(contig_start) + '\tLTR1\t' + db[i.split('\t')[0]] + '\t-\n'
                f.write(toWrite)
            else:
                toWrite = parts[4] + '\t' + str(contig_end) + '\t' + str(contig_start) + '\tinternal\t' + db[i.split('\t')[0]] + '\t-\n'
                f.write(toWrite)
        


#print('lol')
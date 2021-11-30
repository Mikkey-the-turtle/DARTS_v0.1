#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
args = sys.argv

counter = 0
thr = 1000000000
#thr = 2
i = 0
string = ''

with open(args[1],'r') as f:
#with open('test.fa','r') as f:
    for l in f:
        if l.startswith('>'):
            if counter == 0 :
                f1 = open('%s_%d.fa' %(args[1],i), 'w')
            if counter >= thr:
                f1.close()
                i +=1
                counter = 0
                f1 = open('%s_%d.fa' %(args[1],i), 'w')
            string += l
        else:
            string += l
            f1.write(string)
            counter += len(string)
            string = ''   

            
f1.close()
                    

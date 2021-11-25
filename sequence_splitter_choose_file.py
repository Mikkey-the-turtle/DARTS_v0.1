#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
args = sys.argv

counter = 0
thr = 100000000
#thr = 2
i = 0
string = ''
string_before = ''
str_b_counter = 0

with open(args[1],'r') as f:
#with open('test.fa','r') as f:
    for l in f:
        if l.startswith('>'):
            name = l
            i += 1
            f1 = open('%s_part_%d.fa' %(args[1],i), 'w')
            f1.write(name)
        elif counter >= thr:
            f1.close()
            i +=1
            counter = 0
            str_b_counter = 0
            f1 = open('%s_part_%d.fa' %(args[1],i), 'w')
            f1.write(name)
            f1.write(string_before)
            string += l
            f1.write(string)
            counter += len(string)
            string = ''
        else:
            string += l
            f1.write(string)
            counter += len(string)

            if str_b_counter < 1000:
                string_before += string
                str_b_counter += len(string)
            else:
                string_before = string_before[len(string):]
                string_before += string
            string = ''   

            
f1.close()
                    

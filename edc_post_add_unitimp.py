#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/01/24
@author: Naoto Hori
'''

import sys

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: SCRIPT [edc_pre.3i8h.out] [../post_edc_3i8h/edc_post.out] [post_edc.out]'
        print 'assuming at ~/ribo/ribo_015/forDunkleJenner/'
        sys.exit(2)
        
f_out = open(sys.argv[-1],'w')

class mpinfo:
    def __init__(self):
        self.ichain = 0
        self.ires = 0
        self.imp = 0
        self.ires_l = 0
        self.imp_l = 0
        self.res_name = ''
        self.mp_name = ''
        self.ires_aa = 0
        self.ires_u = 0
        self.imp_u = 0
    
cafemp2mpinfo = {}
for line in open(sys.argv[1],'r'):
    if line.find('#') != -1:
        continue
    linesp = line.split()
    mp = mpinfo()
    mp.ichain = int(linesp[0])
    mp.ires = int(linesp[1])
    mp.imp = int(linesp[2])
    mp.ires_l = int(linesp[3])
    mp.imp_l = int(linesp[4])
    mp.res_name = linesp[5].strip()
    mp.mp_name = linesp[6].strip()
    mp.ires_amber = int(linesp[7])
    mp.ires_u = int(linesp[8])
    mp.imp_u =  int(linesp[9])
    cafemp2mpinfo[mp.imp] = mp

for line in open(sys.argv[2],'r'):
    if line.find('#') != -1:
        f_out.write(line)
        continue
    if len(line) == 1:
        f_out.write('\n')
        continue
    linesp = line.split()
    ichain1 = int(linesp[0])
    ires1 = int(linesp[1])
    res1_name = linesp[2].strip()
    imp1 = int(linesp[3])
    mp1_name = linesp[4].strip()
    # linesp[5] is '|'
    ichain2 = int(linesp[6])
    ires2 = int(linesp[7])
    res2_name = linesp[8].strip()
    imp2 = int(linesp[9])
    mp2_name = linesp[10].strip()
    # linesp[11] is '|'
    iaa1 = int(linesp[12])
    iaa2 = int(linesp[13])
    # linesp[14] is '|'
    energy = float(linesp[15])
    
    ires1_u = cafemp2mpinfo[imp1].ires_u
    imp1_u = cafemp2mpinfo[imp1].imp_u
    ires2_u = cafemp2mpinfo[imp2].ires_u
    imp2_u = cafemp2mpinfo[imp2].imp_u
    
    f_out.write('%3i %6i %6i %4s %6i %6i %3s | %3i %6i %6i %4s %6i %6i %3s | %6i %6i | %f\n' 
                % (ichain1,ires1,ires1_u,res1_name,imp1,imp1_u,mp1_name,
                   ichain2,ires2,ires2_u,res2_name,imp2,imp2_u,mp2_name,
                   iaa1,iaa2,energy))
    
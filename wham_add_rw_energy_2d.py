#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys

if len(sys.argv) != 6:
    print('Usage: SCRIPT [IN .eta] [IN opt.ene] [original factor] [new factor] [OUT]')
    sys.exit(2)

f_eta = open(sys.argv[1],'r')
f_ene = open(sys.argv[2],'r')
factor = float(sys.argv[4]) / float(sys.argv[3])
f_out = open(sys.argv[-1],'w')

for l_ene in f_ene:
    ene = float(l_ene.strip())
    
    lsp_eta = f_eta.readline().split()
    irest1d = int(lsp_eta[0])
    step = int(lsp_eta[1])
    eta = float(lsp_eta[2])
    f_out.write('%1i %10i %15.6f %15.6f %15.6f\n' % (irest1d, step, eta, ene, ene*factor))

    lsp_eta = f_eta.readline().split()
    irest1d = int(lsp_eta[0])
    step = int(lsp_eta[1])
    eta = float(lsp_eta[2])
    f_out.write('%1i %10i %15.6f\n' % (irest1d, step, eta))
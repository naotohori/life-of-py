#!/usr/bin/env python

'''
The numerical table 'btmap.dat' contains
values multiplied by 0.6 to the original values 
in the matrix B in Betancourt and Thirumalai.

Amino acid names will be capitalized, e.g. glu -> GLU.
'''

f_out = open('BetancourtThirumalai.dat', 'w')

for l in open('btmap.dat'):
    lsp = l.split()
    aa1 = lsp[0].upper()
    aa2 = lsp[1].upper()
    e = float(lsp[2])

    if e == 0:
        f_out.write('%s %s %5.2f\n' % (aa1, aa2, 0.0))
    else:
        f_out.write('%s %s %5.2f\n' % (aa1, aa2, e/0.6))


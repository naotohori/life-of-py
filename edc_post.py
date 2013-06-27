#!/usr/bin/env python

import sys
from file_io.ninfo import NinfoFile
from elements.ninfo import NinfoSet

if len(sys.argv) < 7 or len(sys.argv)%1 != 0 :
    print ('Usage: SCRIPT [edc residue file] [ninfo file] [edc Amber mdout file] [(ID begin, ID end)...] [out file]')
    sys.exit(2)

f_residue = open(sys.argv[1], 'r')
f_ninfo = NinfoFile(sys.argv[2])
f_ninfo.open_to_read()
ninfo = NinfoSet()
f_ninfo.read_all(ninfo)
f_ninfo.close()
f_amber = open(sys.argv[3] ,'r')
#f_ninfo_out = NinfoFile(sys.argv[-2])
#f_ninfo_out.open_to_write()
f_out = open(sys.argv[-1], 'w')

# Reading arguments (ID begin, ID end)...
id_pairs =[]
for (i,arg) in enumerate(sys.argv[4:-1]) :
    if i%2 == 0:
        pair = (int(arg), )
    else :
        id_pairs.append( pair + (int(arg), ))

class mpinfo :
    def __init__(self):
        self.ichain = 0
        self.ires = 0
        self.imp = 0
        self.ires_l = 0
        self.imp_l = 0
        self.res_name = ''
        self.mp_name = ''
        self.ires_amber = 0
    def show(self):
        print self.ichain, self.ires, self.imp, self.imp, self.ires_l, self.imp_l, self.res_name, self.mp_name, self.ires_amber
        
cafemp2mpinfo = {}
for line in f_residue :
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
    cafemp2mpinfo[mp.imp] = mp
    
eambers = {}
ires_max = 0
for line in f_amber :
    if line.find('TDC') == -1:
        continue
    if line.find("0.000     0.000     0.000     0.000     0.000") != -1:
        continue
    ires1 = int(line[3:11])
    ires2 = int(line[13:20])
    if ires1 > ires_max :
        ires_max = ires1
    if ires2 > ires_max :
        ires_max = ires2
    #======================================================================
    #0         1         2         3         4         5         6         7         8
    #012345678901234567890123456789012345678901234567890123456789012345678901234567890
    #TDC     497->   1375     0.000    -0.005     6.526    -6.518     0.000
    #TDC     497->   1378     0.000     0.000     5.650    -5.645     0.000    
    #    resid1 ->resid2 |internal |vdw      |eel      |pol      |sas
    #======================================================================
    E_local = float(line[20:30])
    E_vdw   = float(line[30:40])
    E_eel   = float(line[40:50])
    E_pol   = float(line[50:60])
    E_sas   = float(line[60:70])
    #print ires1, ires2, E_local, E_vdw, E_eel, E_pol, E_sas
    if ires1 == ires2 or ires1 < ires2 :
        pair = (ires1, ires2)
    else:
        pair = (ires2, ires1)
    if pair in eambers:
        eambers[pair] += E_vdw + E_eel + E_pol + E_sas
    else :
        eambers[pair] = E_vdw + E_eel + E_pol + E_sas
#    eambers[(ires1, ires2)] = E_vdw + E_eel + E_pol + E_sas

#for i in xrange(1, ires_max+1) :
#    for j in xrange(i+3, ires_max+1) :
#        if (i,j) in eambers :
#            print i, j, eambers[(i,j)]
#        if (j,i) in eambers :
#            print j, i, eambers[(j,i)]

#f_out.write('#c1   res1 res1    mp1 mp1 |  c2   res2 res2    mp2 mp2 | amber1 amber2 | energy\n')
edcinfo = {}
for con in ninfo.contacts :
    imp1 = con.imp1
    imp2 = con.imp2
    mpinfo1 = cafemp2mpinfo[imp1]
    mpinfo2 = cafemp2mpinfo[imp2]
    
    ## intra-chain is not target!!
    #if cafemp2mpinfo[imp1].ichain == cafemp2mpinfo[imp2].ichain :
    #    continue
    
    flg_target = False
    for id_pair in id_pairs :
        if (   (imp1 >= id_pair[0] and imp1 <= id_pair[1])
            or (imp2 >= id_pair[0] and imp2 <= id_pair[1])) :
            flg_target = True
    if not flg_target :
        continue
        
    ires_amber1 = cafemp2mpinfo[imp1].ires_amber
    ires_amber2 = cafemp2mpinfo[imp2].ires_amber
    #print imp1, imp2, ires_amber1, ires_amber2
    if ires_amber1 <= ires_amber2 :
        energy = eambers[(ires_amber1, ires_amber2)]
    else:
        energy = eambers[(ires_amber2, ires_amber1)]
    #if energy >= 0.0 :
    #    con.factor = 0.0
    #else :
    #    con.factor = - energy
#    f_out.write('%3i %6i %4s %6i %3s | %3i %6i %4s %6i %3s | %6i %6i | %f\n' %
#          (mpinfo1.ichain, mpinfo1.ires, mpinfo1.res_name, imp1, mpinfo1.mp_name,
#           mpinfo2.ichain, mpinfo2.ires, mpinfo2.res_name, imp2, mpinfo2.mp_name,
#           ires_amber1, ires_amber2, energy) )
    if (con.iunit1, con.iunit2) in edcinfo:
        edcinfo[(con.iunit1, con.iunit2)].append(
                                (mpinfo1.ichain, mpinfo1.ires, mpinfo1.res_name, imp1, mpinfo1.mp_name,
                                 mpinfo2.ichain, mpinfo2.ires, mpinfo2.res_name, imp2, mpinfo2.mp_name,
                                 ires_amber1, ires_amber2, energy) )
    else :
        edcinfo[(con.iunit1, con.iunit2)] = []
        edcinfo[(con.iunit1, con.iunit2)].append(
                                (mpinfo1.ichain, mpinfo1.ires, mpinfo1.res_name, imp1, mpinfo1.mp_name,
                                 mpinfo2.ichain, mpinfo2.ires, mpinfo2.res_name, imp2, mpinfo2.mp_name,
                                 ires_amber1, ires_amber2, energy) )

for unitpair in edcinfo.keys() :
    f_out.write('#unit %i %i\n' % unitpair)
    f_out.write('#c1   res1 res1    mp1 mp1 |  c2   res2 res2    mp2 mp2 | amber1 amber2 | energy\n')
    for info in edcinfo[unitpair] :
        f_out.write('%3i %6i %4s %6i %3s | %3i %6i %4s %6i %3s | %6i %6i | %f\n' % info)
    f_out.write('\n')
        
#f_ninfo_out.write_all(ninfo)
#f_ninfo_out.close()
f_out.close()

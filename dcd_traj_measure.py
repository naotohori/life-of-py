#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17
@author: Naoto Hori
'''

import sys
import math
import copy
import re
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet
from numpy import array, dot, arccos
from numpy.linalg import norm

class DCD_END(Exception):
    pass

if len(sys.argv) != 3:
    print ' Usage: % SCRIPT [input DCD] [command file] '
    print '      command file is...'
    print '      ##########################################'
    print '      #frameskip 5                             #'
    print '      #display 1000                            #'
    print '      #out_time 1                              #'
    print '      #          (0:nothing, 1:frame, 2:time)  #'
    print '      #                                        #'
    print '      #def_mp  mp1  34                         #'
    print '      #def_mp  mp2  189                        #'
    print '      #def_mp   A   200                        #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #def_com com1 2 3 4 5 - 10 12 14 - 20    #'
    print '      #def_com com2 2 3 4 5 - 10 12 14 - 20    #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #def_mol mol1 2 3 4 5 - 10 12 14 - 20    #'
    print '      #def_mol mol2 2 3 4 5 - 10 12 14 - 20    #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #distance  1   2                         #'
    print '      #distance  3   A                         #'
    print '      #distance 100  com2                      #'
    print '      #distance com1 com2                      #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #angle     10  11  12  13                #'
    print '      #angle     com1  11  com2  13            #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #locaxis com1 com2 com3                  #'
    print '      # !location of com1 on axis (com2->com3) #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #def_ninfo ninfo1 filename               #'
    print '      #q all ninfo1                            #'
    print '      #q contact ninfo1                        #'
    print '      #q basepair ninfo1                       #'
    print '      #q basestack ninfo1                      #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #dock mol mol cutdist cutnum cutcom      #'
    print '      #dock mol mol   10.0    4     100.0      #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #enc contact mol mol 10.0 100.0          #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #polar com1                              #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #rg    mol                               #'
    print '      #     :    :   :                         #'
    print '      #                                        #'
    print '      #dt_hb name r theta1 theta2 phi phi1 phi2#'
    print '      ##########################################'
    sys.exit(2)

JUDGE_CONTACT = 1.2
JUDGE_CONTACT_2 = JUDGE_CONTACT ** 2.0
DISPLAY_PROCEDURE_DCD_STEP = 0

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()
nmp = dcd.get_header().nmp_real

f_cmd = open(sys.argv[2], 'r')

# Read command file
cmds = []
defs = {}
ninfo_files = {}
frameskip = 1
display = 0
out_time = 0
for line in f_cmd :
    if line.find('#') != -1:
        continue
    linesp = line.split()
    if len(linesp) == 0:
        continue
        
    if linesp[0] == 'frameskip':
        if int(linesp[1]) > 0:
            frameskip = int(linesp[1])
        else :
            frameskip = 1
            
    elif linesp[0] == 'display' :
        if int(linesp[1]) > 0:
            display = int(linesp[1])
            
    elif linesp[0] == 'out_time' :
        out_time = int(linesp[1])
        if not out_time in (0,1,2):
            print ('out_time is should be 0,1 or 2\n>>>>' + line)
            sys.exit(2)
            
    elif linesp[0] == 'distance':
        re1 = re.match('^\d+$', linesp[1])
        re2 = re.match('^\d+$', linesp[2])
        if re1 :
            linesp[1] = int(linesp[1])
        if re2 :
            linesp[2] = int(linesp[2])
        if re1 and re2 :
            if int(linesp[1]) > nmp or int(linesp[2]) > nmp :
                print ('id is larger than nmp\n in the line:' + line)
                sys.exit(2)
        cmds.append(tuple(linesp))

    elif linesp[0] == 'angle':
        if len(linesp) != 5:
            print ('not enough arguments for angle in the line:' + line)
            sys.exit(2)
        #linesp[1] = int(linesp[1])
        #linesp[2] = int(linesp[2])
        #linesp[3] = int(linesp[3])
        #linesp[4] = int(linesp[4])
        cmds.append(tuple(linesp))

    elif linesp[0] == 'delta_angle':
        if len(linesp) != 3:
            print ('not enough arguments for delta_angle in the line:' + line)
            sys.exit(2)
        linesp[1] = int(linesp[1])
        linesp[2] = int(linesp[2])
        cmds.append(tuple(linesp))

    elif linesp[0] == 'def_com':
        name = linesp[1]
        mps = []
        flg_sequence = False
        for arg in linesp[2:] :
            if flg_sequence :
                mps.extend([imp for imp in xrange(mps[-1] + 1, int(arg) + 1)])
                flg_sequence = False
            if arg == '-' :
                flg_sequence = True
                continue
            mps.append(int(arg))
        defs[name] = ('com', mps)

    elif linesp[0] == 'def_mp':
        name = linesp[1]
        defs[name] = ('mp', [int(linesp[2])])
        
    elif linesp[0] == 'def_mol':
        name = linesp[1]
        mps = []
        flg_sequence = False
        for arg in linesp[2:] :
            if flg_sequence :
                mps.extend([imp for imp in xrange(mps[-1] + 1, int(arg) + 1)])
                flg_sequence = False
            if arg == '-' :
                flg_sequence = True
                continue
            mps.append(int(arg))
        defs[name] = ('mol', mps)

    elif linesp[0] == 'locaxis':
        re1 = re.match('^\d+$', linesp[1])
        re2 = re.match('^\d+$', linesp[2])
        re3 = re.match('^\d+$', linesp[3])
        if re1 :
            linesp[1] = int(linesp[1])
        if re2 :
            linesp[2] = int(linesp[2])
        if re3 :
            linesp[3] = int(linesp[3])
        cmds.append(tuple(linesp))

    elif linesp[0] == 'def_ninfo':
        name = linesp[1]
        filename = linesp[2]
        nf = NinfoFile(filename)
        nf.open_to_read()
        ns = NinfoSet()
        nf.read_all(ns)
        nf.close()
        defs[name] = ('ninfo', filename, ns)

    elif linesp[0] == 'q':
        cmds.append(tuple(linesp))
        
    elif linesp[0] == 'dock':
        cmds.append(tuple(linesp))
    
    elif linesp[0] == 'enc':
        cmds.append(tuple(linesp))
        
    elif linesp[0] == 'polar':
        cmds.append(tuple(linesp))
        
    elif linesp[0] == 'rg':
        cmds.append(tuple(linesp))

    elif linesp[0] == 'dt_hb':
        cmds.append(tuple(linesp))
#    elif linesp[0] == 'distance_com_com' :
#        linesp[1] = int(linesp[1])
#        linesp[2] = int(linesp[2])
#        linesp[3] = int(linesp[3])
#        linesp[4] = int(linesp[4])
#        cmds.append(tuple(linesp))

out_files = []
for cmd in cmds :
    if cmd[0] == 'distance' :
        filename = 'dist'
        if cmd[1] in defs :
            filename += '_%s' % cmd[1]
        else :
            filename += '_%i' % cmd[1]
        if cmd[2] in defs :
            filename += '_%s' % cmd[2]
        else :
            filename += '_%i' % cmd[2]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# distance between 1 and 2\n')
        out_files[-1].write('# 1:')
        if cmd[1] in defs:
            out_files[-1].write(' %s' % defs[cmd[1]][0])
            for imp in defs[cmd[1]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[1])
        out_files[-1].write('# 2:')
        if cmd[2] in defs:
            out_files[-1].write(' %s' % defs[cmd[2]][0])
            for imp in defs[cmd[2]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[2])
    elif cmd[0] == 'angle' :
        #out_files.append(open('%s_%i_%i_%i_%i.out' % cmd, 'w'))
        out_files.append(open('%s_%s_%s_%s_%s.out' % cmd, 'w'))
        out_files[-1].write('# angle between line ij and kl\n')
        out_files[-1].write('# i:')
        if cmd[1] in defs:
            out_files[-1].write(' %s' % defs[cmd[1]][0])
            for imp in defs[cmd[1]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[1])
        out_files[-1].write('# j:')
        if cmd[2] in defs:
            out_files[-1].write(' %s' % defs[cmd[2]][0])
            for imp in defs[cmd[2]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[2])
        out_files[-1].write('# k:')
        if cmd[3] in defs:
            out_files[-1].write(' %s' % defs[cmd[3]][0])
            for imp in defs[cmd[3]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[3])
        out_files[-1].write('# l:')
        if cmd[4] in defs:
            out_files[-1].write(' %s' % defs[cmd[4]][0])
            for imp in defs[cmd[4]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[4])
    elif cmd[0] == 'delta_angle' :
        out_files.append(open('%s_%i_%i.out' % cmd, 'w'))
    elif cmd[0] == 'locaxis' :
        filename = 'locaxis'
        if cmd[1] in defs :
            filename += '_%s_on' % cmd[1]
        else :
            filename += '_%i_on' % cmd[1]
        if cmd[2] in defs :
            filename += '_%s' % cmd[2]
        else :
            filename += '_%i' % cmd[2]
        if cmd[3] in defs :
            filename += '_%s' % cmd[3]
        else :
            filename += '_%i' % cmd[3]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# location of 1 on the axis 1 to 2\n')
        out_files[-1].write('# 1:')
        if cmd[1] in defs:
            out_files[-1].write(' %s' % defs[cmd[1]][0])
            for imp in defs[cmd[1]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[1])
        out_files[-1].write('# 2:')
        if cmd[2] in defs:
            out_files[-1].write(' %s' % defs[cmd[2]][0])
            for imp in defs[cmd[2]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[2])
        out_files[-1].write('# 3:')
        if cmd[3] in defs:
            out_files[-1].write(' %s' % defs[cmd[3]][0])
            for imp in defs[cmd[3]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[3])
    elif cmd[0] == 'q' :
        filename = 'q'
        filename += '_%s_%s' % cmd[1:3]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# q-score based on ninfo file1\n')
        out_files[-1].write('# type: %s' % cmd[1])
        out_files[-1].write('# file1: ')
        if cmd[2] in defs:
            out_files[-1].write(defs[cmd[2]][1])
        else :
            print ('%s is not in definitions\n' % cmd[2])
            sys.exit(2)
        out_files[-1].write('\n')
    elif cmd[0] == 'dock' :
        filename = 'dock'
        if cmd[1] in defs :
            filename += '_%s' % cmd[1]
        else :
            filename += '_%i' % cmd[1]
        if cmd[2] in defs :
            filename += '_%s' % cmd[2]
        else :
            filename += '_%i' % cmd[2]
        filename += '_%s' % cmd[3]
        filename += '_%s' % cmd[4]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# docking between 1 and 2\n')
        out_files[-1].write('# 1:')
        if cmd[1] in defs:
            out_files[-1].write(' %s' % defs[cmd[1]][0])
            for imp in defs[cmd[1]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[1])
        out_files[-1].write('# 2:')
        if cmd[2] in defs:
            out_files[-1].write(' %s' % defs[cmd[2]][0])
            for imp in defs[cmd[2]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[2])
        out_files[-1].write('# 3: %s\n' % cmd[3])
        out_files[-1].write('# 4: %s\n' % cmd[4])
        
    elif cmd[0] == 'enc' :
        filename = 'enc' + '_' + cmd[1]
        if cmd[2] in defs :
            filename += '_%s' % cmd[2]
        else :
            filename += '_%i' % cmd[2]
        if cmd[3] in defs :
            filename += '_%s' % cmd[3]
        else :
            filename += '_%i' % cmd[3]
        filename += '_%s' % cmd[4]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# encounter')
        out_files[-1].write(' %s ' % cmd[1])
        out_files[-1].write('between 1 and 2\n')
        out_files[-1].write('# 1:')
        if cmd[2] in defs:
            out_files[-1].write(' %s' % defs[cmd[2]][0])
            for imp in defs[cmd[2]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[2])
        out_files[-1].write('# 2:')
        if cmd[3] in defs:
            out_files[-1].write(' %s' % defs[cmd[3]][0])
            for imp in defs[cmd[3]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[3])
        out_files[-1].write('# 3 cutoff distance : %s\n' % cmd[4])
        out_files[-1].write('# 4 cutoff distance for skip : %s\n' % cmd[5])
    elif cmd[0] == 'polar':
        filename = 'polar'
        if cmd[1] in defs :
            filename += '_%s' % cmd[1]
        else :
            filename += '_%i' % cmd[1]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# Spherical polar coordinates')
        out_files[-1].write('# targets:')
        if cmd[1] in defs:
            out_files[-1].write(' %s' % defs[cmd[1]][0])
            for imp in defs[cmd[1]][1] :
                out_files[-1].write(' %i' % imp)
            out_files[-1].write('\n')
        else :
            out_files[-1].write('%i\n' % cmd[1])
        out_files[-1].write('%8s %6s %6s %7s %7s\n' % ('# r','theta','phi','theta[deg]','phi[deg]'))
    elif cmd[0] == 'rg':
        filename = 'Rg'
        #if cmd[1] in defs :
        #    filename += '_%s' % cmd[1]
        #else :
        #    filename += '_%i' % cmd[1]
        if not cmd[1] in defs :
            print ("For Rg, you have to define mol by 'def_mol'")
            sys.exit(2)
        filename += '_%s' % cmd[1]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# Rg')
        out_files[-1].write('# targets:')
        #if cmd[1] in defs:
        #    out_files[-1].write(' %s' % defs[cmd[1]][0])
        #    for imp in defs[cmd[1]][1] :
        #        out_files[-1].write(' %i' % imp)
        #    out_files[-1].write('\n')
        #else :
        #    out_files[-1].write('%i\n' % cmd[1])
        out_files[-1].write(' %s' % defs[cmd[1]][0])
        for imp in defs[cmd[1]][1] :
            out_files[-1].write(' %i' % imp)
        out_files[-1].write('\n')
        out_files[-1].write('%8s %6s %6s %7s %7s\n' % ('# r','theta','phi','theta[deg]','phi[deg]'))
        
        # Add com definition 
        name = cmd[1] + '_RgCom'
        mps = defs[cmd[1]][1]
        defs[name] = ('RgCom', mps)
    elif cmd[0] == 'dt_hb':
        filename = 'dthb_%s' % cmd[1]
        out_files.append(open(filename + '.out', 'w'))
        out_files[-1].write('# dt_hb\n')
        out_files[-1].write('# name: %s\n' % cmd[1])
        out_files[-1].write('# r: %s\n' % cmd[2])
        out_files[-1].write('# theta1: %s\n' % cmd[3])
        out_files[-1].write('# theta2: %s\n' % cmd[4])
        out_files[-1].write('# psi: %s\n' % cmd[5])
        out_files[-1].write('# psi1: %s\n' % cmd[6])
        out_files[-1].write('# psi2: %s\n' % cmd[7])
        
iframe = 0
flg_initial_data = True
try:
    while dcd.has_more_data() :
        iframe += 1
        if display != 0 :
            if (iframe % display == 0) :
                print ('%i frame done' % iframe)
            
        # Skip
        if frameskip > 1 and iframe != 0 :
            if iframe % frameskip != 0 :
                dcd.skip(1)
                continue
            
        # Read one step
        data = dcd.read_onestep()
        
        if flg_initial_data:
            data_ini = copy.deepcopy(data)
            flg_initial_data = False
            
        coms = {}
        for key, d in defs.items():
            if d[0] in ('com', 'RgCom'):
                com = [0.0, 0.0, 0.0]
                for imp in d[1] :
                    com[0] += data[imp - 1][0]
                    com[1] += data[imp - 1][1]
                    com[2] += data[imp - 1][2]
                coms[key] = [x / float(len(d[1])) for x in com]
        
        for (icmd, cmd) in enumerate(cmds):
            if cmd[0] == 'distance' :
                if cmd[1] in defs :
                    if defs[cmd[1]][0] == 'com' :
                        xyz_i = coms[cmd[1]]
                    elif defs[cmd[1]][0] == 'mp' :
                        xyz_i = data[defs[cmd[1]][1][0] - 1]
                else:
                    xyz_i = data[cmd[1] - 1]
                if cmd[2] in defs :
                    if defs[cmd[2]][0] == 'com' :
                        xyz_j = coms[cmd[2]]
                    elif defs[cmd[2]][0] == 'mp' :
                        xyz_j = data[defs[cmd[2]][1][0] - 1]
                else:
                    xyz_j = data[cmd[2] - 1]
                d = math.sqrt((xyz_i[0] - xyz_j[0]) ** 2
                              + (xyz_i[1] - xyz_j[1]) ** 2
                              + (xyz_i[2] - xyz_j[2]) ** 2)
                out_files[icmd].write('%12.5f\n' % d)
 
#            elif cmd[0] == 'angle' :
#                i = cmd[1] - 1
#                j = cmd[2] - 1
#                k = cmd[3] - 1
#                l = cmd[4] - 1
#                v1 = [data[j][0] - data[i][0],
#                      data[j][1] - data[i][1],
#                      data[j][2] - data[i][2]]
#                v2 = [data[k][0] - data[i][0],
#                      data[k][1] - data[i][1],
#                      data[k][2] - data[i][2]]
#                v3 = [data[l][0] - data[i][0],
#                      data[l][1] - data[i][1],
#                      data[l][2] - data[i][2]]
#                v21 = array([v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]])
#    #            v32 = array( [v3[0]-v2[0], v3[1]-v2[1], v3[2]-v2[2]] )
#    #            c11 = v21[0]*v21[0] + v21[1]*v21[1] + v21[2]*v21[2]
#    #            c22 = v32[0]*v32[0] + v32[1]*v32[1] + v32[2]*v32[2]
#    #            c21 = v32[0]*v21[0] + v32[1]*v21[1] + v32[2]*v21[2]
#    #            co_theta = - c21 / math.sqrt(c11 * c22)
#    #            theta = math.acos(co_theta)
#                v23 = array([v2[0] - v3[0], v2[1] - v3[1], v2[2] - v3[2]])
#                theta = arccos(dot(v21, v23) / norm(v21) / norm(v23))
#                out_files[icmd].write('%12.5f %12.5f\n' % (theta, math.degrees(theta)))
      
            elif cmd[0] == 'angle' :
                if cmd[1] in defs :
                    if defs[cmd[1]][0] == 'com' :
                        com = [0.0, 0.0, 0.0]
                        n = 0
                        for imp in defs[cmd[1]][1] :
                            n += 1
                            com[0] += data[imp - 1][0]
                            com[1] += data[imp - 1][1]
                            com[2] += data[imp - 1][2]
                        xyz_i = [com[0] / float(n), com[1] / float(n), com[2] / float(n)]
                    elif defs[cmd[1]][0] == 'mp' :
                        xyz_i = data[defs[cmd[1]][1][0] - 1]
                else:
                    xyz_i = data[int(cmd[1]) - 1]
                if cmd[2] in defs :
                    if defs[cmd[2]][0] == 'com' :
                        com = [0.0, 0.0, 0.0]
                        n = 0
                        for imp in defs[cmd[2]][1] :
                            n += 1
                            com[0] += data[imp - 1][0]
                            com[1] += data[imp - 1][1]
                            com[2] += data[imp - 1][2]
                        xyz_j = [com[0] / float(n), com[1] / float(n), com[2] / float(n)]
                    elif defs[cmd[2]][0] == 'mp' :
                        xyz_j = data[defs[cmd[2]][1][0] - 1]
                else:
                    xyz_j = data[int(cmd[2]) - 1]
                if cmd[3] in defs :
                    if defs[cmd[3]][0] == 'com' :
                        com = [0.0, 0.0, 0.0]
                        n = 0
                        for imp in defs[cmd[3]][1] :
                            n += 1
                            com[0] += data[imp - 1][0]
                            com[1] += data[imp - 1][1]
                            com[2] += data[imp - 1][2]
                        xyz_k = [com[0] / float(n), com[1] / float(n), com[2] / float(n)]
                    elif defs[cmd[3]][0] == 'mp' :
                        xyz_k = data[defs[cmd[3]][1][0] - 1]
                else:
                    xyz_k = data[int(cmd[3]) - 1]
                if cmd[4] in defs :
                    if defs[cmd[4]][0] == 'com' :
                        com = [0.0, 0.0, 0.0]
                        n = 0
                        for imp in defs[cmd[4]][1] :
                            n += 1
                            com[0] += data[imp - 1][0]
                            com[1] += data[imp - 1][1]
                            com[2] += data[imp - 1][2]
                        xyz_l = [com[0] / float(n), com[1] / float(n), com[2] / float(n)]
                    elif defs[cmd[4]][0] == 'mp' :
                        xyz_l = data[defs[cmd[4]][1][0] - 1]
                else:
                    xyz_l = data[int(cmd[4]) - 1]
                vij = array([xyz_j[0] - xyz_i[0],
                             xyz_j[1] - xyz_i[1],
                             xyz_j[2] - xyz_i[2]])
                vkl = array([xyz_l[0] - xyz_k[0],
                             xyz_l[1] - xyz_k[1],
                             xyz_l[2] - xyz_k[2]])
                theta = arccos(dot(vij, vkl) / norm(vij) / norm(vkl))
                out_files[icmd].write('%12.5f %12.5f\n' % (theta, math.degrees(theta)))

            elif cmd[0] == 'delta_angle' :
                i = cmd[1] - 1
                j = cmd[2] - 1
                v1 = [data_ini[j][0] - data_ini[i][0],
                      data_ini[j][1] - data_ini[i][1],
                      data_ini[j][2] - data_ini[i][2]]
                v2 = [data[i][0] - data_ini[i][0],
                      data[i][1] - data_ini[i][1],
                      data[i][2] - data_ini[i][2]]
                v3 = [data[j][0] - data_ini[i][0],
                      data[j][1] - data_ini[i][1],
                      data[j][2] - data_ini[i][2]]
                v21 = array([v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]])
                v23 = array([v2[0] - v3[0], v2[1] - v3[1], v2[2] - v3[2]])
                theta = arccos(dot(v21, v23) / norm(v21) / norm(v23))
                out_files[icmd].write('%12.5f %12.5f\n' % (theta, math.degrees(theta)))
    
            elif cmd[0] == 'locaxis' :
                if cmd[1] in defs :
                    if defs[cmd[1]][0] == 'com' :
                        xyz_i = coms[cmd[1]]
                    elif defs[cmd[1]][0] == 'mp' :
                        xyz_i = data[defs[cmd[1]][1][0] - 1]
                else:
                    xyz_i = data[cmd[1] - 1]
                if cmd[2] in defs :
                    if defs[cmd[2]][0] == 'com' :
                        xyz_a = coms[cmd[2]]
                    elif defs[cmd[2]][0] == 'mp' :
                        xyz_a = data[defs[cmd[2]][1][0] - 1]
                else:
                    xyz_a = data[cmd[2] - 1]
                if cmd[3] in defs :
                    if defs[cmd[3]][0] == 'com' :
                        xyz_b = coms[cmd[3]]
                    elif defs[cmd[3]][0] == 'mp' :
                        xyz_b = data[defs[cmd[3]][1][0] - 1]
                else:
                    xyz_b = data[cmd[3] - 1]
                vAB = [xyz_b[0] - xyz_a[0],
                       xyz_b[1] - xyz_a[1],
                       xyz_b[2] - xyz_a[2]]
                vAi = [xyz_i[0] - xyz_a[0],
                       xyz_i[1] - xyz_a[1],
                       xyz_i[2] - xyz_a[2]]
                location = dot(vAB, vAi) / norm(vAB)
                out_files[icmd].write('%12.5f\n' % location)
    
            elif cmd[0] == 'q' :
                type = cmd[1]
                if defs[cmd[2]][0] != 'ninfo' :
                    print ('definition error : %s is not ninfo file\n' % cmd[1])
                    sys.exit(2)
                ns = defs[cmd[2]][2]
                n_denominator = 0
                n_con = 0
                if type in ('contact', 'all') :
                    n_denominator += len(ns.contacts)
                    for ni in ns.contacts :
                        imp1 = ni.imp1
                        imp1 = ni.imp2
                        xyz_i = data[ni.imp1 - 1]
                        xyz_j = data[ni.imp2 - 1]
                        native_dist_2 = ni.native ** 2.0
                        dist2 = ((xyz_i[0] - xyz_j[0]) ** 2
                                + (xyz_i[1] - xyz_j[1]) ** 2
                                + (xyz_i[2] - xyz_j[2]) ** 2)
                        if (dist2 < JUDGE_CONTACT_2 * native_dist_2) :
                            n_con += 1
                if type in ('basepair', 'all') :
                    n_denominator += len(ns.basepairs)
                    for ni in ns.basepairs :
                        imp1 = ni.imp1
                        imp1 = ni.imp2
                        xyz_i = data[ni.imp1 - 1]
                        xyz_j = data[ni.imp2 - 1]
                        native_dist_2 = ni.native ** 2.0
                        dist2 = ((xyz_i[0] - xyz_j[0]) ** 2
                                + (xyz_i[1] - xyz_j[1]) ** 2
                                + (xyz_i[2] - xyz_j[2]) ** 2)
                        if (dist2 < JUDGE_CONTACT_2 * native_dist_2) :
                            n_con += 1
                if type in ('basestack', 'all') :
                    n_denominator += len(ns.basestacks)
                    for ni in ns.basestacks :
                        imp1 = ni.imp1
                        imp1 = ni.imp2
                        xyz_i = data[ni.imp1 - 1]
                        xyz_j = data[ni.imp2 - 1]
                        native_dist_2 = ni.native ** 2.0
                        dist2 = ((xyz_i[0] - xyz_j[0]) ** 2
                                + (xyz_i[1] - xyz_j[1]) ** 2
                                + (xyz_i[2] - xyz_j[2]) ** 2)
                        if (dist2 < JUDGE_CONTACT_2 * native_dist_2) :
                            n_con += 1
                out_files[icmd].write('%20i %12.5f\n' % (n_con, n_con / float(n_denominator)))
                
            elif cmd[0] == 'dock' :
                com_i = [0.0, 0.0, 0.0]
                n = 0
                for imp in defs[cmd[1]][1] :
                    n += 1
                    com_i[0] += data[imp - 1][0]
                    com_i[1] += data[imp - 1][1]
                    com_i[2] += data[imp - 1][2]
                com_i = [x / float(n) for x in com_i]
                com_j = [0.0, 0.0, 0.0]
                n = 0
                for imp in defs[cmd[2]][1] :
                    n += 1
                    com_j[0] += data[imp - 1][0]
                    com_j[1] += data[imp - 1][1]
                    com_j[2] += data[imp - 1][2]
                com_j = [x / float(n) for x in com_j]
                d = math.sqrt((com_i[0] - com_j[0]) ** 2
                             +(com_i[1] - com_j[1]) ** 2
                             +(com_i[2] - com_j[2]) ** 2)
                #if d > 100.0 :
                if d > float(cmd[5]) :
                    break
                mps1 = []
                mps2 = []
                if cmd[1] in defs :
                    if defs[cmd[1]][0] == 'mol' :
                        mps1 = defs[cmd[1]][1]
                    elif defs[cmd[1]][0] == 'mp' :
                        pass #error
                else:
                    pass #error
                if cmd[2] in defs :
                    if defs[cmd[2]][0] == 'mol' :
                        mps2 = defs[cmd[2]][1]
                    elif defs[cmd[2]][0] == 'mp' :
                        pass #error
                else:
                    pass #error
                cutdist = float(cmd[3])
                cutnum = float(cmd[4])
                within = []
                for mp1 in mps1 :
                    for mp2 in mps2 :
                        d = math.sqrt((data[mp1-1][0] - data[mp2-1][0]) ** 2
                                    + (data[mp1-1][1] - data[mp2-1][1]) ** 2
                                    + (data[mp1-1][2] - data[mp2-1][2]) ** 2)
                        if d <= cutdist:
                            within.append((mp1,mp2))
                if len(within) >= cutnum :
                    out_files[icmd].write('#%i\n' % iframe)
                    for pair in within :
                        out_files[icmd].write('%i %i\n' % pair)
                    raise DCD_END
                
            elif cmd[0] == 'enc' :
                if cmd[1] != 'contact':
                    print ('Error: unknown type in "enc"')
                    continue
                cutoff_skip = float(cmd[5])  # 100.0
                com_i = [0.0, 0.0, 0.0]
                n = 0
                for imp in defs[cmd[2]][1] :
                    n += 1
                    com_i[0] += data[imp - 1][0]
                    com_i[1] += data[imp - 1][1]
                    com_i[2] += data[imp - 1][2]
                com_i = [x / float(n) for x in com_i]
                com_j = [0.0, 0.0, 0.0]
                n = 0
                for imp in defs[cmd[3]][1] :
                    n += 1
                    com_j[0] += data[imp - 1][0]
                    com_j[1] += data[imp - 1][1]
                    com_j[2] += data[imp - 1][2]
                com_j = [x / float(n) for x in com_j]
                d = math.sqrt((com_i[0] - com_j[0]) ** 2
                             +(com_i[1] - com_j[1]) ** 2
                             +(com_i[2] - com_j[2]) ** 2)
                if d > cutoff_skip:
                    out_files[icmd].write('%12i %8i %8i %8i\n' % (iframe,0,0,0))
                    break
                mps1 = []
                mps2 = []
                if cmd[2] in defs :
                    if defs[cmd[2]][0] == 'mol' :
                        mps1 = defs[cmd[2]][1]
                    elif defs[cmd[2]][0] == 'mp' :
                        pass #error
                else:
                    pass #error
                if cmd[3] in defs :
                    if defs[cmd[3]][0] == 'mol' :
                        mps2 = defs[cmd[3]][1]
                    elif defs[cmd[3]][0] == 'mp' :
                        pass #error
                else:
                    pass #error
                cutdist2 = float(cmd[4]) ** 2
                within_1 = set()
                within_2 = set()
                n_within = 0
                for mp1 in mps1 :
                    for mp2 in mps2 :
                        d2 = ((data[mp1-1][0] - data[mp2-1][0]) ** 2
                             +(data[mp1-1][1] - data[mp2-1][1]) ** 2
                             +(data[mp1-1][2] - data[mp2-1][2]) ** 2)
                        if d2 <= cutdist2:
                            within_1.add(mp1)
                            within_2.add(mp2)
                            n_within += 1
                            
                out_files[icmd].write('%12i %8i %8i %8i\n' 
                                      % (iframe,len(within_1),len(within_2),
                                         n_within))
            elif cmd[0] == 'polar':
                if cmd[1] in defs :
                    if defs[cmd[1]][0] == 'com' :
                        xyz = coms[cmd[1]]
                    elif defs[cmd[1]][0] == 'mp' :
                        xyz = data[defs[cmd[1]][1][0] - 1]
                else:
                    xyz = data[cmd[1] - 1]
                dist = math.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
                theta = math.acos(xyz[2]/dist)
                phi = math.atan2(xyz[1], xyz[0])
                
                out_files[icmd].write('%8.2f %6.3f %6.3f %7.2f %7.2f\n' 
                                      % (dist,theta,phi,
                                         math.degrees(theta), math.degrees(phi)) )

            elif cmd[0] == 'rg' :
                name = cmd[1] + '_RgCom'
                if name in defs :
                    if defs[name][0] == 'RgCom' :
                        xyz_com = coms[name]
                    else:
                        print ('Error: Rg 1')
                        sys.exit(2)
                else:
                    print ('Error: Rg 2')
                    sys.exit(2)

                if cmd[1] in defs :
                    if defs[cmd[1]][0] == 'mol' :
                        mps1 = defs[cmd[1]][1]
                    else:
                        print ('Error: Rg 4')
                        sys.exit(2)
                else:
                    print ('Error: Rg 3')
                    sys.exit(2)
                s = 0.0
                for mp1 in mps1 :
                    s += ((data[mp1-1][0] - xyz_com[0]) ** 2
                        + (data[mp1-1][1] - xyz_com[1]) ** 2
                        + (data[mp1-1][2] - xyz_com[2]) ** 2)
                s = math.sqrt( s / float(len(mps1)) )
                out_files[icmd].write('%12.5f\n' % s)

            elif cmd[0] == 'dt_hb':
                mp1, mp2 = cmd[2].split('-')
                mp1 = int(mp1)
                mp2 = int(mp2)
                r = math.sqrt((data[mp1-1][0] - data[mp2-1][0]) ** 2
                             +(data[mp1-1][1] - data[mp2-1][1]) ** 2
                             +(data[mp1-1][2] - data[mp2-1][2]) ** 2)

                mp1, mp2, mp3 = cmd[3].split('-')
                xyz_i = data[int(mp1)-1]
                xyz_j = data[int(mp2)-1]
                xyz_k = data[int(mp3)-1]
                vij = array([xyz_j[0] - xyz_i[0],
                             xyz_j[1] - xyz_i[1],
                             xyz_j[2] - xyz_i[2]])
                vjk = array([xyz_k[0] - xyz_j[0],
                             xyz_k[1] - xyz_j[1],
                             xyz_k[2] - xyz_j[2]])
                theta1 = arccos(dot(vij, vjk) / norm(vij) / norm(vjk))

                mp1, mp2, mp3 = cmd[4].split('-')
                xyz_i = data[int(mp1)-1]
                xyz_j = data[int(mp2)-1]
                xyz_k = data[int(mp3)-1]
                vij = array([xyz_j[0] - xyz_i[0],
                             xyz_j[1] - xyz_i[1],
                             xyz_j[2] - xyz_i[2]])
                vjk = array([xyz_k[0] - xyz_j[0],
                             xyz_k[1] - xyz_j[1],
                             xyz_k[2] - xyz_j[2]])
                theta2 = arccos(dot(vij, vjk) / norm(vij) / norm(vjk))

                out_files[icmd].write('%7.3f %6.3f %6.3f\n' % (r, theta1, theta2))
                # math.degrees(theta)
                 
except DCD_END:
    pass

dcd.close()
for ifile in out_files :
    ifile.close()

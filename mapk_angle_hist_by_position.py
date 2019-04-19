#!/usr/bin/env python
#vim:fileencoding=UTF-8



import sys
import math
from numpy import histogram

if len(sys.argv) != 8:
    print('Usage: SCRIPT [angle data] [polar data] [theta from] [theta to] [phi from] [phi to] [output prefix]')
    print(' (theta and phi is in the unit of degree)')
    sys.exit(2)
    
file_in = open(sys.argv[1],'r')
file_pfx = sys.argv[-1]
file_pol = open(sys.argv[2])
theta_from = float(sys.argv[3])
theta_to = float(sys.argv[4])
phi_from = float(sys.argv[5])
phi_to = float(sys.argv[6])

file_out = open(file_pfx+"_hist.out",'w')
file_out.write('#PROGRAM: mapk_angle_hist_by_position.py\n')
file_out.write('#angle data: %s\n' % (sys.argv[1]))
file_out.write('#polar data: %s\n' % (sys.argv[2]))
file_out.write('#theta from %f to %f\n' % (theta_from,theta_to))
file_out.write('#phi from %f to %f\n' % (phi_from,phi_to))

COL_THETA = 2 - 1
COL_POL_THETA = 4 - 1
COL_POL_PHI = 5 - 1

#theta_bins = [x*5.0 for x in xrange(0,37)] #  5度
theta_bins = [x*10.0 for x in range(0,19)] # 10度
#theta_bins = [x*15.0 for x in xrange(0,13)] #  15度


theta = [] 
weight = []
num_ang = 0  # for check number of lines are consistent
num_pol = 0 

''' loop for file reading'''
for l in file_in:
    if l.find('#') != -1:
        continue
    num_ang += 1
    
    '''read next polar position'''
    l_pol = file_pol.readline()
    while l_pol.find('#') != -1:
        l_pol = file_pol.readline()
    num_pol += 1
    l_pol_sp = l_pol.split()
    pol_t = float(l_pol_sp[COL_POL_THETA])
    pol_p = float(l_pol_sp[COL_POL_PHI])
    
    ''' judge'''
    if pol_t < theta_from or theta_to < pol_t:
        continue
    if pol_p < phi_from or phi_to < pol_p:
        continue
    
    '''(accepted) => add angle data '''
    lsp = l.split()
    t = float(lsp[COL_THETA])
    theta.append(t)
    weight.append(1.0/math.sin(math.radians(t)))

if num_ang != num_pol:
    print('Error: angle data and polar coordinate data are inconsistent!')
    print('ABORT')
    sys.exit(2)
        
        
h, theta_edge = histogram(theta,bins=theta_bins)
hw, theta_edge = histogram(theta,weights=weight, bins=theta_bins)
#hd, theta_edge = histogram(theta,bins=theta_bins,density=True)

hsum = float(sum(h))
hwsum = sum(hw)

    
for i,x in enumerate(h):
    ang = (theta_edge[i]+theta_edge[i+1])*0.5
    jcb = math.sin(math.radians(ang))
    file_out.write('%8.3f %10.6f %10.6f %10i\n'
                    % (ang, x/hsum/jcb, hw[i]/hwsum, x))
    
file_out.close()
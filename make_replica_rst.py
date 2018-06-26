#!/usr/bin/env python

import sys
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 6:
    print ('Usage: SCRIPT [Nrep] [nmp] [step] [path/prefix] [output restart file]')
    sys.exit(2)

Nrep = int(sys.argv[1])
nmp  = int(sys.argv[2])
step = int(sys.argv[3])
file_path = sys.argv[4]

f_rst = open(sys.argv[5],'w')


## Replica
count_blanks = 0
f_rep = open(file_path+'.rep')
l_pre = ''
for l in open(file_path+'.rep'):
    if len(l.strip()) == 0:
        count_blanks += 1
        continue
    if count_blanks == 4:
        if l[0] == '#':
            continue
        l_pre = l
        if len(l.strip()) == 0:
            break

print ('# last line of .rep')
print (l_pre)

replica_labels = []
for label in l_pre.split()[1:]:
    replica_labels.append(int(label))

print ('# replica labels')
print (replica_labels)

if len(replica_labels) != Nrep:
    print ('Error: len(replica_labels) != Nrep')
    sys.exit(2)

f_rst.write('# replica\n')
f_rst.write('n_replica_all: %i\n' % (Nrep,))

for i, label in enumerate(replica_labels):
    f_rst.write('rep2lab: %5i %5i\n' % (i+1, label))

#### Step
f_rst.write('# step\n')
f_rst.write('istep_sim:  1\n')
f_rst.write('istep:  %i\n' % (step,))


#### coordinates of each replicas
for irep in range(1, Nrep+1):
    f_rst.write('# xyz_mp_rep\n')
    f_rst.write('grep: %i\n' % (irep,))

    dcd = DcdFile(file_path+'_%04i.dcd' % (irep,))
    dcd.open_to_read()
    header = dcd.read_header()
    while dcd.has_more_data() :
        data = dcd.read_onestep()

    if len(data) != nmp:
        print ('Error len(data) != nmp in replica %i. len(data)=%i' % (irep,len(data)))
        sys.exit(2)

    f_rst.write('nmp_all: %i\n' % (nmp,))
    for imp in range(nmp):
        f_rst.write('%f %f %f\n' % (data[imp][0], data[imp][1], data[imp][2]))

    f_rst.write('# velo_mp\n')
    f_rst.write('grep: %i\n' % (irep,))
    f_rst.write('nmp_real: %i\n' % (nmp,))
    for imp in range(nmp):
        f_rst.write('0.0 0.0 0.0\n')

    f_rst.write('# accel_mp\n')
    f_rst.write('grep: %i\n' % (irep,))
    f_rst.write('nmp_real: %i\n' % (nmp,))
    for imp in range(nmp):
        f_rst.write('0.0 0.0 0.0\n')

#!/usr/bin/env python

import sys

if len(sys.argv) != 6:
    print(' Usage: % SCRIPT [file] [data column (1,2,3,...)] [bin most left] [bin most right] [bin width]')
    sys.exit(2)

for i,arg in enumerate(sys.argv):
    print('#%i %s' % (i,arg))
print('')

f = sys.argv[1]
data_col = int(sys.argv[2])
bin_l = float(sys.argv[3])
bin_r = float(sys.argv[4])
bin_w = float(sys.argv[5])

data = []
for l in open(f,'r'):
    if l.find('#') != -1:
        continue
    data.append(float(l.split()[data_col-1]))

b = []
app = bin_l
while app < bin_r:
    b.append(app)
    app += bin_w
b.append(bin_r)

from numpy import histogram

h,b = histogram(data,b)
dens,b = histogram(data,b,density=True)

for (i,h_x) in enumerate(h):
    print((b[i]+b[i+1])*0.5, h_x, dens[i], b[i],b[i+1])
# print average
print('')
print('')
print('#average= %f' % (sum(data)/len(data),))

    

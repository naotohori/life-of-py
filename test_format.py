#!/usr/bin/env python

f = open('test_format.out','w')

#i = 1
#r = 1.23456789e-8
#
#for k in range(20):
#    f.write(f"{i:8d} {r:8.3f}\n")
#    i *= 10
#    r *= 10
#
#f.write("\n")
#
r = 1.23456789e-8

for k in range(20):
    f.write(f"{r:12.6g}\n")
    f.write(f"{-r:12.6g}\n")
    r *= 10

#!/usr/bin/env python

def maintext(Tc):
    return 87.740 - 0.4008 * Tc + 9.398*0.0001*Tc**2 - 1.410 * 0.000001*Tc**3

def abstract(Tc):
    return 87.740 - 0.40008 * Tc + 9.398*0.0001*Tc**2 - 1.410 * 0.000001*Tc**3

for Tc in [0, 10, 25, 50, 75, 100]:
    m = maintext(Tc)
    a = abstract(Tc)
    print (f'{Tc:3d} {a:6.3f} {m:6.3f} {a-m:6.4f} {100*(a-m)/a:6.5f}%')

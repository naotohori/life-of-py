#!/usr/bin/env python
#vim:fileencoding=UTF-8


#球状の乱数分布

import random
import math

random.seed(3)

r=300.0
out_file = open('sphere_rand.dat','w')
out_file.write('%8s %6s %6s %7s %7s\n' 
               % ('#dist','theta','phi','deg(theta)','deg(phi)'))

for i in range(100000):
    x = random.uniform(-r,r)
    y = random.uniform(-r,r)
    z = random.uniform(-r,r)
    dist = math.sqrt(x**2+y**2+z**2)
    theta = math.acos(z/dist)
    phi = math.atan2(y, x)
    if dist > r:
        continue
    out_file.write('%28.22f %26.23f %26.23f %27.22f %27.22f\n' 
                    % (dist,theta,phi,
                       math.degrees(theta), math.degrees(phi)) )
            
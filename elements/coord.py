#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
@author: Naoto Hori
'''

import math
from numpy import dot

class Coord(object) :
    def __init__(self, x=0.0, y=0.0, z=0.0) :
        self.x = x
        self.y = y
        self.z = z

    def get_as_tuple(self):
        return (self.x, self.y, self.z)
    
    def get_as_list(self):
        return [self.x, self.y, self.z]
    
    def put_as_list(self, xyz):
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
    
    def distance(self, another_coord):
        return math.sqrt((another_coord.x - self.x) ** 2 
                       + (another_coord.y - self.y) ** 2
                       + (another_coord.z - self.z) ** 2)

    def move(self, delta_coord):
        self.x += delta_coord.x
        self.y += delta_coord.y
        self.z += delta_coord.z

    def transform(self, mtx): 
        c = [self.x, self.y, self.z, 1.0] 
        c = dot(mtx,c)
        self.x = c[0]
        self.y = c[1]
        self.z = c[2]
    
    def __add__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self
        
    def __div__(self, n):
        '''重心を求める際などに使う'''
        self.x /= n
        self.y /= n
        self.z /= n
        return self

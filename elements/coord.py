#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import math

class Coord(object) :
    def __init__(self, x=0.0, y=0.0, z=0.0) :
        self.x = x
        self.y = y
        self.z = z

    def get_as_tuple(self):
        return (self.x, self.y, self.z)
    
    def get_as_list(self):
        return [self.x, self.y, self.z]
    
    def distance(self, another_coord):
        return math.sqrt((another_coord.x - self.x) ** 2 
                       + (another_coord.y - self.y) ** 2
                       + (another_coord.z - self.z) ** 2)

    def move(self, delta_coord):
        self.x += delta_coord.x
        self.y += delta_coord.y
        self.z += delta_coord.z

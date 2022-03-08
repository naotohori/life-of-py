#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import os
import struct
from lop.elements.error import MyError

        
class CoordFile :
    def __init__(self, filename, nmp) :
        self._filename = filename
        self._nmp = nmp
        
    def open_to_read(self):
        self._file = open(self._filename, 'rb')
        
    def open_to_write(self):
        self._file = open(self._filename, 'wb')
        
    def open_to_addwrite(self):
        self._file = open(self._filename, 'ab')
        
    def close(self):
        self._file.close()
        
    def flush(self):
        self._file.flush()
        
    def read_onestep(self):
        """return 2-dimensional lists"""
        coord_matrix = []
    
        self._file.read(8)  # 0.000
        x = struct.unpack('d' * self._nmp, self._file.read(8*self._nmp))
        self._file.read(8)  # 0.000
        y = struct.unpack('d' * self._nmp, self._file.read(8*self._nmp))
        self._file.read(8)  # 0.000
        z = struct.unpack('d' * self._nmp, self._file.read(8*self._nmp))
        
        for i in range(self._nmp) :
            xyz = [x[i], y[i], z[i]]
            coord_matrix.append(xyz)
        
        return coord_matrix
    
    def skip_onestep(self):
        self._file.seek(3*8*(self._nmp+1), os.SEED_CUR)
     
    def skip(self, num):
        for i in range(num):
            self.skip_onestep()
       
    def write_onestep(self, coord_matrix):
        # for X
        binary = struct.pack('d', 0.0)
        for xyz in coord_matrix :
            binary += struct.pack('d', xyz[0])
        self._file.write(binary)
        # for Y
        binary = struct.pack('d', 0.0)
        for xyz in coord_matrix :
            binary += struct.pack('d', xyz[1])
        self._file.write(binary)
        # for Z
        binary = struct.pack('d', 0.0)
        for xyz in coord_matrix :
            binary += struct.pack('d', xyz[2])
        self._file.write(binary)
    
    def has_more_data(self):
        """return True or False"""
        char = self._file.read(8)
        if not char :
            return False
        else :
            self._file.seek(-8, os.SEEK_CUR)
            return True

    def _read_at(self, num):
        self._file.seek(0)
        for i in range(num - 1) :
            self.read_onestep()
        return self.read_onestep()
        

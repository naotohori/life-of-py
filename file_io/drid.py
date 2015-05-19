#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2015/05/15
@author: Naoto Hori
'''

import os
import struct
from cafysis.elements.error import MyError

class DridHeader(object):
    def __init__(self):
        self.title = None
        self.centroids = []
        self.atoms = []
        self.bonds = []
        
    def show(self):
        print self.title
        print '#(centroid)', len(self.centroids)
        print '#(atom)', len(self.atoms)
        for i in xrange(len(self.centroids)) :
            print 'centroid %i %i' % (i,self.centroids[i])
        for i in xrange(len(self.atoms)) :
            print 'atom %i %i' % (i,self.atoms[i])

class DridFile:
    def __init__(self, filename):
        self._filename = filename
        self._header = DridHeader()
        self._seek_data = None

    def open_to_read(self):
        self._file = open(self._filename, 'rb')

    def open_to_write(self):
        self._file = open(self._filename, 'wb')

    def close(self):
        self._file.close()
        
    def flush(self):
        self._file.flush()
        
    def read_header(self):
        if not self._file :
            raise MyError('DridFile', 'read_header', 'Logical: _file is None')
        
        self._file.seek(0)
        
        # title block (number of line can be changed)
        b = self._pick_data()
        len_title = struct.unpack('i',b)[0]

        b = self._pick_data()
        self._header.title = struct.unpack('%ds' % (len_title,), b)[0]
            
        # centroids 
        b = self._pick_data()
        n_centroids = struct.unpack('i',b)[0]

        b = self._pick_data()
        self._header.centroids = struct.unpack('i'*n_centroids, b)

        # atoms 
        b = self._pick_data()
        n_atoms = struct.unpack('i',b)[0]

        b = self._pick_data()
        self._header.atoms = struct.unpack('i'*n_atoms, b)

        self._seek_data = self._file.tell()
        
    def write_header(self):
        if not self._header :
            raise MyError('DcdFile', 'write_header', 'Logical: _header is None')
        
        self._file.seek(0)
        
        # title block
        binary = struct.pack('i', len(self._header.title))
        self._put_data(binary, 4)

        form = '%ds' % (len(self._header.title),)
        binary = struct.pack(form, self._header.title)
        self._put_data(binary, struct.calcsize(form))
        
        # centroids 
        binary = struct.pack('i', len(self._header.centroids))
        self._put_data(binary, 4)

        binary = struct.pack('i' * len(self._header.centroids), *(self._header.centroids))
        self._put_data(binary, 4 * len(self._header.centroids))

        # atoms
        binary = struct.pack('i', len(self._header.atoms))
        self._put_data(binary, 4)

        binary = struct.pack('i' * len(self._header.atoms), *(self._header.atoms))
        self._put_data(binary, 4 * len(self._header.atoms))
        
            
    def show_header(self):
        """print header information"""
        self._header.show()
        
    def set_header(self, header):
        import copy
        self._header = copy.deepcopy(header)
    
    def get_header(self):
        return self._header
    
    def read_onestep(self):
        """return 2-dimensional lists"""
        b = self._pick_data()
        return struct.unpack('f' * 3 * len(self._header.centroids), b)
    
    def skip_onestep(self):
        num = struct.unpack('i', self._file.read(4))[0]
        self._file.seek(4+num, os.SEEK_CUR)
     
    def skip(self, num):
        for i in xrange(num):
            self.skip_onestep()
       
    def write_onestep(self, data):
        n = len(self._header.centroids)
        binary = struct.pack('%df' % (3*n,), *data)
        self._put_data(binary, 4 * 3 * n)
    
    def has_more_data(self):
        """return True or False"""
        char = self._file.read(4)
        if not char :
            return False
        else :
            self._file.seek(-4, os.SEEK_CUR)
            return True

    def rewind(self):
        self._file.seek( self._seek_data, os.SEEK_SET)

    def _pick_data(self):
        """return binary data between 'integer' and 'integer'. 'integer' indicates the number of bytes"""
        num = struct.unpack('i', self._file.read(4))[0]
        b = self._file.read(num)
        self._file.seek(4, os.SEEK_CUR)
        return b

    def _put_data(self, binary, size):
        self._file.write(struct.pack('<L', size))
        self._file.write(binary)
        self._file.write(struct.pack('<L', size))
        # '<' indicates little endian, 'L' stands 4 byte data
        
    def _read_at(self, num):
        self._file.seek(0)
        self.read_header()
        for i in xrange(num - 1) :
            self.read_onestep()
        return self.read_onestep()
        

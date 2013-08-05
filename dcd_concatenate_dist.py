#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/11/09
@author: Naoto Hori
'''

import sys
import os
import struct
from cafysis.file_io.dcd import DcdFile
from copy import copy

############################
# cafysis.elements.error
############################
class MyError(Exception) :
    def __init__(self, classname='unknown', funcname='unknown', title='unknown'):
        self._class = classname
        self._func = funcname
        self._title = title
        
    def show(self):
        print 'class:', self._class
        print 'function:', self._func
        print 'matter:', self._title

############################
# cafysis.file_io.dcd
############################
class DcdHeader(object):
    def __init__(self):
        self.nset = None
        self.istart = None
        self.nstep_save = None
        self.nstep = None
        self.nunit_real = None
        self.delta = None
        self.title = None
        self.tempk = None
        self.lunit2mp = None
        self.nmp_real = None
        
    def show(self):
        for line in self.title :
            print line
        print 'nset', self.nset
        print 'istart', self.istart
        print 'nstep_save', self.nstep_save
        print 'nstep', self.nstep
        print 'nunit_real', self.nunit_real
        print 'delta', self.delta
        print 'tempk', self.tempk
        for i in xrange(self.nunit_real) :
            print 'lunit2mp[', i, ']', self.lunit2mp[i]
        print 'nmp_real', self.nmp_real
        
class DcdFile :
    def __init__(self, filename) :
        self._filename = filename
        self._header = DcdHeader()
        
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
            raise MyError('DcdFile', 'read_header', 'Logical: _file is None')
        
        self._file.seek(0)
        
        # first block 
        b = self._pick_data()
        bdata = struct.unpack('4siii5iifi9i', b)
        #bdata = struct.unpack('4siii5iid9i',b)
        if bdata[0] != 'CORD' :
            raise MyError('DcdFile', 'read_header', '%s is not coordinate file' % (self._filename,))
        self._header.nset = bdata[1]
        self._header.istart = bdata[2]
        self._header.nstep_save = bdata[3]
        self._header.nstep = bdata[4]
        self._header.nunit_real = bdata[5]
        self._header.delta = bdata[10]
        
        # title block (number of line can be changed)
        b = self._pick_data()
        bdata = struct.unpack(('i' + '80s' * (3 + self._header.nunit_real)), b)
        self._header.title = (bdata[1], bdata[2])
        #self._header.tempk = float(bdata[3].strip('\0 '))
        self._header.tempk = float(bdata[3])
        self._header.lunit2mp = []
        for i in xrange(self._header.nunit_real) :
#            self._header.lunit2mp.append(int(bdata[i + 4].strip('\0 ')))
            self._header.lunit2mp.append(int(bdata[i + 4]))
            
        # nmp_real
        b = self._pick_data()
        self._header.nmp_real = struct.unpack('i', b)[0]
        
    def write_header(self):
        if not self._header :
            raise MyError('DcdFile', 'write_header', 'Logical: _header is None')
        
        self._file.seek(0)
        
        #first block
        format = '4siii5iifi9i'
        binary = struct.pack(format,
                             'CORD',
                             self._header.nset,
                             self._header.istart,
                             self._header.nstep_save,
                             self._header.nstep,
                             self._header.nunit_real,
                             0, 0, 0, 0,
                             self._header.delta,
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        self._put_data(binary, struct.calcsize(format))
        
        #title block
        binary = struct.pack('i' + '80s' * 2,
                             (3 + len(self._header.lunit2mp)),
                             self._header.title[0],
                             self._header.title[1])
        
        import re
        re_null = re.compile('\0')
        
        p = struct.pack('80s', str(self._header.tempk))
        binary += re_null.sub(' ', p)
        for i in xrange(self._header.nunit_real) :
            p = struct.pack('80s', str(self._header.lunit2mp[i]))
            binary += re_null.sub(' ', p)
            
        format = 'i' + '80s' * (3 + self._header.nunit_real)
        self._put_data(binary, struct.calcsize(format))
        
        # nmp_real
        binary = struct.pack('i', self._header.nmp_real)
        self._put_data(binary, 4)
        
            
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
        coord_matrix = []
        b = self._pick_data()
        x = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        y = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        z = struct.unpack('f' * self._header.nmp_real, b)
        
        for i in xrange(self._header.nmp_real) :
            xyz = [x[i], y[i], z[i]]
            coord_matrix.append(xyz)
        
        return coord_matrix
    
    def skip_onestep(self):
        num = struct.unpack('i', self._file.read(4))[0]
        self._file.seek(4+num, os.SEEK_CUR)
        num = struct.unpack('i', self._file.read(4))[0]
        self._file.seek(4+num, os.SEEK_CUR)
        num = struct.unpack('i', self._file.read(4))[0]
        self._file.seek(4+num, os.SEEK_CUR)
     
    def skip(self, num):
        for i in xrange(num):
            self.skip_onestep()
       
    def write_onestep(self, coord_matrix):
        # for X
        binary = ''
        for xyz in coord_matrix :
            binary += struct.pack('f', xyz[0])
        self._put_data(binary, 4 * self._header.nmp_real)
        # for Y
        binary = ''
        for xyz in coord_matrix :
            binary += struct.pack('f', xyz[1])
        self._put_data(binary, 4 * self._header.nmp_real)
        # for Z
        binary = ''
        for xyz in coord_matrix :
            binary += struct.pack('f', xyz[2])
        self._put_data(binary, 4 * self._header.nmp_real)
    
    def has_more_data(self):
        """return True or False"""
        char = self._file.read(4)
        if not char :
            return False
        else :
            self._file.seek(-4, os.SEEK_CUR)
            return True

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
        
############################
# cafysis.dcd_concatenate
############################
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: % SCRIPT [DCD1] [DCD2] ([DCD3] ...) [output DCD]'
        sys.exit(2)
        
    # Number of input DCD files
    num_dcd = len(sys.argv) - 2
    
    # New DCD file
    f_out  = DcdFile(sys.argv[-1])
    f_out.open_to_write()
    
    # Count the total frame number
    num_frame = 1
    for i in xrange(1,num_dcd) :
        f_in = DcdFile(sys.argv[i])
        f_in.open_to_read()
        f_in.read_header()
        num_frame += f_in.get_header().nset - 1
        f_in.close()
        
    # Get the total step number from final DCD file
    f_in = DcdFile(sys.argv[num_dcd])
    f_in.open_to_read()
    f_in.read_header()
    num_frame += f_in.get_header().nset - 1
    num_step = f_in.get_header().nstep
    f_in.close()
        
    f_in = DcdFile(sys.argv[1])
    f_in.open_to_read()
    f_in.read_header() 
    header = copy(f_in.get_header())
    header.nset = num_frame
    header.nstep = num_step
    f_out.set_header(header)
    f_out.write_header()
    print sys.argv[1], f_in.get_header().nset
    while f_in.has_more_data() :
        f_out.write_onestep(f_in.read_onestep())
    f_in.close()
    
    for i in xrange(2,num_dcd+1) :
        f_in = DcdFile(sys.argv[i])
        f_in.open_to_read()
        f_in.read_header()
        f_in.skip_onestep()  # skip the first step
        print sys.argv[i], f_in.get_header().nset - 1
        while f_in.has_more_data() :
            f_out.write_onestep(f_in.read_onestep())
        f_in.close()

    f_out.close()

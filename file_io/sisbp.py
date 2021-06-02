#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import os
import struct
from cafysis.elements.error import MyError

class SisbpFile :
    def __init__(self, filename) :
        self._filename = filename
        self._seek_data = None
        self._seek_mark = None
        self._kind_str = None
        self._kind = None
        
    def open_to_read(self):
        self._file = open(self._filename, 'rb')
        self.read_header()
        
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

        # Read kind
        num = struct.unpack('i', self._file.read(4))[0]
        if num == 2:
            self._kind = 2
            self._kind_str = 'h'
        elif num == 4:
            self._kind = 4
            self._kind_str = 'i'
        elif num == 8:
            self._kind = 8
            self._kind_str = 'q'
        else:
            raise MyError('SisbpFile', 'read_header', 'unknown kind type: '+'%i'%num)

        self._seek_data = self._file.tell()
        
    def read_onestep(self):

        pairs = []

        while True:

            i = self._pick_data()
            if i == 0:
                break

            j = self._pick_data()

            pairs.append((i,j))

        return pairs

    def skip_onestep(self):
        try:
            while(True):
                if (self._pick_data() == 0):
                    break
        except:
            raise EOFError('EOF in a middle of dcd.skip_onestep().')

    ''' This will throw EOFError exception if there are not enough frames.'''
    def skip(self, num):
        for _ in range(num):
            self.skip_onestep()

    ''' This does NOT throw EOFError exception if there are not enough frames.'''
    ''' Instead, this function will return the number of frames skipped successfully.'''
    def skip_as_many_as_possible_upto(self, num):
        for i in range(num):
            try:
                self.skip_onestep()
            except EOFError:
                return i
        return num
            
    def has_more_data(self):
        """return True or False"""
        char = self._file.read(self._kind)
        if not char :
            return False
        else :
            self._file.seek(-self._kind, os.SEEK_CUR)
            return True

    def count_frame(self):
        self.set_mark()

        if self._seek_data is None:
            self.read_header()

        self.rewind()
        
        n = 0
        while self.has_more_data():
            try:
                self.skip_onestep()
            except EOFError:
                break
            n += 1
        
        self.go_mark()

        return n

    def rewind(self):
        self._file.seek( self._seek_data, os.SEEK_SET)

    def set_mark(self):
        self._seek_mark = self._file.tell()

    def go_mark(self):
        self._file.seek( self._seek_mark, os.SEEK_SET)

    def _pick_data(self):
        return struct.unpack(self._kind_str, self._file.read(self._kind))[0]

    def _read_at(self, num):
        self._file.seek(0)
        self.read_header()
        for i in range(num - 1) :
            self.read_onestep()
        return self.read_onestep()
        

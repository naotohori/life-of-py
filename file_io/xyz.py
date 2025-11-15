#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

from lop.elements.coord import Coord
from lop.elements.error import MyError
import numpy as np

class XyzFile(object) :
    def __init__(self, filename, openmode=None):
        self._filename = filename
        self._status = 'Closed'

        if openmode is not None:
            if openmode == 'r':
                self.open_to_read()
            elif openmode == 'w':
                self.open_to_write()
            else:
                print("openmode has to be either 'r' or 'w'")
                raise MyError('PdbFile', '__init__', "openmode has to be either 'r' or 'w'")

    def open_to_read(self):
        if self._status != 'Closed' :
            raise MyError('XyzFile', 'open_for_read', 'file is not closed')
        self._file = open(self._filename, 'r')

    def open_to_write(self):
        if self._status != 'Closed' :
            raise MyError('XyzFile', 'open_for_read', 'file is not closed')
        self._file = open(self._filename, 'w')

    def close(self):
        self._file.close()

    def read(self, np_array=False, close=False):
        seq = []
        coords = []
        l = self._file.readline()
        try:
            n = int(l)
        except:
            raise MyError('XyzFile', 'read', 'could not retrieve the number of coordinates.')

        # Empty line
        l = self._file.readline()

        for i in range(n):
            l = self._file.readline()
            lsp = l.split()
            if len(lsp) != 4:
                raise MyError('XyzFile', 'read', 'Format error :\n' + l)

            seq.append(lsp[0])
            coords.append((float(lsp[1]), float(lsp[2]), float(lsp[3])))

        if close:
            self.close()

        if np_array:
            return seq, np.array(coords)
        else:
            return seq, coords

    def write(self, seq, coords, close=False):

        self._file.write(f'{len(coords)}\n')
        self._file.write('\n')

        for s, xyz in zip(seq, coords):
            self._file.write(f'{s} {xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f}\n')

        if close:
            self.close()

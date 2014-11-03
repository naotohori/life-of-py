#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

from cafysis.elements.error import MyError

class TsHeader(object):
    def __init__(self) :
        self.step = None
        self.temp = None
        self.label = None
        
        self.rg = None
        self.rmsd = None
        self.q = None
        
        # 現状、mgo未対応
        #self.mgo_kai = None
        #self.mgo_coef = None
        #self.mgo_q = None
        #self.mgo_state = None
        
        self.e_tot = None
        self.e_velo = None
        self.e_local = None
        self.e_go = None
        self.e_morse = None
        self.e_repul = None
        self.e_dna_solv = None
        self.e_dna_base = None
        self.e_ele = None
        self.e_ion_hyd = None
        self.e_hp = None
        self.e_lip_tail = None
        self.e_lip_int = None
        self.e_lip_solv = None
        self.e_box = None
        self.e_cap = None
        self.e_bridge = None
        self.e_pull = None
        self.e_anch = None
        self.e_rest1d = None
        self.e_implig = None
        self.e_window = None
        self.e_stack = None
        self.e_hbond = None
        
    def show(self):
        print 'step', self.step
        print 'temp', self.temp
        print 'label', self.label
        
        print 'rg', self.rg
        print 'rmsd', self.rmsd
        print 'q', self.q
        
        #print 'mgo_kai', self.mgo_kai
        #print 'mgo_coef', self.mgo_coef
        #print 'mgo_q', self.mgo_q
        #print 'mgo_state', self.mgo_state
        
        print 'e_tot', self.e_tot
        print 'e_velo', self.e_velo
        print 'e_local', self.e_local
        print 'e_go', self.e_go
        print 'e_morse', self.e_morse
        print 'e_repul', self.e_repul
        print 'e_dna_solv', self.e_dna_solv
        print 'e_dna_base', self.e_dna_base
        print 'e_ele', self.e_ele
        print 'e_ion_hyd', self.e_ion_hyd
        print 'e_hp', self.e_hp
        print 'e_lip_tail', self.e_lip_tail
        print 'e_lip_int', self.e_lip_int
        print 'e_lip_solv', self.e_lip_solv
        print 'e_box', self.e_box
        print 'e_cap', self.e_cap
        print 'e_bridge', self.e_bridge
        print 'e_pull', self.e_pull
        print 'e_anch', self.e_anch
        print 'e_rest1d', self.e_rest1d
        print 'e_implig', self.e_implig
        print 'e_window', self.e_window
        print 'e_stack', self.e_stack
        print 'e_hbond', self.e_hbond
                
class TsFile(object):
    def __init__(self, filename) :
        self._filename = filename
        self.header_lines = []
        self.head_str = None
        self.head_col = None
        self.flg_u_u = False
        self.num_unit = 0
        
    def open_to_read(self):
        self._file = open(self._filename, 'rb')
        
    def open_to_write(self):
        self._file = open(self._filename, 'wb')
        
    def close(self):
        self._file.close()
        
    def flush(self):
        self._file.flush()
        
    def read_header(self):
        ''' store lines to self.header_lines
            store information to self.head_str and self.head_col
        '''
        if not self._file :
            raise MyError('TsFile', 'read_header', 'Logical: _file is None')
        
        self._file.seek(0)
        
        # skip 5 lines
        for i in xrange(5):
            self.header_lines.append( self._file.readline() )
            
        # read header
        self.header_lines.append( self._file.readline() )
        line = self._file.readline()
        self.header_lines.append( line )

        self.head_str = line.split()[1:]
        self.head_col = self._head_str2col(self.head_str)
        
        # skip 2 lines
        self.header_lines.append( self._file.readline() )
        self.header_lines.append( self._file.readline() )
        
        #### Count the number of units ####
        # skip the first line of step 0
        self._file.readline()
        
        if self._file.readline().split()[0] != '#all':
            raise MyError('TsFile','read_header','Bad format: #all')
        
        self.num_unit = 0
        for l in self._file:
            if l[0] == ' ':
                break
            if l.split()[1][0:1] == '-':
                ''' unit-unit lines exist'''
                self.flg_u_u = True
                break
            self.num_unit = int(l.split()[0][1:])
        if self.num_unit < 0:
            raise MyError('TsFile','read_header','Bad format: num_unit <= 0')
        ###########################################
            
        # skip 9 lines
        self._file.seek(0)
        for i in xrange(9):
            self._file.readline()
        
    def write_header(self):
        for l in self.header_lines:
            self._file.write(l)

    def copy_header(self, ts):
        import copy
        self.header_lines = copy.deepcopy(ts.header_lines)
        self.head_str = copy.deepcopy(ts.head_str)
        self.head_col = copy.deepcopy(ts.head_col)
        self.flg_u_u = copy.deepcopy(ts.flg_u_u)
        self.num_unit = copy.deepcopy(ts.num_unit)

    def read_onestep(self):
        lines = []
        lines.append( self._file.readline() )
        
        ts_list = []
        for i in xrange(self.num_unit+1):
            l = self._file.readline()
            lines.append(l)
            ts_list.append(l.split()[1:])

        if self.flg_u_u:
            for i in xrange(self.num_unit*(self.num_unit-1)/2):
                l = self._file.readline()
                lines.append(l)
                ts_list.append(l.split()[1:])
            
        return (ts_list, lines)

    def write_onestep(self, lines):
        for l in lines:
            self._file.write(l)
    
    def skip_onestep(self):
        for i in xrange(self.num_unit +2):
            self._file.readline()
     
    def skip(self, num):
        for i in xrange((self.num_unit +2)*num):
            self._file.readline()
       
    #def write_onestep(self, coord_matrix):
    
    def has_more_data(self):
        """return True or False"""
        s = self._file.tell()
        try:
            if self._file.readline() == '':
                self._file.seek(s)
                return False
            else:
                self._file.seek(s)
                return True
        except:
            self._file.seek(s)
            return False

    #def _pick_data(self):
    #    """return binary data between 'integer' and 'integer'. 'integer' indicates the number of bytes"""
    #    #num = struct.unpack('i', self._file.read(4))[0]
    #    b = self._file.read(num)
    #    self._file.seek(4, os.SEEK_CUR)
    #    return b

    #def _put_data(self, binary, size):
    #    #self._file.write(struct.pack('<L', size))
    #    self._file.write(binary)
    #    #self._file.write(struct.pack('<L', size))
    #    # '<' indicates little endian, 'L' stands 4 byte data
        
    #def _read_at(self, num):
    #    self._file.seek(0)
    #    self.read_header()
    #    for i in xrange(num - 1) :
    #        self.read_onestep()
    #    return self.read_onestep()

    def _head_str2col(self, strlist):
        th = TsHeader()
        for i, str in enumerate(strlist):
            if str == 'step':
                th.step = i
            elif str == 'tempk':
                th.temp = i
            elif str == 'label':
                th.label = i
            elif str == 'radg':
                th.rg = i
            elif str == 'rmsd':
                th.rmsd = i
            elif str == 'qscore':
                th.q = i
            # 現状、mgo未対応
            #elif str == '':
            #    th.mgo_kai = i
            #elif str == '':
            #    th.mgo_coef = i
            #elif str == '':
            #    th.mgo_q = i
            #elif str == '':
            #    th.mgo_state = i
            elif str == 'etot':
                th.e_tot = i
            elif str == 'velet':
                th.e_velo = i
            elif str == 'local':
                th.e_local = i
            elif str == 'go':
                th.e_go = i
            elif str == 'morse':
                th.e_morse = i
            elif str == 'repul':
                th.e_repul = i
            elif str == 'solv_dna':
                th.e_dna_solv = i
            elif str == 'base':
                th.e_dna_base = i
            elif str == 'elect':
                th.e_ele = i
            elif str == 'hyd_ion':
                th.e_ion_hyd = i
            elif str == 'hp':
                th.e_hp = i
            elif str == 'tail':
                th.e_lip_tail = i
            elif str == 'int':
                th.e_lip_int = i
            elif str == 'solv_lip':
                th.e_lip_solv = i
            elif str == 'box':
                th.e_box = i
            elif str == 'cap':
                th.e_cap = i
            elif str == 'bridge':
                th.e_bridge = i
            elif str == 'pulling':
                th.e_pull = i
            elif str == 'anchor':
                th.e_anch = i
            elif str == 'rest1d':
                th.e_rest1d = i
            elif str == 'imp_lig':
                th.e_implig = i
            elif str == 'window':
                th.e_window = i
            elif str == 'stack':
                th.e_stack = i
            elif str == 'hbond':
                th.e_hbond = i
            else:
                raise MyError('','','')
        return th

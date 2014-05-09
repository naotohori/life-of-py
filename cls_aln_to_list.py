#!/usr/bin/env python

        
class Aln :
    def __init__(self):
        self.names = []
        self.seq = []
#        self.seq_idnt = []
        self.length = 0
        self.idnt2loc = []
        self.loc2idnt = []
        self.loc2res = []
        self.res2loc = []
        self.gap = []
        
from my_error import Error

class AlnFile :
    def __init__(self, filename, status):
        if status in ('r','w') :
            self._file = open(filename, status)
        else :
            raise Error('PdbFile', 'open_for_read', 'file is not closed')
        
    def close(self):
        self._file.close()
        
    def read_all(self):
        al = Aln()
        
        line = self._file.readline()
#        flg_read_aln = False
        while line :
            if line[6:27] == '|PDBID|CHAIN|SEQUENCE':
                self._add_aln(al, line)
#                flg_read_aln = True
#            elif flg_read_aln :
#                self._add_idnt(al, line)
#                flg_read_aln = False
            line = self._file.readline()
                
        length = len(al.seq[0])
        for idx in xrange(len(al.names)) :
            if len(al.seq[idx]) != length :
                raise Error('AlnFile', 'read_all', 'al.seq[idx] != length')
            
            al.loc2res.append([])
            al.res2loc.append([-1])
            res_num = 0
            # res -- 1 start
            # loc -- 0 start
            for i,s in enumerate(al.seq[idx]) :
                if s == '-' :
                    al.loc2res[idx].append(-1)
                else :
                    res_num += 1
                    al.loc2res[idx].append(res_num)
                    al.res2loc[idx].append(i)
        al.length = length
#                    
#        for i,s in enumerate(al.seq_idnt) :
#            if i == length :
#                break
#            if s == ' ' :
#                al.loc2idnt.append(False)
#            elif s == '*' :
#                al.loc2idnt.append(True)
#                al.idnt2loc.append(i)
#            else :
#                raise Error('AlnFile', 'read_all', 'unknown letter in seq_idnt')
        
        for i in xrange(length) :
            ref = al.seq[0][i]
            flg = True 
            for idx in xrange(len(al.names)) :
                if al.seq[idx][i] != ref :
                    flg = False
                    break
            al.loc2idnt.append(flg)
                
        for i in xrange(length) :
            flg = True
            for idx in xrange(len(al.names)) :
                if al.seq[idx][i] == '-' :
                    flg = False
                    break
            if not flg :
                al.gap.append(True)
            else :
                al.gap.append(False)
                        
        return al
            
    def _add_aln(self, al, line):
        if line[6:27] != '|PDBID|CHAIN|SEQUENCE':
            raise Error('AlnFile', '_add_aln', 'line[6:27] != |PDBID|CHAIN|SEQUENCE')
        name = line[0:6]
        seqall = line[33:83]
        if seqall.find(' ') != -1:
            seq = seqall[0: seqall.find(' ')]
        else :
            seq = seqall[:]
        
#        print name
#        print seq
        
        if name in al.names :
            idx = al.names.index(name)
        else :
            al.names.append(name)
            idx = len(al.names) - 1
        
        if idx < len(al.seq) :
            al.seq[idx].extend(list(seq))
        else :
            al.seq.append(list(seq))
        
#    def _add_idnt(self, al, line):
#        al.seq_idnt.extend(list(line[33:83]))

import sys

if len(sys.argv) != 3:
    print ('\n Usage: [input aln file] [output list file]\n')
    sys.exit(2)
    
f_aln = AlnFile(sys.argv[1], 'r')
f_out = open(sys.argv[2], 'w')

al = f_aln.read_all()

for i in xrange(al.length) :
    print al.seq[0][i], al.seq[1][i], al.loc2res[0][i], al.loc2res[1][i], al.loc2idnt[i], al.gap[i]
print 'al.loc2res'
print al.loc2res[0]
print 'al.res2loc'
print al.res2loc[0]
print 'loc2idnt'
print al.loc2idnt
print 'idnt2loc'
print al.idnt2loc

# write header
f_out.write('#')
for name in al.names :
    f_out.write('%6s' % name)
f_out.write('\n')

# write data
for i in xrange(al.length) :
    if al.gap[i] :
        continue
    f_out.write(' ')
    for idx in xrange(len(al.names)) :
        f_out.write(' %5i' % al.loc2res[idx][i])
    f_out.write('\n')



#### Aln file sample
'''
0         1         2         3         4         5         6         7         8         9
0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
3R8O_A|PDBID|CHAIN|SEQUENCE      AAUUGAAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAA
3I8H_A|PDBID|CHAIN|SEQUENCE      ---UGGAGAGUUUGAUCCUGGCUCAGGGUGAACGCUGGCGGCGUGCCUAA
                                    ** *********** ********  **************  ******
 
3R8O_A|PDBID|CHAIN|SEQUENCE      CACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGA
3I8H_A|PDBID|CHAIN|SEQUENCE      GACAUGCAAGUCGUGCGGG---CCGCGGGGUUUUACUCCGU----GGUCA
                                  ************  ***      *  *    ** ** * *    *   *
'''

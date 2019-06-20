#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile

def count(path):

    dcd = DcdFile(path)
    dcd.open_to_read()

    dcd.read_header()

    icount = 0
    while dcd.has_more_data():
        try:
            dcd.skip(1)
            icount += 1
        except EOFError:
            icount += 1
            print ('There is another frame at the end, that is not written completely.')
            break
        except:
            break

    dcd.close()

    return icount


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print ('Usage: SCRIPT [dcd file]')
        sys.exit(2)

    print(('# frames = %i' % (count(sys.argv[1]))))

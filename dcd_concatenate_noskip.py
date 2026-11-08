#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/11/09
@author: Naoto Hori

Thin wrapper kept for backward compatibility.
Equivalent to: dcd_concatenate.py --noskip [DCD1] [DCD2] ... [output DCD]
'''

import sys
from lop.dcd_concatenate import dcd_concatenate


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: % SCRIPT [DCD1] [DCD2] ([DCD3] ...) [output DCD]')
        print('(Equivalent to: dcd_concatenate.py --noskip ...)')
        sys.exit(2)
    n_frames = dcd_concatenate(sys.argv[1:], skip_first=False)
    for i,n in enumerate(n_frames):
        print('file ',i+1,': ', n)

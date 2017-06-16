#!/usr/bin/env python

import sys
import argparse
from CalcRMSD import calcrmsd
from cafysis.file_io.dcd import DcdFile
from cafysis.dcd_frame_count import count

parser = argparse.ArgumentParser(description='Compute auto-correation of RMSD from a dcd trajectory',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#parser.add_argument('--gap', dest='gap', default=1,
#                    action='store', type=int, help='Gap')
parser.add_argument('--skip', dest='skip', default=0,
                    action='store', type=int, help='number of frames to be skipped')
parser.add_argument('--first_mp', dest='first_mp', default=0,
                    action='store', type=int, help='first particle ID')
parser.add_argument('--last_mp', dest='last_mp', default=-1,
                    action='store', type=int, help='last particle ID')
parser.add_argument('--logscale', dest='flg_log', default=False,
                    action='store_true', help='log scale time')

parser.add_argument('dcd', help='target DCD file')
parser.add_argument('out', help='output file')

args = parser.parse_args()

################################################################################

n_frame = count(args.dcd)   # e.g. 501   (i=0, 1, 2,..., 500)
dcd = DcdFile(args.dcd)
dcd.open_to_read()
dcd.read_header()

f_out = open(args.out,'w')

sum_rmsd = [0.0] * (n_frame-args.skip+1)
count_rmsd = [0.0] * (n_frame-args.skip+1)

iframe_max = n_frame - 1  # 500 (= 501 - 1)

for iframe_ref in range(args.skip, iframe_max):  # iframe_ref = 5, 6, ...., 499

    dcd.rewind()

    dcd.skip( iframe_ref )
    data_ref = dcd.read_onestep_npF()
    iframe = iframe_ref
    
    for interval in range(1, iframe_max-iframe_ref+1):   # interval = 1, 2, 3, ... , 495 (= 501-5-1)
        '''
        inteval =   1:   iframe_ref = 5, 6, ... , 499
        inteval = 495:   iframe_ref = 5
        '''
        print iframe_ref, interval, iframe
        data = dcd.read_onestep_npF()
        iframe += 1

        rmsd = calcrmsd(data, data_ref)
        print iframe_ref, iframe, interval, rmsd

        sum_rmsd[ interval ] += rmsd
        count_rmsd[ interval ] += 1

dcd.close()

for interval in range(1, n_frame-args.skip):
    f_out.write('%i %f %i\n' % (interval, sum_rmsd[interval]/float(count_rmsd[interval]), count_rmsd[interval]))

f_out.close()


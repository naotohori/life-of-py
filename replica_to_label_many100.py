#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
import argparse
from lop.file_io.dcd import DcdFile
from lop.file_io.ts import TsFile
#import resource

NUM_OUTFILE_OPEN = 100

#resource.setrlimit(resource.RLIMIT_NOFILE, (2000,2000))

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Usage: %SCRIPT [input Dir] [name] [final ID] [output DIR]')
        print('       ID shoule be from 1 to final ID.')
        sys.exit(2)

    ###########################################
    ## Define the Parser and check errors    ##
    ###########################################
    parser = argparse.ArgumentParser(
             description='Convert REMD simulation results',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('indir', help='Input directory')
    parser.add_argument('nrep', help='Number of replicas')
    parser.add_argument('outdir',help='Output directory')

    group_seq = parser.add_mutually_exclusive_group(required=True)
    group_seq.add_argument('--pdb', type=PdbFile, help='PDB file')
    group_seq.add_argument('--seq', help='Sequence')
    group_seq.add_argument('--seqfile', type=argparse.FileType('r'), help='Sequence FASTA file')

    parser.add_argument('--hbfile', type=argparse.FileType('r'), help='HB list file')

    parser.add_argument('--circ', action="store_true", help='Flag for circRNA')
    
    parser.add_argument('--end5', default='S', help="5'-end P or S")

    parser.add_argument('--exvfile', type=argparse.FileType('w'), help='output exv pair file')

    args = parser.parse_args()

    if args.pdb is not None:
        #pdb = PdbFile(sys.argv[1])
        args.pdb.open_to_read()
        chains = args.pdb.read_all()
        args.pdb.close()
        
        chain = chains[0]
        n_nt = chain.num_res()
        
        seq = []
        for r in chain.residues:
            # "RA " ---> "A"
            seq.append(r.atoms[0].res_name.strip()[1])

    elif args.seq is not None:
        seq = args.seq
        n_nt = len(seq)

    elif args.seqfile is not None:
        seq = ''
        for l in args.seqfile:
            if l[0] == '>' or l[0] == '#' or l[0] == ';':
                continue
            seq += l.strip()
        n_nt = len(seq)

    else:
        print ('Error: this error should be detected by mutually_exclusive_group of the parser')
        sys.exit(2)

    if args.circ:
        if args.end5 == 'S':
            print ('end5 option is ignored for circRNA')
            args.end5 = 'P'

    if args.end5 == 'P':
        offset_end5 = 1
    elif args.end5 == 'S':
        offset_end5 = 0
    else:
        print ('end5 option has to be Either S or P')
        sys.exit(2)


    dir_in = sys.argv[1]
    name = sys.argv[2]
    id_end = int(sys.argv[3])
    dir_out = sys.argv[-1]

    #Open input files
    in_ts_files = []
    in_dcd_files = []
    for id_rep in range(1, id_end+1):
        ts = TsFile('%s/%s_%04i.ts' % (dir_in, name, id_rep))
        ts.open_to_read()
        #ts.read_header()   ## Commented out since the header will be read later
        in_ts_files.append(ts)

        dcd = DcdFile('%s/%s_%04i.dcd' % (dir_in, name, id_rep))
        dcd.open_to_read()
        #dcd.read_header()  ## Commented out since the header will be read later
        in_dcd_files.append(dcd)

    id_finish = 0
    while id_finish < id_end:
        
        id_begin_now = id_finish + 1
        id_end_now = id_finish + NUM_OUTFILE_OPEN
        if id_end_now > id_end:
            id_end_now = id_end
        id_finish = id_end_now

        # Prepare output files
        list_id_out = list(range(id_begin_now, id_end_now+1))
        out_ts_files = {}
        out_dcd_files = {}
        for id_lab in list_id_out:
            tsout = TsFile('%s/%s_%04i.ts' % (dir_out, name, id_lab))
            tsout.open_to_write()
            out_ts_files[id_lab] = tsout
            dcdout = DcdFile('%s/%s_%04i.dcd' % (dir_out, name, id_lab))
            dcdout.open_to_write()
            out_dcd_files[id_lab] = dcdout

        # Rewind (by reading header again)
        for id_rep in range(1, id_end+1):
            idx = id_rep - 1

            in_ts_files[idx].read_header()
            in_dcd_files[idx].read_header()
 
            if (id_rep in list_id_out):
                out_ts_files[id_rep].copy_header(in_ts_files[idx])
                out_ts_files[id_rep].write_header()
                out_dcd_files[id_rep].set_header( in_dcd_files[idx].get_header() )
                out_dcd_files[id_rep].write_header()
    
        while in_dcd_files[0].has_more_data():
            for idx in range(id_end):
                (tsdata, tslines) = in_ts_files[idx].read_onestep()
                coord_matrix = in_dcd_files[idx].read_onestep()
    
                step  = int( tsdata[0][in_ts_files[idx].head_col.step] )
                label = int( tsdata[0][in_ts_files[idx].head_col.label] )

                if step == 1:
                    (tsdata, tslines) = in_ts_files[idx].read_onestep()
                    label = int( tsdata[0][in_ts_files[idx].head_col.label] )

                if (label in list_id_out):
                    out_ts_files[label].write_onestep( tslines )
                    out_dcd_files[label].write_onestep( coord_matrix )


        # Close output files
        for f in list(out_ts_files.values()):
            f.close()
        for f in list(out_dcd_files.values()):
            f.close()
    

    # Close input files
    for f in in_ts_files:
        f.close()
    for f in in_dcd_files:
        f.close()

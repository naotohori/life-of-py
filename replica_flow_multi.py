#!/usr/bin/env python

import sys
import os
import argparse

def read_replica_table(f):

    flg_started = False
    ndim = None
    nrep = 0
    table = {}

    for l in open(f):
        if l.startswith('# Table of replica variable'):
            flg_started = True
            continue

        if not flg_started:
            continue

        if l.startswith('#'):
            continue

        if flg_started and len(l.strip()) == 0:
            break

        lsp = l.split()

        if ndim is None:
            ndim = len(lsp) - 1

        irep = int(lsp[0])
        table[irep] = list(map(float, lsp[1:ndim+1]))
        if irep > nrep:
            nrep = irep

    ilabel = [1]*ndim
    var_label = [[1,]]*ndim
    previous = table[1][:]
    nrep_var = [0]*ndim

    for irep in range(2, nrep+1):
        for idim in range(ndim):
            if table[irep][idim] > previous[idim]:
                ilabel[idim] += 1
                previous[idim] = table[irep][idim]
            elif table[irep][idim] < previous[idim]:
                ilabel[idim] = 1
                previous[idim] = table[irep][idim]
            var_label[idim].append(ilabel[idim])

            if ilabel[idim] > nrep_var[idim]:
                nrep_var[idim] = ilabel[idim]

    #print(len(var_label))
    #print(var_label[0])
    #print(var_label[1])
    #print(var_label[2])
    #print(nrep_var)

    return ndim, nrep, nrep_var, var_label


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
             description='Measure replica flows based on .rep files',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', nargs='+', help='input .rep file(s)')
    parser.add_argument('output', default='replica_flow.out', help='output filename')
    parser.add_argument('--verbose', default=False, action="store_true", help="print out some additional info.")
    parser.add_argument('--overwrite', default=False, action="store_true", help="overwrite the output file if exists.")

    args = parser.parse_args()


    ndim, nrep, nrep_var, var_label = read_replica_table(args.input[0])

    if args.verbose:
        print('ndim = ', ndim)
        print('nrep = ', nrep)
        print('nrep_var = ', nrep_var)

    if args.overwrite:
        f_out = open(args.output, 'w')
    else:
        if not os.path.isfile(args.output):
            f_out = open(args.output, 'w')
        else:
            print('Error: the output file exsits. Please delete it first.')
            sys.exit(2)

    #filenames = sys.argv[2:len(sys.argv)-1]

    #flows = [0]*(nrep+1)
    #up = [0]*(nrep+1)
    #up_steps = [0]*(nrep+1)
    #down = [0]*(nrep+1)
    #down_steps = [0]*(nrep+1)

    #last_visit = [0]*(nrep+1)
    #memory = [0]*(nrep+1)
    #exchange = [0]*nrep
    #steps = 0

    flows = [[0]*(nrep+1)]*ndim
    up = [[0]*(nrep+1)]*ndim
    up_steps = [[0]*(nrep+1)]*ndim
    down = [[0]*(nrep+1)]*ndim
    down_steps = [[0]*(nrep+1)]*ndim

    last_visit = [[0]*(nrep+1)]*ndim
    memory = [[0]*(nrep+1)]*ndim
    #exchange = []
    #for idim in range(ndim):
    #    exchange.append([0]*nrep_var[idim])

    steps = 0

    for ifile, filename in enumerate(args.input):

        flg_read = False
        #flg_skip = False

        for l in open(filename):

            if l.startswith('# History'):
                flg_read = True
                #if ifile > 0:
                #    # To skip the first step
                #    flg_skip = True
                continue

            if not flg_read:
                continue

            #if flg_skip:
            #    # To skip the first step
            #    flg_skip = False
            #    continue

            if len(l) < 3:
                break

            ## The step column will be labels[0] which won't be used.
            labels = list(map(int, l.split()))
            steps += 1

            if steps == 1:
                for irep, lbl in enumerate(labels):
                    if irep == 0:
                        continue
                    for idim in range(ndim):
                        l = var_label[idim][lbl-1]
                        #memory[idim][irep] = l

                        if l == 1:
                            flows[idim][irep] = 1
                            last_visit[idim][irep] = steps
                        elif l == nrep_var[idim]:
                            flows[idim][irep] = -1
                            last_visit[idim][irep] = steps

            else:
                for irep, lbl in enumerate(labels):
                    if irep == 0:
                        continue

                    for idim in range(ndim):
                        l = var_label[idim][lbl-1]
                        #if memory[idim][irep] == l - 1:
                        #    exchange[idim][l-1] += 1
                        #memory[idim][irep] = l

                        if l == 1:
                            if flows[idim][irep] == -1:
                                down[idim][irep] += 1
                                down_steps[idim][irep] += steps - last_visit[idim][irep]

                            flows[idim][irep] = 1
                            last_visit[idim][irep] = steps

                        elif l == nrep_var[idim]:
                            if flows[idim][irep] == 1:
                                up[idim][irep] += 1
                                up_steps[idim][irep] += steps - last_visit[idim][irep]

                            flows[idim][irep] = -1
                            last_visit[idim][irep] = steps


#f_out.write('#Exchange probabilities\n')
#for lbl in range(1, Nrep):
#    f_out.write('%5i %5i %20i %8.4f\n' % (lbl, lbl+1, exchange[lbl], exchange[lbl] / (0.5 * float(steps))))
#
#f_out.write('\n\n')

f_out.write('# Flow\n')
for idim in range(ndim):
    f_out.write(f'# dimension {idim+1}\n')

    for irep in range(1, nrep+1):
        f_out.write('%4i'  % irep)

        f_out.write(' %10i' % up[idim][irep])
        if up[idim][irep] > 0:
            f_out.write(' %10.1f'  % (up_steps[idim][irep]/float(up[idim][irep]),))
        else:
            f_out.write(' Inf')

        f_out.write(' %10i' % down[idim][irep])
        if down[idim][irep] > 0:
            f_out.write(' %10.1f'  % (down_steps[idim][irep]/float(down[idim][irep]),))
        else:
            f_out.write(' Inf')

        f_out.write(' %10i' % (up[idim][irep] + down[idim][irep],))
        if up[idim][irep] + down[idim][irep] > 0:
            f_out.write(' %10.1f\n'  % ((up_steps[idim][irep]+down_steps[idim][irep])/float(up[idim][irep]+down[idim][irep]),))
        else:
            f_out.write(' Inf\n')
    f_out.write('\n')


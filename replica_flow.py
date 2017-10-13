#!/usr/bin/env python

import sys

if len(sys.argv) < 3:
    print 'Usage: SCRIPT [rep file1] [[rep file2] [rep file3] ...] [output]'
    sys.exit(2)

filenames = sys.argv[1:len(sys.argv)-1]
f_out = open(sys.argv[-1], 'w')

Nrep = 64

flows = [0]*(Nrep+1)
travels = [0]*(Nrep+1)
travel_times = [0]*(Nrep+1)
last_visit = [0]*(Nrep+1)

memory = [0]*(Nrep+1)
exchange = [0]*Nrep
steps = 0


for ifile, filename in enumerate(filenames):

    flg_read = False
    flg_skip = False

    for l in open(filename):

        if l[0:9] == '# History':
            flg_read = True
            if ifile > 0:
                # To skip the first step
                flg_skip = True
            continue

        if not flg_read:
            continue

        if flg_skip:
            # To skip the first step
            flg_skip = False
            continue

        if len(l) < 3:
            break

        ## The step column will be labels[0] which won't be used.
        labels = map(int, l.split())
        steps += 1

        if steps == 1:
            for irep, lbl in enumerate(labels):
                if irep == 0:
                    continue
                memory[irep] = lbl

                if lbl == 1:
                    flows[irep] = 1
                    last_visit[irep] = steps
                elif lbl == Nrep:
                    flows[irep] = -1
                    last_visit[irep] = steps

        else:
            for irep, lbl in enumerate(labels):

                if memory[irep] == lbl - 1:
                    exchange[lbl-1] += 1
                memory[irep] = lbl
                
                if lbl == 1:
                    if flows[irep] == -1:
                        travels[irep] += 1
                        travel_times[irep] += steps - last_visit[irep]

                    flows[irep] = 1
                    last_visit[irep] = steps

                elif lbl == Nrep:
                    if flows[irep] == 1:
                        travels[irep] += 1
                        travel_times[irep] += steps - last_visit[irep]

                    flows[irep] = -1
                    last_visit[irep] = steps


f_out.write('#Exchange probabilities\n')
for lbl in range(1, Nrep):
    f_out.write('%i %i %i %f\n' % (lbl, lbl+1, exchange[lbl], exchange[lbl] / (0.5 * float(steps))))

f_out.write('\n\n')
f_out.write('# Flow\n')
for irep in range(1, Nrep+1):
    if travels[irep] > 0:
        f_out.write('%i %i %f\n'  % (irep, travels[irep], travel_times[irep]/travels[irep]))
    else:
        f_out.write('%i %i Inf\n'  % (irep, travels[irep]))

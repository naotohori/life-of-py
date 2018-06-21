#!/usr/bin/env python

import sys

if len(sys.argv) < 4:
    print 'Usage: SCRIPT [#replica] [rep file1] [[rep file2] [rep file3] ...] [output]'
    sys.exit(2)

filenames = sys.argv[2:len(sys.argv)-1]
f_out = open(sys.argv[-1], 'x')

Nrep = int(sys.argv[1])

flows = [0]*(Nrep+1)
up = [0]*(Nrep+1)
up_steps = [0]*(Nrep+1)
down = [0]*(Nrep+1)
down_steps = [0]*(Nrep+1)

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
                        down[irep] += 1
                        down_steps[irep] += steps - last_visit[irep]

                    flows[irep] = 1
                    last_visit[irep] = steps

                elif lbl == Nrep:
                    if flows[irep] == 1:
                        up[irep] += 1
                        up_steps[irep] += steps - last_visit[irep]

                    flows[irep] = -1
                    last_visit[irep] = steps


f_out.write('#Exchange probabilities\n')
for lbl in range(1, Nrep):
    f_out.write('%5i %5i %20i %8.4f\n' % (lbl, lbl+1, exchange[lbl], exchange[lbl] / (0.5 * float(steps))))

f_out.write('\n\n')
f_out.write('# Flow\n')
for irep in range(1, Nrep+1):
    f_out.write('%4i'  % irep)

    f_out.write(' %10i' % up[irep])
    if up[irep] > 0:
        f_out.write(' %10.1f'  % (up_steps[irep]/float(up[irep]),))
    else:
        f_out.write(' Inf')

    f_out.write(' %10i' % down[irep])
    if down[irep] > 0:
        f_out.write(' %10.1f'  % (down_steps[irep]/float(down[irep]),))
    else:
        f_out.write(' Inf')

    f_out.write(' %10i' % (up[irep] + down[irep],))
    if up[irep] + down[irep] > 0:
        f_out.write(' %10.1f\n'  % ((up_steps[irep]+down_steps[irep])/float(up[irep]+down[irep]),))
    else:
        f_out.write(' Inf\n')


#!/usr/bin/env python

def read_bpseq(bpseq_file):
    with open(bpseq_file, 'r') as f:
        lines = f.readlines()

    sequence = []
    pairs = []
    for line in lines:
        if line.startswith('#'):
            continue
        pos, nt, partner = line.strip().split()
        pos = int(pos)
        partner = int(partner)
        sequence.append(nt)
        pairs.append(partner)

    return sequence, pairs

def write_ct(ct_file, sequence, pairs):
    with open(ct_file, 'w') as f:
        f.write(f'{len(sequence)}\tRNA structure\n')
        for i, (nt, partner) in enumerate(zip(sequence, pairs), start=1):
            prev_pos = i - 1
            next_pos = i + 1 if i < len(sequence) else 0
            f.write(f'{i}\t{nt}\t{prev_pos}\t{next_pos}\t{partner}\t{i}\n')

def convert_bpseq_to_ct(bpseq_file, ct_file):
    sequence, pairs = read_bpseq(bpseq_file)
    write_ct(ct_file, sequence, pairs)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print ('Usage: [bpseq file] [ct file]')
        sys.exit(0)

    bpseq_file = sys.argv[1]
    ct_file = sys.argv[2]
    convert_bpseq_to_ct(bpseq_file, ct_file)


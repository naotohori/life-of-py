#!/usr/bin/env python

def find_max_nt_in_pairs(pairs):
    Nnt = 0
    for pair in pairs:
        if Nnt < max(pair):
            Nnt = max(pair)
    return Nnt

def pairs2Mountain(pairs, Nnt=None, scaled=True):

    """ If Nnt is not given, find the largest nucleotide ID in the pairs."""
    if Nnt is None:
        Nnt = find_max_nt_in_pairs(pairs)

    f = [0.]*Nnt
    for (nt1, nt2) in pairs:
        if nt2 < nt1:
            nt1, nt2 = nt2, nt1
        if scaled:
            l = nt2 - nt1
            f[nt1-1] += 1./float(l)
            f[nt2-1] += -1./float(l)
        else:
            f[nt1-1] += 1
            f[nt2-1] += -1

    s = 0.
    for i in range(Nnt):
        s += f[i]
        f[i] = s

    return f

def Mountain_distance(pairs1, pairs2, Nnt=None, p=1, scaled=True):

    """ If Nnt is not given, find the largest nucleotide ID in the pairs."""
    if Nnt is None:
        Nnt = find_max_nt_in_pairs(pairs1 + pairs2)

    m1 = pairs2Mountain(pairs1, Nnt, scaled=scaled)
    m2 = pairs2Mountain(pairs2, Nnt, scaled=scaled)

    d = 0.0
    for i in range(Nnt):
        d += (abs(m1[i] - m2[i]))**p

    return d**(1./float(p))


if __name__ == "__main__":

    import sys
    from rna_ss_convert import parse_input_file

    if len(sys.argv) != 3:
        print ("Usage: SCRIPT [input 1 (.db/.ct/.bpseq)] [input 2 (.db/.ct/.bpseq)]")
        sys.exit(2)
    else:
        seq1, pairs1 = parse_input_file(sys.argv[1])
        seq2, pairs2 = parse_input_file(sys.argv[2])

        print(Mountain_distance(pairs1, pairs2))


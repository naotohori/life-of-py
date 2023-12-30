#!/usr/bin/env python

from lop.file_io.sisbp import SisbpFile
import os
import sys
import argparse

def interpret_dot_bracket(dot_bracket):
    bras = {1: '(', 2: '[', 3: '{', 4: '<'}
    kets = {')': 1, ']': 2, '}': 3, '>': 4}
    stack = {'(': [], '[': [], '{': [], '<': []}
    pairs = []

    for i, char in enumerate(dot_bracket):
        if char in bras.values():
            stack[char].append(i)
        elif char in kets.keys():
            opening = bras[kets[char]]
            if len(stack[opening]) == 0:
                raise ValueError("Unbalanced brackets in the input string.")
            opening_index = stack[opening].pop()
            pairs.append((opening_index+1, i+1))

    for s in stack.values():
        if len(s) != 0:
            raise ValueError("Unbalanced brackets in the input string.")

    return pairs


def parse_input_file(input_file):
    extension = os.path.splitext(input_file)[1].lower()

    if extension == '.db':
        with open(input_file, 'r') as f:
            sequence = f.readline().strip()
            dot_bracket = f.readline().strip()
            pairs = interpret_dot_bracket(dot_bracket)
    elif extension in ['.ct', '.bpseq']:
        pairs = []
        sequence = ''

        with open(input_file, 'r') as f:
            if extension == '.ct':
                f.readline()  # Skip the header line for ct format

            for line in f:
                fields = line.strip().split()
                i, nt, partner = int(fields[0]), fields[1], int(fields[4])

                sequence += nt
                if extension == '.ct' and partner > i:
                    pairs.append((i, partner))  # ID starting from 1
                elif extension == '.bpseq' and partner != 0 and partner > i:
                    pairs.append((i, partner))  # ID starting from 1

    return sequence, pairs


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
             description='Calculate number of base pairs from bp file',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bpfile', type=SisbpFile, help='input bp file')
    parser.add_argument('outfile', default='nbp.out', help='output filename')

    parser.add_argument('--ss', help='Secondary structure (.db/.bpseq/.ct)')
    parser.add_argument('--ene', type=float, help='Energy threshold')

    args = parser.parse_args()

    if args.ene is None:
        args.ene = 0.0
    else:
        if args.ene > 0.0:
            print ('Warning: a positive value was specified for --ene. Note that BP energies are normally negative.')

    native_pairs = None
    sequence = None
    n_native = 0
    if args.ss is not None:
        sequence, native_pairs = parse_input_file(args.ss)
        n_native = len(native_pairs)

    f = args.bpfile

    f.open_to_read()
    fout = open(args.outfile, 'w')

    while f.has_more_data():

        pairs, energies = f.read_onestep()

        n = 0
        if native_pairs is not None:
            for pair, e in zip(pairs, energies):
                if pair in native_pairs:
                    if e < args.ene:
                        n += 1

            fout.write(f"{n} {n/n_native:6.4f}\n")
        else:
            for pair, e in zip(pairs, energies):
                if e < args.ene:
                    n += 1
            fout.write(f"{n}\n")

    f.close()

#!/usr/bin/env python

from lop.file_io.sisbp import SisbpFile
import os
import sys
import argparse

def interpret_dot_bracket(dot_bracket):
    stack = []
    brackets = {'(': ')', '[': ']', '{': '}', '<': '>'}
    pairs = []

    for i, char in enumerate(dot_bracket):
        if char in brackets.keys():
            stack.append((char, i))
        elif char in brackets.values():
            if not stack:
                raise ValueError("Unbalanced brackets in the input string.")
            opening_bracket, opening_index = stack.pop()
            if brackets[opening_bracket] != char:
                raise ValueError("Mismatched brackets in the input string.")
            pairs.append((opening_index+1, i+1))  # ID starting from 1

    if stack:
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
                i, nt, partner = int(fields[0]), fields[1], int(fields[-1])

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

    args = parser.parse_args()

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
            for pair in pairs:
                if pair in native_pairs:
                    n += 1

            fout.write(f"{n} {n/n_native:6.4f}\n")
        else:
            fout.write(f"{n}\n")

    f.close()

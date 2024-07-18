#!/usr/bin/env python

import os

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
            pairs.append((opening_index, i))

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
                if extension == '.ct':
                    i, nt, partner = int(fields[0]), fields[1], int(fields[4])
                else:
                    i, nt, partner = int(fields[0]), fields[1], int(fields[-1])

                sequence += nt
                if extension == '.ct' and partner > i:
                    pairs.append((i - 1, partner - 1))
                elif extension == '.bpseq' and partner != 0 and partner > i:
                    pairs.append((i - 1, partner - 1))

    return sequence, pairs

def write_output_file(sequence, pairs, output_file):
    extension = os.path.splitext(output_file)[1].lower()

    if extension == '.db':
        dot_bracket = ['.'] * len(sequence)
        for i, j in pairs:
            dot_bracket[i] = '('
            dot_bracket[j] = ')'

        with open(output_file, 'w') as f:
            f.write(sequence + '\n')
            f.write(''.join(dot_bracket) + '\n')

    elif extension == '.ct':
        with open(output_file, 'w') as f:
            f.write(f"{len(sequence)}\tConverted_sequence\n")
            for i, nt in enumerate(sequence):
                partner = 0
                for pair in pairs:
                    if i in pair:
                        partner = pair[0] if pair[1] == i else pair[1]
                        partner += 1
                        break

                f.write(f"{i + 1}\t{nt}\t{i}\t{i + 2}\t{partner}\t{i + 1}\n")

    elif extension == '.bpseq':
        with open(output_file, 'w') as f:
            for i, nt in enumerate(sequence):
                partner = 0
                for pair in pairs:
                    if i in pair:
                        partner = pair[0] if pair[1] == i else pair[1]
                        partner += 1
                        break

                f.write(f"{i + 1}\t{nt}\t{partner}\n")

def convert_format(input_file, output_file):
    sequence, pairs = parse_input_file(input_file)
    write_output_file(sequence, pairs, output_file)

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print ("Usage: SCRIPT [input (.db/.ct/.bpseq)] [output (.db/.ct/.bpseq)]")
        sys.exit(2)
    else:
        convert_format(sys.argv[1], sys.argv[2])

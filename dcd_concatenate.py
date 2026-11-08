#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/11/09
@author: Naoto Hori
'''

import os
import sys
import struct
import argparse
from lop.file_io.dcd import DcdFile
from copy import copy


CHUNK_SIZE = 64 * 1024 * 1024  # bytes per read/write in the raw copy


def _scan_input(path):
    '''Read the header of an input DCD and return its layout.

    Frames are copied as raw bytes, so every frame must have the same size.
    The number of complete frames is derived from the file size; an
    incomplete trailing frame is reported and ignored.'''
    f = DcdFile(path)
    f.open_to_read()
    try:
        f.read_header()
        header = f.get_header()
        seek_data = f._seek_data
        bpf = f._bpf
        data_size = os.path.getsize(path) - seek_data
        n_frame, remainder = divmod(data_size, bpf)

        # Sanity check of the record markers of the first and the last complete frames
        first_record = 48 if header.with_unit_cell else 4 * header.nmp_real
        for iframe in ({0, n_frame - 1} if n_frame > 0 else ()):
            f._file.seek(seek_data + iframe * bpf)
            marker = struct.unpack('i', f._file.read(4))[0]
            if marker != first_record:
                raise ValueError('%s: unexpected record size %d at frame %d (expected %d). '
                                 'Frames are not of uniform size.'
                                 % (path, marker, iframe + 1, first_record))
    finally:
        f.close()

    if remainder != 0:
        print('Warning: %s has an incomplete frame at the end (%d bytes), which is ignored.'
              % (path, remainder), file=sys.stderr)

    return {'header': header, 'seek_data': seek_data, 'bpf': bpf, 'n_frame': n_frame}


def _copy_bytes(fin, fout, nbytes):
    while nbytes > 0:
        buf = fin.read(min(CHUNK_SIZE, nbytes))
        if not buf:
            raise EOFError('Unexpected end of file while copying')
        fout.write(buf)
        nbytes -= len(buf)


def dcd_concatenate(filepaths, skip_first=True):
    '''Concatenate DCD files. filepaths = [input1, input2, ..., output].
    If skip_first is True, the first frame of the second and later inputs
    is skipped (it duplicates the last frame of the previous run on restart).
    Frame data (coordinates and unit cell) are copied as raw bytes.
    Returns the list of the number of frames taken from each input.'''

    inputs = filepaths[:-1]
    output = filepaths[-1]

    # The output is opened with 'wb', so never let it be one of the inputs
    for path in inputs:
        if os.path.exists(output) and os.path.samefile(path, output):
            raise ValueError('Output %s is the same file as an input' % (output,))

    # Check all inputs and count the frames
    layouts = [_scan_input(path) for path in inputs]
    h0 = layouts[0]['header']
    for path, lay in zip(inputs[1:], layouts[1:]):
        h = lay['header']
        if h.nmp_real != h0.nmp_real:
            raise ValueError('%s has %d atoms, but %s has %d'
                             % (path, h.nmp_real, inputs[0], h0.nmp_real))
        if h.with_unit_cell != h0.with_unit_cell:
            raise ValueError('%s and %s differ in the presence of unit cell information'
                             % (path, inputs[0]))

    n_skip = [0] + [1 if skip_first else 0] * (len(inputs) - 1)
    n_frames = [max(lay['n_frame'] - s, 0) for lay, s in zip(layouts, n_skip)]

    # Header: taken from the first input, total frames, and nstep of the last input
    header = copy(h0)
    header.nset = sum(n_frames)
    header.nstep = layouts[-1]['header'].nstep

    f_out = DcdFile(output)
    f_out.open_to_write()
    try:
        f_out.set_header(header)
        f_out.write_header()

        for path, lay, s, n in zip(inputs, layouts, n_skip, n_frames):
            with open(path, 'rb') as fin:
                fin.seek(lay['seek_data'] + s * lay['bpf'])
                _copy_bytes(fin, f_out._file, n * lay['bpf'])
    finally:
        f_out.close()

    return n_frames


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Concatenate DCD files. By default, the first frame of '
                    'the second and later inputs is skipped, as it duplicates '
                    'the last frame of the previous run (restart).')
    parser.add_argument('inputs', nargs='+', help='input DCD files (two or more)')
    parser.add_argument('output', help='output DCD file')
    parser.add_argument('--noskip', action='store_true',
                        help='keep the first frame of every input DCD')
    args = parser.parse_args()
    if len(args.inputs) < 2:
        parser.error('at least two input DCD files are required')
    n_frames = dcd_concatenate(args.inputs + [args.output],
                               skip_first=not args.noskip)
    for i,n in enumerate(n_frames):
        print('file ',i+1,': ', n)

#!/usr/bin/env python

import pickle
import pathlib
import math
import argparse

import alignment.conv as conv


def getrotation(line):
    x, y, z, phi = map(float, line.split()[15:19])
    return conv.ev2mat(phi, [x, y, z])


def main(pklfile, protdir, outfile):
    with open(pklfile, 'rb') as fh:
        guesses = pickle.load(fh)
    results = {}
    for guess in guesses:
        results[guess[0]] = guess[1:]

    protdir = pathlib.Path(protdir)
    for protf in protdir.iterdir():
        protid = protf.stem
        with protf.open() as f:
            lines = f.readlines()
        for lineno, line in enumerate(lines):
            key = (protid, lineno+1)
            if key in results:
                guessrot = getrotation(line)
                truerot = conv.ev2mat(results[key][4][0],
                                      results[key][4][1:])
                diff = conv.so3dist(guessrot, truerot)
                results[key] = (*results[key],
                                guessrot, truerot, diff)

    output = []
    for k in results:
        t = '{}\t{}\t{}'.format(k[0], k[1], results[k][-1])
        output.append(t)

    with open(outfile, 'w') as fh:
        fh.write('\n'.join(output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pklfile',
                        help='.pkl file with alignment results.')
    parser.add_argument('protdir',
                        help='Directory containing test '
                        'protein files.')
    parser.add_argument('outfile',
                        help='Output file')
    args = parser.parse_args()
    main(args.pklfile, args.protdir, args.outfile)

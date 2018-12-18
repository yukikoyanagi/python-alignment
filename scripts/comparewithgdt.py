#!/usr/bin/env python3

from argparse import ArgumentParser
from operator import itemgetter
from collections import OrderedDict
import pathlib

def loaddata(filepath):
    data = {}
    fpath = pathlib.Path(filepath).resolve()
    with open(fpath) as fh:
        for line in fh:
            cols = line.split()
            target = cols[0].split('_')[0]
            try:
                data[target].append((cols[0], float(cols[1])))
            except KeyError:
                data[target] = [(cols[0], float(cols[1]))]
    return data


def output(data):
    text = []
    for target in data:
        text.append('Target: {}'.format(target))
        subdata = data[target]
        text.append('\tHighest GDT score: '
                    '\n\t\t{}, {}'.format(subdata[0][0],
                                          subdata[0][1]))
        text.append('\tHighest alignment score:')
        for cand in subdata[1:]:
            text.append('\t\t{}, {}, GDT={}, '
                        'diff={:.2f}, rank={}'.format(*cand))
        text.append('\n')
    print('\n'.join(text))

        
def compare(gdt, other, howmany=1):
    #Ensure gdt and other share a common target set
    target1 = set(gdt.keys())
    target2 = set(other.keys())
    for key in target1 - target2:
        gdt.pop(key)
    for key in target2 - target1:
        other.pop(key)

    results = {}
    for target in gdt:
        first = sorted(gdt[target], key=itemgetter(1), reverse=True)
        firstod = OrderedDict(first)
        second = sorted(other[target], key=itemgetter(1), reverse=True)
        best1 = first[0]
        data = [(best1[0], best1[1], best1[1], 0 , 1)]
        i = 0
        for cand in second:
            candname = cand[0]
            candscore = cand[1]
            try:
                gdtscore = firstod[candname]
                diff = best1[1] - firstod[candname]
                rank = first.index((candname, gdtscore)) + 1
            except KeyError:
                continue
            data.append((
                candname, candscore, gdtscore, diff, rank
            ))
            i += 1
            if i == howmany:
                break
        results[target] = data

    return results


def run(gdtfile, otherfile, howmany=1):
    data = compare(loaddata(gdtfile), loaddata(otherfile), howmany)
    output(data)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('gdtfile', help='GDT score file')
    parser.add_argument('otherfile',
                        help='File to be compared to GDT file')
    parser.add_argument('-m', '--howmany', type=int, default=1,
                        help='How many resuts from the other file '
                        'to display')
    args = parser.parse_args()
    run(args.gdtfile, args.otherfile, args.howmany)

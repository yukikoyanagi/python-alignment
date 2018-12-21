#!/usr/bin/env python

import argparse
import pickle
from itertools import islice, tee
from operator import itemgetter

import pp

import alignment.sequenceloader
import alignment.alignmenthelper
import alignment.substmatrices


def chunk(data, n=24):
    '''Splits input data (dict) into n equal chunks.'''
    it = iter(data)
    for i, iterator in enumerate(tee(it, n)):
        yield {k: data[k] for k in islice(iterator, i, None, n)}


def findbestalignment(target, candidates, targetid, window,
                      matrix, gap):
    '''Find the candidate that best aligns with the target'''
    lpsl = alignment.sequenceloader.localpatternsequenceloader()
    tseq = lpsl.load(target, window)
    cscores = []
    for candidate in candidates:
        cseq = lpsl.load(candidate, window)
        helper = alignment.alignmenthelper.AlignmentHelper()
        score = helper.alignStrict(tseq, cseq,
                                   substMatrix=matrix, gapScore=gap,
                                   btrace=False)
        cscores.append(
            (candidate, cseq, score, candidates[candidate]))
    cscores = sorted(cscores, key=operator.itemgetter(2), reverse=True)
    top = cscores[0]
    return tseq, top[0], top[1], top[3]


def alignmany(tpats, patrots, window, matrix, gap):
    '''Run alignment sequentially'''
    results = []
    for k in tpats:
        tseq, cpat, cseq, crot = findbestalignment(
            tpats[k], patrots, k, window, matrix, gap)
        results.append((k, tpats[k], tseq, cpat, cseq, crot))
    return results


def loadpatrot(patrotfile):
    patrots = {}
    with open(patrotfile) as fh:
        for line in fh:
            pat, rot = line.split()
            p, v = rot.split(';')
            x, y, z = v.split(',')
            p, x, y, z = map(float, [p, x, y, z])
            patrots[pat] = (p, x, y, z)
    return patrots

def main(targetfile, candidatesfile, window, matrix, gap, outfile):
    #Load targetfile
    with open(targetfile, 'rb') as fh:
        tpats = pickle.load(fh)
    patrots = loadpatrot(candidatesfile)

    matrices = alignment.substmatrices.SubstitutionMatrices()
    matrix = matrices.getmatrix(matrix)

    #Set up ppservers
    try:
        ppservers = open('/tmp/nodelist').read().strip().split()
    except IOError:
        # It wasn't called by slrm bashscript. Run locally.
        jobserver = pp.Server()
    else:
        ppservers = tuple(pp + ':2048' for pp in ppservers)
        jobserver = pp.server(ppservers=ppservers)
    ncpus = jobserver.get_ncpus()

    jobs = []
    for subset in chunk(tpats, ncpus):
        job = jobserver.submit(
            alignmany,
            args=(subset, patrots, window, matrix, gap),
            depfuncs=(findbestalignment,),
            modules=('operator',
                     'alignment.sequenceloader',
                     'alignment.alignmenthelper'))
        jobs.append(job)

    jobserver.wait()

    results = []
    for job in jobs:
        results.extend(job())
    with open(outfile, 'wb') as fh:
        pickle.dump(results, fh)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('targetfile',
                        help='.pkl file containing target protein '
                        'patterns.')
    parser.add_argument('candidatesfile',
                        help='File containing pattern and rotation '
                        'information, extracted from _assess file.')
    parser.add_argument('outfile',
                        help='Output .pkl file')
    parser.add_argument('--window', type=int,
                        help='Window size used for pattern generation')
    parser.add_argument('--matrix',
                        help='Substitution matrix.',
                        default='bondlength')
    parser.add_argument('--gap', type=float,
                        default=-1.0,
                        help='Score for each gap')
    args = parser.parse_args()

    main(args.targetfile, args.candidatesfile, args.window,
         args.matrix, args.gap, args.outfile)
    

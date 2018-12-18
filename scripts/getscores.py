#!/usr/bin/env python

import argparse
import pathlib

import pp

import alignment.sequenceloader
import alignment.alignmenthelper
import alignment.substmatrices

SUBST_MATRIX = alignment.substmatrices.simplematrix
DEFAULT_GAP_SCORE = -1


def alignmany(target, cands, length,
              match=None, mismatch=None, mat=None, gap=-1):
    """Return alignment score for many candidate structures."""
    """Requires pathlib module"""
    csl = alignment.sequenceloader.clustersequenceloader()
    helper = alignment.alignmenthelper.AlignmentHelper()
    tseq = csl.load(target, length)
    targetname = pathlib.Path(target).stem
    results = {}
    for cand in cands:
        cname = pathlib.Path(cand).stem
        cseq = csl.load(cand, length)
        score = helper.alignStrict(tseq, cseq,
                           matchScore=match, mismatchScore=mismatch,
                           substMatrix=mat, gapScore=gap,
                           btrace=False)
        results[cname] = score
    return results


def run(tdir, cdir, seqf,
        match=None, mismatch=None, mat=None,
        gap=DEFAULT_GAP_SCORE, output=None):
    targetdir = pathlib.Path(tdir)
    canddir = pathlib.Path(cdir)
    if match is None and mat is None:
        mat=SUBST_MATRIX
    if gap is None:
        gap=DEFAULT_GAP_SCORE
    files = {}
    
    for candfile in canddir.iterdir():
        try:
            files[candfile.stem.split('_')[0]].append(candfile)
        except KeyError:
            files[candfile.stem.split('_')[0]] = [candfile]

    seqlen = {}
    with open(seqf) as fh:
        for line in fh:
            target, length, __ = line.split()
            seqlen[target] = int(length)

    # Set up ppservers
    try:
        ppservers = open('tmp/nodelist').read().strip().split()
        ppservers = tuple(pp + ':2048' for pp in ppservers)
        jobserver = pp.server(ppservers=ppservers)
    except IOError:
        # It wasn't called by slrm bashscript. Run locally.
        jobserver = pp.Server()
    ncpus = jobserver.get_ncpus()
    if len(files) < ncpus:
        jobserver.set_ncpus(len(files))

    jobs = []
    for targetfile in targetdir.iterdir():
        targetname = targetfile.stem
        try:
            cands = [str(f.resolve()) for f in files[targetname]]
        except KeyError:
            continue
        job = jobserver.submit(alignmany,
                               args=(str(targetfile.resolve()),
                                     cands,
                                     seqlen[targetfile.stem],
                                     match,
                                     mismatch,
                                     mat,
                                     gap),
                               modules=('pathlib',
                                        'alignment.sequenceloader',
                                        'alignment.alignmenthelper',
                                        'alignment.substmatrices'))
        jobs.append(job)

    jobserver.wait()

    out = []
    for job in jobs:
        result = job()
        out += ['{}\t{}'.format(k, result[k]) for k in result]
    if output is None:
        print('\n'.join(out))
    else:
        outfile = pathlib.Path(output)
        outfile.write_text('\n'.join(out))
        

if __name__ == '__main__':
    matrices = alignment.substmatrices.SubstitutionMatrices()
    alist = matrices.available
    available = ['{}: {}'.format(k, alist[k]) for k in alist]
    available = ['Available substitution matrices:'] + available
    parser = argparse.ArgumentParser(
        epilog='\n  '.join(available),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('targetdir',
                        help='Directory containing target files')
    parser.add_argument('canddir',
                        help='Directory containing candidate files')
    parser.add_argument('seqf',
                        help='File with length of each target')
    parser.add_argument('--match', type=float,
                        help='Score for each match')
    parser.add_argument('--mismatch', type=float,
                        help='Score for each mismatch')
    parser.add_argument('--matrix', type=str,
                        help='Substitution matrix to use')
    parser.add_argument('--gap', type=float,
                        help='Score for each gap')
    parser.add_argument('--output',
                        help='Output file')
    args = parser.parse_args()
    if args.matrix:
        matrices = alignment.substmatrices.SubstitutionMatrices()
        mat = matrices.getmatrix(args.matrix)
    else:
        mat = None

    run(args.targetdir, args.canddir, args.seqf,
        args.match, args.mismatch, mat, args.gap, args.output)

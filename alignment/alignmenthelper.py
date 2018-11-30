from .sequence import Sequence
from .vocabulary import Vocabulary
from . import sequencealigner as seqal

class AlignmentHelper(object):
    '''Helper class to simplify alignment.'''
    def __init__(self):
        pass

    def alignStrict(self, first, second,
                    matchScore=None, mismatchScore=None,
                    substMatrix=None, gapScore=-1, btrace=True):
        voc = Vocabulary()
        fe = voc.encodeSequence(Sequence(first))
        se = voc.encodeSequence(Sequence(second))
        if matchScore is not None and mismatchScore is not None:
            scoring = seqal.SimpleScoring(matchScore, mismatchScore)
        elif substMatrix is not None:
            scoring = seqal.MatrixScoring(
                voc.encodeScoreMatrix(substMatrix))
        else:
            raise TypeError('No score provided.')
        aligner = seqal.StrictGlobalSequenceAligner(
            scoring, gapScore
        )
        if btrace:
            score, alignments = aligner.align(fe, se, backtrace=btrace)
            return score, [voc.decodeSequenceAlignment(g)
                           for g in alignments]
        else:
            score = aligner.align(fe, se, backtrace=btrace)
            return score

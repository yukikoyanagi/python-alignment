from .vocabulary import Vocabulary
from .sequence import Sequence
from .sequencealigner import MatrixScoring
from .sequencealigner import GlobalSequenceAligner
from .sequencealigner import StrictGlobalSequenceAligner
from .sequencealigner import LocalSequenceAligner
from .alignmenthelper import AlignmentHelper


DEFAULT_SUBST_MATRIX = {
            ('a', 'a'): 3, ('a', 'b'): -1, ('a', 'c'): -1,
            ('b', 'b'): 1, ('b', 'c'): -1,
            ('c', 'c'): 1
        }
DEFAULT_MATCH_SCORE = 3
DEFAULT_MISMATCH_SCORE = -1
DEFAULT_GAP_SCORE = -2


def alignstrict(first, second, matrix):
    voc = Vocabulary()
    fe = voc.encodeSequence(Sequence(first))
    se = voc.encodeSequence(Sequence(second))
    submat = voc.encodeScoreMatrix(matrix)
    aligner = StrictGlobalSequenceAligner(MatrixScoring(submat),
                                          DEFAULT_GAP_SCORE)
    score, alignments = aligner.align(fe, se, backtrace=True)
    return score, [voc.decodeSequenceAlignment(g) for g in alignments]

def test_substitution_matrix():
    s1 = 'abc'
    s2 = 'bac'
    score, alignments = alignstrict(s1, s2, DEFAULT_SUBST_MATRIX)

    assert len(alignments) == 1
    assert str(alignments[0].first) == '- a b c'
    assert str(alignments[0].second) == 'b a - c'        
    assert score == 4 + DEFAULT_GAP_SCORE *2

    score, alignments = alignstrict(s2, s1, DEFAULT_SUBST_MATRIX)

    assert str(alignments[0].second) == '- a b c'
    assert str(alignments[0].first) == 'b a - c'        
    assert score == 4 + DEFAULT_GAP_SCORE * 2

def test_alignmenthelper():
    helper = AlignmentHelper()
    score, alignments = helper.alignStrict(
        'abc', 'bac', substMatrix=DEFAULT_SUBST_MATRIX,
        gapScore=DEFAULT_GAP_SCORE)

    assert len(alignments) == 1
    assert str(alignments[0].first) == '- a b c'
    assert str(alignments[0].second) == 'b a - c'        
    assert score == 4 + DEFAULT_GAP_SCORE *2

    score, alignments = helper.alignStrict(
        'xaby', 'ab',
        matchScore=DEFAULT_MATCH_SCORE,
        mismatchScore=DEFAULT_MISMATCH_SCORE,
        gapScore=DEFAULT_GAP_SCORE
    )
    assert len(alignments) == 1
    assert str(alignments[0].first) == 'x a b y'
    assert str(alignments[0].second) == '- a b -'
    assert alignments[0].score == score
    assert alignments[0].percentIdentity() == 2.0 / 4.0 * 100.0
    assert alignments[0].percentSimilarity() == 2.0 / 4.0 * 100.0
    assert alignments[0].percentGap() == 2.0 / 4.0 * 100.0
    assert score == DEFAULT_MATCH_SCORE * 2 + DEFAULT_GAP_SCORE * 2

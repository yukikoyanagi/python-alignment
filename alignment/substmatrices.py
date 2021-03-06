# File: substmatrices.py
#
# Description: List of substitution matrices
#
# Author: Yuki Koyanagi
#

from itertools import combinations_with_replacement as cwr
from pkg_resources import resource_filename
import pickle
import math


import alignment.conv as conv


class SubstitutionMatrices(object):
    alphabet_cluster = [
        'F2a', 'F2x',
        'F3a', 'F3b', 'F3c', 'F3x',
        'F4a', 'F4x',
        'F5a', 'F5x',
        'F6a', 'F6x',
        'LRa', 'LRb', 'LRc', 'LRd', 'LRe', 'LRf', 'LRx',
        'NC-',
        'R2a', 'R2b', 'R2c', 'R2x',
        'R3a', 'R3b', 'R3c', 'R3d', 'R3e', 'R3x',
        'R4a', 'R4b', 'R4x',
        'R5a', 'R5b', 'R5c', 'R5d', 'R5e', 'R5x',
        'R6a', 'R6b', 'R6x',
        'UNB'
    ]
    
    alphabet_length = [
        'U', 'A', 'L+', 'L-',
        '+2+', '+2-', '-2+', '-2-',
        '+3+', '+3-', '-3+', '-3-',
        '+4+', '+4-', '-4+', '-4-',
        '+5+', '+5-', '-5+', '-5-',
        '+6+', '+6-', '-6+', '-6-',
    ]

    alphabet_protein = list('ARNDCQEGHILKMFPSTWYVBJZXOU*')

    available = {'simplematrix': 'Simple matrix with match=1 '
                 'and mismatch=-1',
                 'clusterdistance': 'Mismatch scores are '
                 'given by -(distance between clusters).',
                 'bondlength': 'Scores based on hydrogen bond '
                 'lengths.',
                 'bondlength2': 'Similar to bondlength, but '
                 'short-to-short substitution is scaled '
                 'exponentially.',
                 'bondlength3': 'Similar to bondlength, but '
                 'short-to-short substitution is scaled '
                 'logarithmically.',
                 'bondlength4': 'Similar to bondlength, but '
                 'twist has a larger penalty (-0.8).',
                 'blosum62': 'blosum62 substitution matrix.',
                 'blosum80': 'blosum80 substitution matrix.',
                 'blosum90': 'blosum90 substitution matrix.',
                 'betapairs': 'blosum-like substitution matrix '
                 'for beta-strand pairing.',
                 'blosum62plus': 'blosum62, plus @ for helix '
                 'and = for beta strands.'}

    def __init__(self, alphabet=None):
        if alphabet is None:
            alphabet = self.alphabet_cluster
        elif alphabet == 'cluster':
            alphabet = self.alphabet_cluster
        elif alphabet == 'length':
            alphabet = self.alphabet_length
        self.alphabet = alphabet

    @property
    def simplematrix(self):
        a2 = cwr(self.alphabet, 2)
        matrix = {k: (-1)**int(k[0] != k[1]) for k in a2}
        return matrix

    @property
    def clusterdistance(self):
        fname = resource_filename('alignment',
                                  'config/clustermode.pkl')
        with open(fname, 'rb') as f:
            modes = pickle.load(f)
        a2 = cwr(self.alphabet_cluster, 2)
        matrix = {}
        for k in a2:
            if k[0] == k[1]:
                dist = -1  # this gives match score = 1
            else:
                try:
                    first = modes[k[0]]
                    second = modes[k[1]]
                    dist = conv.so3dist(
                        conv.ev2mat(first[0], first[1]),
                        conv.ev2mat(second[0], second[1])
                    )
                except KeyError:
                    # The given key is not found in the dict of modes.
                    # Most likely because the key is of the type '??x'
                    # or 'NC-'/'UNB'
                    dist = 3.14
            matrix[k] = -1 * dist
        return matrix

    def loadproteinmatrix(self, fname):
        with open(fname) as fh:
            lines = fh.readlines()
        rmat = [l for l in lines if l[0] != '#']
        ks = rmat[0].strip().split()
        matrix = {}
        for row in rmat[1:]:
            cs = row.strip().split()
            thiskey = cs[0]
            for k, v in zip(ks, cs[1:]):
                matrix[(thiskey, k)] = int(v)
        # key '*' corresponds to O or U in actual sequences
        append = {}
        for kpair in matrix.keys():
            if kpair[0] == '*' and kpair[1] == '*':
                rks = ['OO', 'OU', 'UO', 'UU']
                for rk in rks:
                    append[tuple(rk)] = matrix[kpair]
            elif kpair[0] == '*':
                append[('O', kpair[1])] = matrix[kpair]
                append[('U', kpair[1])] = matrix[kpair]
            elif kpair[1] == '*':
                append[(kpair[0], 'O')] = matrix[kpair]
                append[(kpair[0], 'U')] = matrix[kpair]
        matrix.update(append)
        return matrix
        
    @property
    def blosum62(self):
        fname = resource_filename('alignment',
                                  'config/blosum62.txt')
        return self.loadproteinmatrix(fname)

    @property
    def blosum62plus(self):
        fname = resource_filename('alignment',
                                  'config/blosum62.txt')
        mat = self.loadproteinmatrix(fname)
        append = {}
        append[('@', '@')] = 4
        for k in self.alphabet_protein:
            append[(k, '@')] = -4
        append[('=', '=')] = 4
        append[('@', '=')] = -4
        for k in self.alphabet_protein:
            append[(k, '=')] = -4
        mat.update(append)
        return mat

    def blosum62plus2(self, match=4, mismatch=-4):
        fname = resource_filename('alignment',
                                  'config/blosum62.txt')
        mat = self.loadproteinmatrix(fname)
        append = {}
        append[('@', '@')] = match
        for k in self.alphabet_protein:
            append[(k, '@')] = mismatch
        append[('=', '=')] = match
        append[('@', '=')] = mismatch
        for k in self.alphabet_protein:
            append[(k, '=')] = mismatch
        mat.update(append)
        return mat
    
    @property
    def blosum80(self):
        fname = resource_filename('alignment',
                                  'config/blosum80.txt')
        return self.loadproteinmatrix(fname)
    
    @property
    def blosum90(self):
        fname = resource_filename('alignment',
                                  'config/blosum90.txt')
        return self.loadproteinmatrix(fname)

    @property
    def betapairs(self):
        fname = resource_filename('alignment',
                                  'config/bstrands.txt')
        return self.loadproteinmatrix(fname)

    @property
    def bondlength(self):
        a2 = cwr(self.alphabet_length, 2)
        matrix = {}
        for k in a2:
            if k[0] == k[1]:
                score = 1
            elif k[0] == 'U' or k[0] == 'A':
                score = -1
            elif k[0][:-1] == k[1][:-1]:
                # Same lengths
                score = 0
            elif k[0][:-1] == 'L':
                # Long to short substitution
                score = -0.75
            else:
                # Short to short substitution
                d = abs(int(k[0][:-1]) - int(k[1][:-1]))
                score = -d/20
            # Twist penalty
            if len(k[0]) > 1 and k[0][-1] != k[1][-1]:
                score += -0.1
            matrix[k] = score
        return matrix

    @property
    def bondlength2(self):
        a2 = cwr(self.alphabet_length, 2)
        matrix = {}
        for k in a2:
            if k[0] == k[1]:
                score = 1
            elif k[0] == 'U' or k[0] == 'A':
                score = -1
            elif k[0][:-1] == k[1][:-1]:
                # Same lengths
                score = 0
            elif k[0][:-1] == 'L':
                # Long to short substitution
                score = -0.75
            else:
                # Short to short substitution
                d = abs(int(k[0][:-1]) - int(k[1][:-1]))
                score = math.expm1(d) / math.expm1(12) * -0.6
            # Twist penalty
            if len(k[0]) > 1 and k[0][-1] != k[1][-1]:
                score += -0.1
            matrix[k] = score
        return matrix

    @property
    def bondlength3(self):
        a2 = cwr(self.alphabet_length, 2)
        matrix = {}
        for k in a2:
            if k[0] == k[1]:
                score = 1
            elif k[0] == 'U' or k[0] == 'A':
                score = -1
            elif k[0][:-1] == k[1][:-1]:
                # Same lengths
                score = 0
            elif k[0][:-1] == 'L':
                # Long to short substitution
                score = -0.75
            else:
                # Short to short substitution
                d = abs(int(k[0][:-1]) - int(k[1][:-1]))
                score = math.log1p(d) / math.log1p(12) * -0.6
            # Twist penalty
            if len(k[0]) > 1 and k[0][-1] != k[1][-1]:
                score += -0.1
            matrix[k] = score
        return matrix

    @property
    def bondlength4(self):
        a2 = cwr(self.alphabet_length, 2)
        matrix = {}
        for k in a2:
            if k[0] == k[1]:
                score = 1
            elif k[0] == 'U' or k[0] == 'A':
                score = -1
            elif k[0][:-1] == k[1][:-1]:
                # Same lengths
                score = 0
            elif k[0][:-1] == 'L':
                # Long to short substitution
                score = -0.75
            else:
                # Short to short substitution
                d = abs(int(k[0][:-1]) - int(k[1][:-1]))
                score = -d/20
            # Twist penalty
            if len(k[0]) > 1 and k[0][-1] != k[1][-1]:
                score += -0.8
            matrix[k] = score
        return matrix

    def getmatrix(self, matrixname):
        if matrixname == 'simplematrix':
            return self.simplematrix
        elif matrixname == 'clusterdistance':
            return self.clusterdistance
        elif matrixname == 'bondlength':
            return self.bondlength
        elif matrixname == 'bondlength2':
            return self.bondlength2
        elif matrixname == 'bondlength3':
            return self.bondlength3
        elif matrixname == 'bondlength4':
            return self.bondlength4
        elif matrixname == 'blosum62':
            return self.blosum62
        elif matrixname == 'blosum80':
            return self.blosum80
        elif matrixname == 'blosum90':
            return self.blosum90
        elif matrixname == 'betapairs':
            return self.betapairs
        elif matrixname == 'blosum62plus':
            return self.blosum62plus
        else:
            raise KeyError('{} not found.'.format(matrixname))

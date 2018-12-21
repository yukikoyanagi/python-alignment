=========
Alignment
=========

Alignment is a native Python library for generic sequence alignment. It is
useful in cases where your alphabet is arbitrarily large and you cannot use
traditional biological sequence analysis tools. It supports global and local
pairwise sequence alignment. I also plan to add support for profile-profile
alignments, but who knows when.

Fork
====

The original alignment library was forked to create package specifically for
the local pattern/GDT analyses. It allows the use of substitution matrices
and class:: alignmenthelper class was created to simpligy alignment.

For local pattern alignment, sequence is based on the local pattern around
a hydrogen bond. Each letter in a sequence represents a peptide unit, i.e.
a sequence of O, N, C-alpha atoms. If there are partial peptide units at
either end of the pattern window they are also represented by a letter.
Each peptide unit is assinged a letter according to the following rule;
- If the peptide unit's N atom is a donor, assign a string representing
  the length (and the direction) of the bond, together with whether the
  bond is twisted. This is done by either 2- or 3-letter string, with the
  first one or two letters representing the length (+2, -4, L, etc) and the
  last letter indicating twistedness (+: twisted or -;not twisted).
- Of the remaining peptide units, assign the letter A to any whose O atom
  is a acceptor of a hydrogen bond.
- All remaining peptide units are assigned the letter U.

If the local pattern is disconnected, the two segments are concatenated to
form a single word.

Each sequence associated to a hydrogen bond in the test protein dataset is
then aligned with all sequences associated with cluster results. Then the
rotation value associated with sequence with the highest alignment score
is chosed as the prediction for a given hydrogen bond.

Some scripts are provided for the analysis. The two relevant for the local
pattern analysis are alignlocpats.py and computediff.py. They are installed
by pip. alignlocpats takes two input data; one is the local pattern for
the test protein data, in pickle format. The data should be a dict of the
form {(proteinname, lineno): localpattern}. The other is the cluster result
text file, extracted from an _assess file. This should contain two columns,
first column with the local pattern text, and the second with the rotation.
slrm- bash script is also provided (but not installed) to be used on abacus
(or any other slurm-based) cluster.


Installation
============

You can install the most recent release using pip:

    pip install alignment

Usage
=====

Typical usage looks like this::

    from alignment.sequence import Sequence
    from alignment.vocabulary import Vocabulary
    from alignment.sequencealigner import SimpleScoring, GlobalSequenceAligner

    # Create sequences to be aligned.
    a = Sequence('what a beautiful day'.split())
    b = Sequence('what a disappointingly bad day'.split())

    # Create a vocabulary and encode the sequences.
    v = Vocabulary()
    aEncoded = v.encodeSequence(a)
    bEncoded = v.encodeSequence(b)

    # Create a scoring and align the sequences using global aligner.
    scoring = SimpleScoring(2, -1)
    aligner = GlobalSequenceAligner(scoring, -2)
    score, encodeds = aligner.align(aEncoded, bEncoded, backtrace=True)

    # Iterate over optimal alignments and print them.
    for encoded in encodeds:
        alignment = v.decodeSequenceAlignment(encoded)
        print alignment
        print 'Alignment score:', alignment.score
        print 'Percent identity:', alignment.percentIdentity()
        print

TODO List
=========

- Profile-profile alignment is not working yet.

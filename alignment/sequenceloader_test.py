from . import sequenceloader

class Testlocalpatternsequenceloader(object):

    _loader = sequenceloader.localpatternsequenceloader()
    
    def test_findwindowsize(self):
        seqs = [[0,1,2,3,4],[100,101,102,103,104]]
        assert self._loader.findwindowsize(seqs,  None, None) == 2
        seqs = [list(range(12))]
        edges = [(9,2,'-'), (3,8,'+')]
        assert self._loader.findwindowsize(seqs, edges, None) == 3
        seqs = [list(range(17))]
        edges = [(4,12,'-'), (10,6,'-'), (13,9,'-')]
        length = -2
        assert self._loader.findwindowsize(seqs, edges, length) == 4

    def test_connectseqs(self):
        seqs = [[0,1,2,3,4],list(range(100,105))]
        edges = [(2,102,'-')]
        seqs, edges = self._loader.connectseqs(
            seqs, edges, -2, edges[0])
        assert seqs == [list(range(13))]
        assert edges == [(2,10,'-')]
        seqs = [list(range(7)), list(range(100,107))]
        edges = [(4,102,'-'), (100,105,'-'), (103,3,'-'), (106,'R','+')]
        seqs, edges = self._loader.connectseqs(
            seqs, edges, 4, edges[-2]
        )
        assert seqs == [list(range(17))]
        assert edges == [(4,12,'-'), (10,15,'-'), (13,3,'-'), (16,'R','+')]

    def test_load(self):
        pat = '0-10_H7:3-HR:6-_2'
        seq = self._loader.load(pat)
        assert seq == ['U', 'A', '-2-', 'U']
        pat = '0-4:100-104_H2:102-H100:4-_LLLL_-2'
        seq = self._loader.load(pat, 2)
        assert seq == ['U', '+3-', 'A', '-2-', 'A']
        pat = '0-4:100-104_H2:102-H100:4+_LLLL_L'
        seq = self._loader.load(pat, 2)
        assert seq == ['U', 'L-', 'A', 'L+', 'A']
        pat = '0-6:100-106_H3:103+H6:R-H104:2+HR:5-_LLLL_-5'
        seq = self._loader.load(pat, 3)
        assert seq == ['U', '+6+', 'L-', 'U', 'L+', 'U']

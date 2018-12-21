import re
import math

class clustersequenceloader(object):
    """sequenceloader class"""

    def __init__(self):
        pass

    def load(self, filepath, length):
        """Load cluster sequence from the file specified."""
        src_col = 12
        dst_col = 13
        flag_col = 21
        cluster_col = 23

        regex = re.compile("__[US][US]$")
        
        seq = ["UNB"]*length
        with open(filepath) as fh:
            for line in fh:
                fields = line.strip().split()
                if not regex.search(fields[flag_col]):
                    continue
                seq[int(fields[src_col])-1] = fields[cluster_col]

        return seq

class localpatternsequenceloader(object):
    '''Sequence loader for local pattern'''

    def __init__(self):
        pass

    def patlen2delta(self, length):
        '''Convert pattern length to delta.'''
        '''This is a hack until it is fixed in Protein.findpattern2'''
        if length > 0:
            return -length
        else:
            return -length+1

    def computedelta(self, donor, acceptor):
        '''Compute delta, length of the given bond.'''
        if donor < acceptor:
            return int((acceptor-donor+1)/3)
        else:
            return int((acceptor-donor-2)/3)

    def findandreplace(self, stack, old, new):
        '''Replace old in stack with new. Stack is a (possibly 
        multi-level) list.'''
        try:
            idx = stack.index(old)
        except ValueError:
            for substack in stack:
                if substack == stack: #needed if stack is, e.g. str
                    continue
                self.findandreplace(substack, old, new)
        except (AttributeError, TypeError) as e:
            pass
        else:
            stack[idx] = new

    def findwindowsize(self, seqs, edges, length):
        '''Find window size from backbone segments and bonds.'''
        if len(seqs) == 2:
            return (len(seqs[0])-1)/2
        seq = seqs[0]
        #Find all potential central bonds
        candidates = []
        for i in range(math.ceil((len(seq)-1)/4), (len(seq)-1)//2):
            if (
                    (seq[i], seq[-i-1], '-') in edges or
                    (seq[i], seq[-i-1], '+') in edges
            ):
                candidates.append((seq[i], seq[-i-1]))
            elif (
                    (seq[-i-1], seq[i], '-') in edges or
                    (seq[-i-1], seq[i], '+') in edges
            ):
                candidates.append((seq[-i-1], seq[i]))
        if len(candidates) == 1:
            #There is only one candidate for the central bond
            return min(candidates[0][0], candidates[0][1])
        if type(length) is int:
            #Check which candidate matches the length
            for cand in candidates:
                d, a = cand
                if self.computedelta(d,a)==self.patlen2delta(length):
                    return min(d,a)
        #Can't do more. Raise error.
        raise ValueError('Cannot determine window size '
                         'from the inputs given.')
        
    def connectseqs(self, seqs, edges, length, centre):
        '''Connect two backbone sequences and translate edges.'''
        #Get diff=acceptor-donor from length
        delta = self.patlen2delta(length)
        if delta > 0:
            diff = 3*delta-1
        else:
            diff = 3*delta+2
        if centre[0] < centre[1]: #donor < acceptor
            x = centre[0] + diff - centre[1]
        else: #acceptor < donor
            x = centre[1] - diff - centre[0]
        #Now translate everything in seqs[1] by x
        newseqs = [list(range(seqs[1][-1] + x + 1))]
        newedges = []
        for edge in edges:
            newedge = tuple(e + x
                            if type(e) is int and e not in seqs[0]
                            else e for e in edge)
            newedges.append(newedge)
        return newseqs, newedges
        
    
    def load(self, pattern, window=None):
        '''Load local pattern sequence from pattern string.'''
        length = pattern.split('_')[-1]
        try:
            length = int(length)
        except ValueError:
            pass

        segs = re.findall(r'\d+-\d+', pattern)
        seqs = []
        for seg in segs:
            start, end = map(int, seg.split('-'))
            seqs.append(list(range(start, end+1)))

        bonds = re.findall(r'H\w+:\w+[+-]', pattern)
        edges = []
        for bond in bonds:
            don, acc = bond.strip('H').split(':')
            twist = acc[-1]
            acc = acc[:-1]
            if don != 'R':
                don = int(don)
            if acc != 'R':
                acc = int(acc)
            edges.append((don, acc, twist))

        #If we have two segments and the length is not L, connect
        #the two segments so that we can calculate the lengths of
        #other bonds in the window.
        if window is None:
            window = self.findwindowsize(seqs, edges, length)
        if len(seqs)==2 and type(length) is int:
            for edge in edges:
                if edge[0] == window or edge[1] == window:
                    centre = edge
            seqs, edges = self.connectseqs(
                seqs, edges, length, centre)

        #Populate seqs with bond length identifiers
        for edge in edges:
            if type(edge[0]) is int and type(edge[1]) is int:
                delta = self.computedelta(edge[0], edge[1])
                if abs(delta) > 6:
                    delta = 'L'
                else:
                    delta = format(delta, '+')
            else:
                delta = 'L'
            letter = delta + edge[2]
            self.findandreplace(seqs, edge[0], letter)
            self.findandreplace(seqs, edge[1], 'A')

        #If required, split the pattern again
        if len(seqs)==1 and len(seqs[0])-2 > 4*window:
            seq = seqs[0]
            seqs = [seq[0:2*window+1], seq[-2*window-1:len(seq)]]

        #Produce sequence
        words = []
        for seq in seqs:
            try:
                start = seq.index('A')
            except ValueError:
                for i, x in enumerate(seq):
                    if type(x) is str:
                        start = i - 1
                        break
            #start is now the index of an O atom
            #We now split the seq into chunks, such that each chunk
            #consists of atoms O, N, X
            m = 0 if start%3==0 else start%3-3
            chunks = [seq[max(i,0):i+3]
                      for i in range(m, len(seq), 3)]
            word = []
            for chunk in chunks:
                alist = [s=='A' for s in chunk]
                slist = [type(s) is str for s in chunk]
                nlist = [a^s for a,s in zip(alist, slist)]
                if any(nlist):
                    letter = chunk[nlist.index(True)]
                elif any(alist):
                    letter = 'A'
                else:
                    letter = 'U'
                word.append(letter)
            words.extend(word)
            
        return words

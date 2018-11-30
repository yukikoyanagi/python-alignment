#!/usr/bin/env python
#
# File: ClusterAlign.py
#
# Description: Class describing local pattern. Produced by running
#  locpat.Protein.findpattern()
#
# Author: Yuki Koyanagi
#

import re

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
                seq[int(fields[src_col])] = fields[cluster_col]

        return seq


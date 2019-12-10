__author__ = 'mahajrod'
import os
import sys
from collections import OrderedDict

if sys.version_info[0] == 3:
    izip = zip
else:
    from itertools import izip

import pandas as pd
from Bio import SearchIO, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from RouToolPa.Routines.Sequence import SequenceRoutines
#from Data.Nucleotides import back_degenerate_nucleotides
from RouToolPa.Collections.General import SynDict, IdList

import numpy as np


class WGARoutines(SequenceRoutines):
    def __init__(self):
        self.psl_query_id_column = 9           # zero-based
        self.psl_target_id_column = 13         # zero-based
        pass

    def label_maf(self, input, output, label_list):
        with open(input, "r") as in_fd, open(output, "r") as out_fd:
            for line in in_fd:
                if line[0] == "a":
                    out_fd.write(line)
                    pass
                else:
                    out_fd.write(line)
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
        SequenceRoutines.__init__(self)
        self.psl_query_id_column = 9           # zero-based
        self.psl_target_id_column = 13         # zero-based
        pass

    def label_maf(self, input, output, label_list, separator="."):
        with self.metaopen(input, "r") as in_fd, open(output, "w") as out_fd:
            label_id = 0

            for line in in_fd:
                if line[0] == "s":
                    out_fd.write("s %s%s%s" % (label_list[label_id], separator, line[2:]))
                    label_id += 1
                elif line[0] == "a":
                    out_fd.write(line)
                    label_id = 0
                elif (line[0] == "q") or (line[0] == "i"):
                    out_fd.write("%s %s%s%s" % (line[0], label_list[label_id-1], separator, line[2:]))
                else:
                    out_fd.write(line)

    def replace_ambigious_nucleotides_in_maf(self, input_maf, output_maf):

        with self.metaopen(input_maf, "r") as in_fd, self.metaopen(output_maf, "w") as out_fd:
            for line in in_fd:
                if line[0] == "s":
                    tmp = line.split(" ")
                    tmp[-1] = self.replace_ambiguous_nucleotides(tmp[-1])
                    out_fd.write(" ".join(tmp))
                else:
                    out_fd.write(line)


__author__ = 'mahajrod'
import os
import re
import sys
import math
import pickle

if sys.version_info[0] == 2:
    from string import maketrans

from random import randint
from copy import deepcopy
from collections import OrderedDict, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from RouToolPa.Collections.General import TwoLvlDict, SynDict, IdList, IdSet
from RouToolPa.Routines.Sequence import SequenceRoutines
from RouToolPa.Tools.Clustering import CDHit


class HaplotypeRoutines(SequenceRoutines):

    def __init__(self):
        SequenceRoutines.__init__(self)

    @staticmethod
    def find_indistinguishable_haplotypes(seq_file, haplotype_syn_file, output_prefix,
                                          threads=1, cdhit_dir=None):
        from RouToolPa.Routines import SequenceClusterRoutines

        sequence_fam = "%s.fam" % output_prefix
        haplotypes_fam = "%s.haplotypes.fam" % output_prefix
        indistinguishable_haplotypes_fam = "%s.haplotypes.indistinguishable.fam" % output_prefix

        CDHit.threads = threads
        CDHit.path = cdhit_dir
        CDHit.find_duplicates(seq_file, output_prefix)
        SequenceClusterRoutines.rename_elements_in_clusters(sequence_fam,
                                                            haplotype_syn_file,
                                                            haplotypes_fam,
                                                            keep_only_unique_elements=True)
        SequenceClusterRoutines.extract_clusters_by_size_from_file(haplotypes_fam,
                                                                   min_cluster_size=2,
                                                                   max_cluster_size=None,
                                                                   white_list_ids=None,
                                                                   out_file=indistinguishable_haplotypes_fam)






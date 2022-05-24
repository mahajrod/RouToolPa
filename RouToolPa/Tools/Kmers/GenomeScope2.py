#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import numpy as np
from scipy.signal import argrelextrema
import matplotlib


matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines import MatplotlibRoutines, MathRoutines


class GenomeScope2(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "genomescope2", path=path, max_threads=max_threads)

    def get_genome_size(self, input_histo, kmer_length, out_dir, output_prefix, ploidy=2, initial_haploid_coverage=None,
                        draw_fitted_hist=True, testing=True, max_kmer_coverage=100000000, cmd=None):
        # IMPORTANT! Not all options were implemented

        options = " -i %s " % input_histo
        options += " -p %i " % ploidy
        options += " -l %f " % initial_haploid_coverage if initial_haploid_coverage else ""
        options += " -k %i " % kmer_length
        options += " -n %s " % output_prefix
        options += " -m %i " % max_kmer_coverage
        options += " --fitted_hist " if draw_fitted_hist else ""
        options += " --testing " if testing else ""
        options += " -o %s " % out_dir

        self.execute(options, cmd=cmd if cmd else self.cmd)

#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import numpy as np
from scipy.signal import argrelextrema
import matplotlib


matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines import MatplotlibRoutines, MathRoutines


class Jellyfish(Tool):
    """
    Several subcommands are not implemented: qhisto, qdump, qmerge, cite
    Not all options were implemented for count subcommand

    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "jellyfish", path=path, max_threads=max_threads)

    def count(self, in_file, out_file, kmer_length=23, hash_size=1000000, count_both_strands=False,
              lower_count=None, upper_count=None):
        # IMPORTANT! Not all options were implemented
        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = "-m %i" % kmer_length
        options += " -s %s" % str(hash_size)
        options += " -t %i" % self.threads
        options += " -o %s" % out_file
        options += " -C" if count_both_strands else ""
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""

        input_files = [in_file] if isinstance(in_file, str) else in_file
        filetypes_list = []
        for filename in input_files:
            filetype = self.detect_filetype_by_extension(filename)
            if filetype not in filetypes_list:
                filetypes_list.append(filetype)

        if len(filetypes_list) > 1:
            raise ValueError("Mix of filetypes in input files: %s" % ",".join(filetypes_list))

        if filetypes_list[0] == "fasta" or filetypes_list[0] == "fastq" or filetypes_list[0] is None:
            if filetypes_list[0] is None:
                print("Warning!!! Type of input files was not recognized. Treating them as fasta...")
            cmd = "jellyfish count"
            options += " %s" % " ".join(input_files)
        elif filetypes_list[0] == "bzip":
            file_options = "bzcat %s" % " ".join(input_files)
            cmd = "%s | jellyfish count" % file_options
            options += " /dev/fd/0"
        elif filetypes_list[0] == "gz":
            file_options = "zcat %s" % " ".join(input_files)
            cmd = "%s | jellyfish count" % file_options
            options += " /dev/fd/0"

        self.execute(options, cmd=cmd)

    def stats(self, in_file, out_file, lower_count=None, upper_count=None):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -o %s" % out_file
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""
        options += " %s" % in_file

        self.execute(options, cmd="jellyfish stats")

    def histo(self, in_file, out_file, bin_width=1, lower_count=1, upper_count=100000000,
              include_absent_kmers=False):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -o %s" % out_file
        options += " -t %i" % self.threads
        options += " -l %i" % lower_count
        options += " -h %i" % upper_count
        options += " -f" if include_absent_kmers else ""
        options += " -i %i" % bin_width
        options += " %s" % in_file

        self.execute(options, cmd="jellyfish histo")

    def dump(self, in_file, out_file, lower_count=None, upper_count=None,
             column_format=True, tab_separator=True):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -o %s" % out_file
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""
        options += " -c" if column_format else ""
        options += " -t" if tab_separator else ""
        options += " %s" % in_file

        self.execute(options, cmd="jellyfish dump")

    def dump_kmers(self, in_file, output_prefix, lower_count=None, upper_count=None):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -t"     # use tab separator
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""
        options += " -c"    # column_format
        options += " %s" % in_file
        options += " | tee %s.counts | cut -f 1 " % output_prefix
        options += " > %s.kmers" % output_prefix

        self.execute(options, cmd="jellyfish dump")

    def query(self, sequence_file, jf_database, output_file):

        options = " -s %s" % sequence_file
        options += " -o %s" % output_file
        options += " %s" % jf_database

        self.execute(options, cmd="jellyfish query")

    def parallel_query(self, sequence_file, jf_database, output_file, splited_input_dir="splited_fasta",
                       per_sequence_analysis=False, parsing_mode="parse", splited_output_dir="splited_output",
                       retain_intermediate_file=False):

        self.split_fasta(sequence_file, splited_input_dir, num_of_recs_per_file=1 if per_sequence_analysis else None,
                         num_of_files=None if per_sequence_analysis else 5 * self.threads,
                         output_prefix="splited_fasta", parsing_mode=parsing_mode)
        options_list = []

        for filename in os.listdir(splited_input_dir):
            input_file = "%s/%s" % (splited_input_dir, filename)
            output_file = "%s/%s.kmer.counts" % (splited_output_dir, filename)

            options = " -s %s" % input_file
            options += " -o %s" % output_file
            options += " %s" % jf_database

            options_list.append(options)

        self.parallel_execute(options_list, cmd="jellyfish query")

        os.system("cat %s/* > %s" %(splited_output_dir, output_file))

        if not retain_intermediate_file:
            shutil.rmtree(splited_input_dir)
            shutil.rmtree(splited_output_dir)

    def get_kmer_list(self, in_file, out_prefix, kmer_length=23, hash_size=1000000, count_both_strands=False,
                      lower_count=None, upper_count=None):
        base_file = "%s.%i.jf" % (out_prefix, kmer_length)
        self.count(in_file, base_file, kmer_length=kmer_length, hash_size=hash_size,
                   count_both_strands=count_both_strands,
                   lower_count=lower_count, upper_count=upper_count)
        self.dump_kmers(in_file, out_prefix, lower_count=lower_count, upper_count=upper_count)

#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from aindex import *

from RouToolPa.Tools.Abstract import Tool

from RouToolPa.Tools.EMPHF import EMPHF
from RouToolPa.Tools.Kmers import Jellyfish


class AIndex(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "aindex", path=path, max_threads=max_threads)

    def compute_index(self, counts_file, pf_file, index_prefix,
                      full_hash=True):
        """

        :param counts_file:
        :param pf_file:
        :param index_prefix:
        :param full_hash:
        :return:

        creates following files: <index_prefix>.kmers.bin
                                 <index_prefix>.tf.bin
        """
        options = " %s" % counts_file
        options += " %s" % pf_file
        options += " %s" % index_prefix
        options += " %i" % self.threads
        options += " 0" if full_hash else ""

        self.execute(options=options, cmd="compute_index.exe")

    def compute_reads(self, forward_file, reverse_file, output, filetype="fastq"):
        """

        :param forward_file:
        :param reverse_file:
        :param output:
        :param filetype:
        :return:

        creates <output_file>
        """
        options = " %s" % forward_file
        options += " %s" % reverse_file
        options += " %s" % filetype
        options += " %s" % output

        self.execute(options=options, cmd="compute_reads.exe")

    def compute_aindex(self, reads_file, pf_file, freq_index_prefix,
                       read_index_prefix, kmer_length, index_scheme_file):

        options = " %s" % reads_file
        options += " %s " % pf_file
        options += " %s " % freq_index_prefix
        options += " %s " % read_index_prefix
        options += " %i " % self.threads
        options += " %i " % kmer_length
        options += " %s " % index_scheme_file

        self.execute(options=options, cmd="compute_aindex.exe")

    def create_index_from_jf(self, jf_db, kmer_length, output_prefix=None,
                             lower_count=1, upper_count=None,
                             forward_file=None, reverse_file=None,
                             filetype="fastq",
                             create_aindex=False,):

        out_pref = output_prefix if output_prefix else jf_db

        counts_file = "%s.counts" % out_pref
        kmer_file = "%s.kmers" % out_pref
        pf_file = "%s.pf" % out_pref
        reads_file = "%s.reads" % out_pref
        index_scheme_file = "%s.tf.bin" % out_pref

        print("Extracting kmers from jf database...")
        Jellyfish.dump_kmers(jf_db, output_prefix, lower_count=lower_count,
                             upper_count=upper_count)
        print("Creating index...")
        EMPHF.compute_mphf_seq(kmer_file, pf_file)
        self.compute_index(counts_file, pf_file, output_prefix)

        if create_aindex:
            print("Creating AIndex...")
            self.compute_reads(forward_file, reverse_file, reads_file, filetype=filetype)
            self.compute_aindex(reads_file, pf_file, output_prefix,
                                output_prefix, kmer_length, index_scheme_file)

    def create_index_from_fastq(self, forward_file, reverse_file,
                                kmer_length, output_prefix,
                                lower_count=2, upper_count=None,
                                filetype="fastq",
                                create_aindex=False,
                                hash_size="10G"):

        out_pref = "%s.%i" % (output_prefix, kmer_length)
        jf_db = "%s.jf" % out_pref

        Jellyfish.threads = self.threads
        print("Creating jf database...")
        Jellyfish.count([forward_file, reverse_file], jf_db,
                        kmer_length=kmer_length, hash_size=hash_size,
                        count_both_strands=True,
                        lower_count=lower_count, upper_count=upper_count)

        self.create_index_from_jf(jf_db, kmer_length,
                                  output_prefix=out_pref,
                                  lower_count=lower_count,
                                  upper_count=upper_count,
                                  forward_file=forward_file,
                                  reverse_file=reverse_file,
                                  filetype=filetype,
                                  create_aindex=create_aindex)


if __name__ == "__main__":
    pass

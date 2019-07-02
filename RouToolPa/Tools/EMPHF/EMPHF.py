#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from RouToolPa.Tools.Abstract import Tool


class EMPHF(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "mphf", path=path, max_threads=max_threads)

    def compute_mphf_seq(self, kmer_file, pf_file):
        """

        :param kmer_file: input file with kmer list
        :param pf_file: pf file to create
        :return:

        creates pf file

        """
        options = " %s" % kmer_file
        options += " %s" % pf_file

        self.execute(options=options, cmd="compute_mphf_seq")

#!/usr/bin/env python

import os
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines import DrawingRoutines


class Sambamba(Tool):
    """
    Class for sambamba
    Several subcommands are not implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "sambamba", path=path, max_threads=max_threads)

        # bam/sam flag values:
        self.bam_flags = {
                          "read_paired": 1,
                          "read_mapped_in_proper_pair": 2,
                          "read_unmapped": 4,
                          "mate_unmapped": 8,
                          "read_reverse_strand": 16,
                          "mate_reverse_strand": 32,
                          "first_in_pair": 64,
                          "second_in_pair": 128,
                          "not_primary_alignment": 256,
                          "read_fails_platform/vendor_quality_checks": 512,
                          "read_is_PCR_or_optical_duplicate": 1024,
                          "supplementary_alignment": 2048
                          }

    def mkdup(self, input_bam, output_bam, verbose=None, tmp_dir=None, hash_table_size=None,
              overflow_list_size=None, io_buffer=3000):
        """

        :param input_bam:
        :param output_bam:
        :param verbose:
        :param tmp_dir:
        :param hash_table_size:
        :param overflow_list_size:
        :param io_buffer: size of buffers for input and output in megabytes
        :return:
        """

        options = " -t %i" % self.threads
        options += " --show-progress" if verbose else ""
        options += " --tmpdir %s" % tmp_dir if tmp_dir else ""
        options += " --hash-table-size %i" % hash_table_size if hash_table_size else ""
        options += " --overflow-list-size %i" % overflow_list_size if overflow_list_size else ""
        options += "--io-buffer-size %i" % io_buffer if io_buffer else ""

        options += " %s" % input_bam
        options += " %s" % output_bam

        self.execute(options=options, cmd="sambamba markdup")

    def index(self, input_bam, output_bam):
        options = " -t %i" % self.threads

        options += " %s" % input_bam
        options += " %s" % output_bam

        self.execute(options=options, cmd="sambamba index")

#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from RouToolPa.Tools.Abstract import Tool


class Supernova(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "supernova", path=path, max_threads=max_threads)

    def generate_fasta(self, assembly_dir, output_prefix, minsize=1000, headers="full"):

        options = " mkoutput"
        options += " --asmdir %s" % assembly_dir
        options += " --minsize %i" % minsize
        options += " --headers %s" % headers

        options_list = []

        for out_type in "raw", "megabubbles", "pseudohap", "pseudohap2":
            tmp_options = options
            tmp_options += " --style %s" % out_type
            tmp_options += " --outprefix %s.%s" % (output_prefix, out_type)

            options_list.append(tmp_options)

        self.parallel_execute(options_list=options_list)
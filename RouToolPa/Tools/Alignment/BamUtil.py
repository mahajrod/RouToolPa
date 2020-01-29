#!/usr/bin/env python
import os
import shutil
from RouToolPa.Tools.Abstract import Tool


class BamUtil(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bam clipOverlap", path=path, max_threads=max_threads)

    def parse_options(self, input, output, poolsize=None):

        options = " --in %s" % input
        options += " --out %s" % output
        options += " --poolSize %i" % poolsize if poolsize else ""

        return options

    def clipoverlap(self, input, output, poolsize=None):
        from RouToolPa.Tools.Samtools import SamtoolsV1
        options = self.parse_options(input, output, poolsize=poolsize)

        self.execute(options=options, cmd="bam clipOverlap")
        SamtoolsV1.index(output)

    def parallel_clipoverlap(self, input_dir, output_dir, samples_list, bam_suffix="", poolsize=None):
        from RouToolPa.Tools.Samtools import SamtoolsV1
        samples_to_handle = samples_list if samples_list else self.get_sample_list(input_dir)

        self.safe_mkdir(output_dir)

        options_list = []
        samtools_option_list = []

        for sample in samples_to_handle:
            sample_dir = "%s/%s/" % (output_dir, sample)
            self.safe_mkdir(sample_dir)
            input_bam = "%s/%s/%s%s.bam" % (sample_dir, sample, sample, bam_suffix)
            output_bam = "%s/%s.bam" % (output_dir, sample)

            options_list.append(self.parse_options(input_bam, output_bam, poolsize=poolsize))
            samtools_option_list.append(output_bam)

        self.parallel_execute(options_list=options_list)
        self.parallel_execute(options_list=samtools_option_list, cmd="samtools index")


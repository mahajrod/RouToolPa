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

        options = self.parse_options(input, output, poolsize=poolsize)

        self.execute(options=options, cmd="bam clipOverlap")

#!/usr/bin/env python
import os
from RouToolPa.Tools.Abstract import Tool


class CombineVariants4(Tool):

    def __init__(self, max_threads=4, max_memory=None, timelog=None):
        Tool.__init__(self,
                      "gatk --java-options %s CombineVariants" % max_memory if max_memory else "gatk CombineVariants",
                      max_threads=max_threads, max_memory=max_memory,
                      timelog=timelog)
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html

    def combine_from_same_source(self, reference_file, input_filelist, output_vcf):

        options = " -R %s" % reference_file
        options += " -nt %i" % self.threads
        options += " --variant %s" % (input_filelist if isinstance(input_filelist, str) else " --variant ".join(input_filelist))
        options += " --genotypemergeoption UNSORTED"
        options += " -o %s" % output_vcf

        self.execute(options=options)

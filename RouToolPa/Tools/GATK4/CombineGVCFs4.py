#!/usr/bin/env python

__author__ = 'mahajrod'

from RouToolPa.Tools.Abstract import Tool


class CombineGVCFs4(Tool):
    def __init__(self, max_threads=4, max_memory=None, timelog=None):
        Tool.__init__(self,
                      "gatk CombineGVCFs",
                      max_threads=max_threads, max_memory=max_memory,
                      timelog=timelog)

    def parse_options(self, reference, gvcf_list, output, extension_list=["g.vcf",]):

        options = " -O %s" % output
        options += " -R %s" % reference

        for gvcf in self.make_list_of_path_to_files_by_extension(gvcf_list,
                                                                 extension_list=extension_list,
                                                                 recursive=False, return_absolute_paths=True):
            options += " --variant %s" % gvcf

        return options

    def combine(self, reference, gvcf_list, output, extension_list=["g.vcf",]):

        options = self.parse_options(reference,
                                     gvcf_list, output,
                                     extension_list=extension_list)

        self.execute(options=options)

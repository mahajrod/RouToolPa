#!/usr/bin/env python

__author__ = 'mahajrod'

from RouToolPa.Tools.Abstract import Tool


class GenomicsDBImport4(Tool):
    def __init__(self, max_threads=4, max_memory=None, timelog=None):
        Tool.__init__(self,
                      "gatk GenomicsDBImport",
                      max_threads=max_threads, max_memory=max_memory,
                      timelog=timelog)

    def parse_options(self, gvcf_list, output_dbi_dir, interval_list, extension_list=["g.vcf",]):

        options = " --genomicsdb-workspace-path %s" % output_dbi_dir
        options += " --intervals %s" % ",".join(interval_list)

        for gvcf in self.make_list_of_path_to_files_by_extension(gvcf_list,
                                                                 extension_list=extension_list,
                                                                 recursive=False, return_absolute_paths=True):
            options += " --variant %s" % gvcf

        return options

    def create_db(self, gvcf_list, output_dbi_dir, interval_list, extension_list=["g.vcf",]):

        options = self.parse_options(gvcf_list, output_dbi_dir,
                                     interval_list=interval_list,
                                     extension_list=extension_list)

        self.execute(options=options)

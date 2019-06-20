#!/usr/bin/env python

__author__ = 'mahajrod'

from RouToolPa.Tools.Abstract import Tool


class BaseRecalibrator4(Tool):

    def __init__(self, max_threads=4, max_memory=None, timelog=None):
        Tool.__init__(self,
                      "gatk --java-options %s BaseRecalibrator" % max_memory if max_memory else "gatk BaseRecalibrator",
                      max_threads=max_threads, max_memory=max_memory,
                      timelog=timelog)

    def get_recalibration_table(self, reference, alignment, output_table, known_sites_vcf, BQSR=None,
                                include_region_id_file=None, exclude_region_id_file=None):

        # TODO: add rest of  options
        options = ""
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " -nct %i" % self.threads
        options += " -BQSR %s" % BQSR if BQSR else ""
        options += " -knownSites %s" % known_sites_vcf if isinstance(known_sites_vcf, str) else " -knownSites ".join(known_sites_vcf)
        options += " -L %s" % include_region_id_file if include_region_id_file else ""
        options += " -XL %s" % exclude_region_id_file if exclude_region_id_file else ""
        options += " -o %s" % output_table

        self.execute(options)

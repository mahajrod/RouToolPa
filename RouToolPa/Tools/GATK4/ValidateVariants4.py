#!/usr/bin/env python
import os
from RouToolPa.Tools.Abstract import Tool


class ValidateVariants4(Tool):
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_filters_VariantFiltration.html
    # default filters for indel and snp filtration were taken from GATK BestPractice
    def __init__(self, max_threads=4, max_memory=None, timelog=None):
        Tool.__init__(self, "GenomeAnalysisTK.jar -T ValidateVariants",
                      max_threads=max_threads, max_memory=max_memory,
                      timelog=timelog)

    @staticmethod
    def parse_common_options(reference, input_vcf, exclude_list=[], dbsnp=None, input_is_gvcf=False,
                             create_index=False):
        options = " -R %s" % reference
        options += " -V %s" % input_vcf
        options += " --dbsnp %s" % dbsnp if dbsnp else ""
        options += " --create-output-variant-index" if create_index else False
        for entry in exclude_list:
            options += " --validation-type-to-exclude %s" % entry

        options += " --validate-GVCF" if input_is_gvcf else ""

        return options

    def test_vcf_format(self, reference, input_vcf, input_is_gvcf=False, create_index=False):

        options = self.parse_common_options(reference, input_vcf, exclude_list=["ALL", ], input_is_gvcf=input_is_gvcf,
                                            create_index=create_index)

        self.execute(options=options)

    def index_vcf(self, reference, input_vcf, input_is_gvcf=False):
        self.test_vcf_format(reference, input_vcf, input_is_gvcf=input_is_gvcf)


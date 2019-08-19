import os
from RouToolPa.Tools.Abstract import Tool


class SelectVariants4(Tool):

    def __init__(self,  path="", max_threads=4,  max_memory=None, timelog=None):
        Tool.__init__(self, "gatk SelectVariants", path=path,
                      max_threads=max_threads,max_memory=max_memory,
                      timelog=timelog)
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_SelectVariants.html
    # selectType
    # INDEL
    # SNP
    # MIXED
    # MNP
    # SYMBOLIC
    # NO_VARIATION

    def select_variants(self, reference_file, input_vcf, output_vcf, vartype=None, varfilter=None):

        options = " -R %s" % reference_file
        options += " -V %s" % input_vcf
        options += " --select-type-to-include \'%s\'" % vartype if vartype else ""
        options += " -select \'%s\'" % varfilter if varfilter else ""
        options += " -O %s" % output_vcf

        self.execute(options=options)

    def get_SNP(self, reference_file, input_vcf, output_vcf):
        self.select_variants(reference_file, input_vcf, output_vcf, vartype="SNP")

    def get_indel(self, reference_file, input_vcf, output_vcf):
        self.select_variants(reference_file, input_vcf, output_vcf, vartype="INDEL")

    def remove_entries_with_filters(self, reference_file, input_vcf, output_vcf):
        options = " -R %s" % reference_file
        options += " -V %s" % input_vcf
        options += " --exclude-filtered"
        options += " -O %s" % output_vcf
        self.execute(options=options)

    def get_STR(self, reference_file, input_vcf, output_vcf):
        #Use STR filtration based on GATK prediction very carefully as it counts even AAA -> AA as STR
        self.select_variants(reference_file, input_vcf, output_vcf, varfilter="STR")

    def get_nonSTR(self, reference_file, input_vcf, output_vcf):
        self.select_variants(reference_file, input_vcf, output_vcf, varfilter="vc.hasAttribute(\"STR\") == 0")


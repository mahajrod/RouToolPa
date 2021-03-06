#!/usr/bin/env python
__author__ = 'mahajrod'

from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Tools.Abstract import JavaTool


class RealignerTargetCreator(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T RealignerTargetCreator", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def create(self, reference, alignment, output="forIndelRealigner.intervals", known_indels_vcf=None,
               max_interval_size=None, min_reads_cov=None, mismatch_fraction=None, window_size=None,
               default_base_qualities=None):

        options = " -nt %i" % self.threads
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " -o %s" % output
        options += " --known %s" % known_indels_vcf if known_indels_vcf else ""
        options += " -maxInterval %i" % max_interval_size if max_interval_size else ""
        options += " -minReads %i" % min_reads_cov if min_reads_cov else ""
        options += " -mismatch %i" % mismatch_fraction if mismatch_fraction else ""
        options += " -window %i" % window_size if window_size else ""
        options += " --defaultBaseQualities %i" % default_base_qualities if default_base_qualities else ""

        self.execute(options)

if __name__ == "__main__":
    import os
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc"
    os.chdir(workdir)
    #gatk_dir =
    reference = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.fasta"
    alignment = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc/001_trimmed_sorted_rm_pcr_chrom.bam"
    RealignerTargetCreator = RealignerTargetCreator(jar_path="/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0")
    RealignerTargetCreator.create(reference, alignment)
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
            input_bam = "%s/%s%s.bam" % (sample_dir, sample, bam_suffix)
            output_bam = "%s/%s/%s.clipped.bam" % (output_dir, sample, sample)

            options_list.append(self.parse_options(input_bam, output_bam, poolsize=poolsize))
            samtools_option_list.append(output_bam)

        self.parallel_execute(options_list=options_list)
        self.parallel_execute(options_list=samtools_option_list, cmd="samtools index")

    def get_stats_from_bam(self):
        """
         --in  --pBaseQC  --cBaseQC
        -in : the SAM/BAM file to calculate stats for
        Types of Statistics that can be generated:
                --basic         : Turn on basic statistic generation
                --qual          : Generate a count for each quality (displayed as non-phred quality)
                --phred         : Generate a count for each quality (displayed as phred quality)
                --pBaseQC       : Write per base statistics as Percentages to the specified file. (use - for stdout)
                                  pBaseQC & cBaseQC cannot both be specified.
                --cBaseQC       : Write per base statistics as Counts to the specified file. (use - for stdout)
                                  pBaseQC & cBaseQC cannot both be specified.
        Optional Parameters:
                --maxNumReads   : Maximum number of reads to process
                                  Defaults to -1 to indicate all reads.
                --unmapped      : Only process unmapped reads (requires a bamIndex file)
                --bamIndex      : The path/name of the bam index file
                                  (if required and not specified, uses the --in value + ".bai")
                --regionList    : File containing the regions to be processed chr<tab>start_pos<tab>end_pos.
                                  Positions are 0 based and the end_pos is not included in the region.
                                  Uses bamIndex.
                --excludeFlags  : Skip any records with any of the specified flags set
                                  (specify an integer representation of the flags)
                --requiredFlags : Only process records with all of the specified flags set
                                  (specify an integer representation of the flags)
                --noeof         : Do not expect an EOF block on a bam file.
                --params        : Print the parameter settings.
        Optional phred/qual Only Parameters:
                --withinRegion  : Only count qualities if they fall within regions specified.
                                  Only applicable if regionList is also specified.
        Optional BaseQC Only Parameters:
                --baseSum       : Print an overall summary of the baseQC for the file to stderr.
                --bufferSize    : Size of the pileup buffer for calculating the BaseQC parameters.
                                  Default: 1024
                --minMapQual    : The minimum mapping quality for filtering reads in the baseQC stats.
                --dbsnp         : The dbSnp file of positions to exclude from baseQC analysis.
        """
        pass


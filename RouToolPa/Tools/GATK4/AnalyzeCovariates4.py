#!/usr/bin/env python
import os
from RouToolPa.Tools.Abstract import Tool


class AnalyzeCovariates4(Tool):

    def __init__(self,  max_threads=4, max_memory=None, timelog=None):
        Tool.__init__(self,
                      "gatk --java-options %s AnalyzeCovariates" % max_memory if max_memory else "gatk AnalyzeCovariates",
                      max_threads=max_threads, max_memory=max_memory,
                      timelog=timelog)

    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html

    def parse_options(self, reference, BQSR_table=None, output_plots=None, before_table=None, after_table=None,
                      ignoreLMT=False, csv_out=None): #, known_sites=None):

        options = " -R %s" % reference
        options += " -BQSR %s" % BQSR_table if BQSR_table else ""
        options += " -plots %s" % output_plots if output_plots else ""
        options += " -before %s" % before_table if before_table else ""
        options += " -after %s" % after_table if after_table else ""
        options += " -ignoreLMT" if ignoreLMT else ""
        options += " -csv %s" % csv_out if csv_out else ""

        #if known_sites:
        #    if isinstance(known_sites, str):
        #        options += " -knownSites %s" % known_sites
        #    else:
        #        for entry in known_sites:
        #            options += " -knownSites %s" % entry

        return options

    def plot_single_recall_table(self, reference, BQSR_table, output_plots, csv_out=None):
        """
        gatk AnalyzeCovariates \
          -R myrefernce.fasta \
          -BQSR myrecal.table \
          -plots BQSR.pdf
        """

        options = self.parse_options(reference, BQSR_table=BQSR_table, output_plots=output_plots, csv_out=csv_out)

        self.execute(options=options)

    def plot_two_recall_table(self, reference, before_table, after_table, output_plots, csv_out=None):
        """
        Plot before (first pass) and after (second pass) recalibration tables to compare them
         gatk AnalyzeCovariates \
              -R myrefernce.fasta \
              -before recal2.table \
              -after recal3.table \
              -plots recalQC.pdf
        """

        options = self.parse_options(reference, before_table=before_table, csv_out=csv_out,
                                     after_table=after_table, output_plots=output_plots)

        self.execute(options=options)

    def plot_three_recall_table(self, reference, BQSR_table, before_table, after_table, output_plots,
                                ignoreLMT=False, csv_out=None):
        """
        Plot up to three recalibration tables for comparison

         # You can ignore the before/after semantics completely if you like (if you do, add -ignoreLMT
         # to avoid a possible warning), but all tables must have been generated using the same parameters.

         gatk AnalyzeCovariates \
              -R myrefernce.fasta \
              -ignoreLMT \
              -BQSR recal1.table \   # you can discard any two
              -before recal2.table \
              -after recal3.table \
              -plots myrecals.pdf
        """

        options = self.parse_options(reference, BQSR_table=BQSR_table, before_table=before_table, csv_out=csv_out,
                                     after_table=after_table, output_plots=output_plots, ignoreLMT=ignoreLMT)

        self.execute(options=options)


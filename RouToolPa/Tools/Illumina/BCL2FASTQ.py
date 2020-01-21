#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
from RouToolPa.Tools.Abstract import Tool


class BCL2FASTQ(Tool):
    """
    TODO: NOT ALL options were implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bcl2fastq", path=path, max_threads=max_threads)

    def parse_options(self, run_dir, output_dir, sample_sheet, intensities_dir=None, barcode_mismatches=None,
                      loading_threads=None, processing_threads=None, writing_threads=None):

        options = " -r %i" % loading_threads if loading_threads else self.threads
        options += " -p %i" % processing_threads if processing_threads else self.threads
        options += " -w %i" % writing_threads if writing_threads else self.threads

        options += "-i %s" % intensities_dir if intensities_dir else ""# Ex 191228_M02435_0069_000000000-CP7BV/Data/Intensities/BaseCalls/
        options += " --output-dir %s" % output_dir
        options += " -R %s" % run_dir # Ex 191228_M02435_0069_000000000-CP7BV
        options += " --sample-sheet %s" % sample_sheet

        options += " --barcode-mismatches %i" % barcode_mismatches if barcode_mismatches else "" #Allowed: 0,1(default),2

        return options

    def demultiplex(self, run_dir, output_dir, sample_sheet, intensities_dir=None, barcode_mismatches=None,
                    loading_threads=None, processing_threads=None, writing_threads=None):

        options = self.parse_options(run_dir, output_dir, sample_sheet,
                                     intensities_dir=intensities_dir,
                                     barcode_mismatches=barcode_mismatches,
                                     loading_threads=loading_threads,
                                     processing_threads=processing_threads,
                                     writing_threads=writing_threads)

        self.execute(options=options)


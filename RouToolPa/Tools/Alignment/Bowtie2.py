#!/usr/bin/env python
import os
from RouToolPa.Tools.Abstract import Tool


class Bowtie2(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bowtie2", path=path, max_threads=max_threads)

    def index(self, reference, index_name):
        options = "%s %s" % (reference, index_name)
        self.execute(options, cmd="bowtie2-build")

    def align(self,
              bowtie2_index,
              forward_reads_list=None,
              reverse_reads_list=None,
              unpaired_reads_list=None,
              quality_score="phred33",
              alignment_mode="very-sensitive",
              find_discordant_alignments=True,
              find_separated_alignments=True,
              max_insert_size=None,
              output_prefix="alignment",
              output_format="bam",
              read_group_name="reads",
              PU="x",
              SM="sample",
              platform="Illumina",
              LB="x",
              mark_duplicates=True,
              sort_by_coordinate=False,
              sort_by_name=False,
              max_per_sorting_thread_memory=None,
              softclipping_penalty=None, local_alignment=False):

        options = " -p %i" % self.threads
        options += " --local" if local_alignment else ""
        options += " --%s" % alignment_mode
        options += " --%s" % quality_score
        options += " -x %s" % bowtie2_index
        options += " -X %i" % max_insert_size if max_insert_size else ""
        options += " --rg-id %s" % read_group_name
        options += " --rg \'PU:%s	SM:%s	PL:%s	LB:%s\'" % (PU, SM, platform, LB)
        options += " --no-discordant" if not find_discordant_alignments else ""
        options += " --no-mixed" if not find_separated_alignments else ""
        options += " -1 %s -2 %s" % (forward_reads_list if isinstance(forward_reads_list, str) else ",".join(forward_reads_list),
                                     reverse_reads_list if isinstance(reverse_reads_list, str) else",".join(reverse_reads_list)) \
            if forward_reads_list and reverse_reads_list else ""
        options += " -U %s" % ",".join(unpaired_reads_list) if unpaired_reads_list else ""
        options += " 2>%s.stats" % output_prefix
        if sort_by_coordinate or sort_by_name:
            if sort_by_coordinate and sort_by_name:
                raise ValueError("Sorting by both coordinate and read name was requested")
            options += " | samtools view -b | samtools sort"
            if sort_by_name:
                options += " -n"
            options += " -@ %i" % self.threads
            options += " -m %s" % max_per_sorting_thread_memory if max_per_sorting_thread_memory else self.max_per_thread_memory
            options += " -O %s" % output_format.upper()

        if mark_duplicates:
            options += " | samtools fixmate"
            options += " -@ %i" % self.threads
            options += " -m - - "

            options += " | samtools sort"
            options += " -@ %i" % self.threads
            options += " -m %s" % max_per_sorting_thread_memory if max_per_sorting_thread_memory else self.max_per_thread_memory

            options += " | samtools markdup"
            options += " -@ %i" % self.threads
            options += " - %s.%s" % (output_prefix, output_format)
        else:
            options += " > %s.%s" % (output_prefix, output_format)

        self.execute(options)


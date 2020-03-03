#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from RouToolPa.Tools.Abstract import Tool


class RaGOO(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "ragoo.py", path=path, max_threads=max_threads)

    def parse_options(self, target, reference, dont_merge_unscafollded_contigs=True):
        # do not use full pathes as ragoo don't understand them
        options = " -t %i" % self.threads
        options += " -C" if dont_merge_unscafollded_contigs else dont_merge_unscafollded_contigs

        options += " %s" % target
        options += " %s" % reference

        return options

    def scaffold(self, target, reference, stat_prefix="stats", target_label="target",
                 reference_label="reference", ragoo_label="target_ragoo",
                 dont_merge_unscafollded_contigs=True):

        ragoo = "ragoo_output/ragoo.fasta"

        options = self.parse_options(target, reference, dont_merge_unscafollded_contigs=dont_merge_unscafollded_contigs)

        self.execute(options=options)

        self.get_stats_from_assemblies([reference, target, ragoo], [reference_label, target_label, ragoo_label],
                                       stat_prefix, thresholds_list=(0, 100, 250, 500, 1000))


if __name__ == "__main__":
    pass

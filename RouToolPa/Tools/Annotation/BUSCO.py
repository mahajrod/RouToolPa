#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from collections import OrderedDict
from RouToolPa.Collections.General import TwoLvlDict, IdSet

import matplotlib.pyplot as plt

from RouToolPa.Parsers.BUSCO import BUSCOtable
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines import MatplotlibRoutines



class BUSCO(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "run_BUSCO.py", path=path, max_threads=max_threads)
        self.kingdoms = ["euk", "bac", "mito", "arc"]
        self.status_list = ["Complete",
                            "Duplicated",
                            "Fragmented",
                            "Missing"]

    def parse_options(self, input_fasta, label, busco_db, augustus_species, mode=None):

        options = " -c %i" % self.threads
        options += " -l %s" % busco_db
        options += " -m %s" % mode if mode else ""
        options += " -o %s" % label if label else ""
        options += " -sp %s" % augustus_species if augustus_species else ""
        options += " -i %s" % input_fasta

        return options

    def assess(self, input_fasta, label, busco_db, augustus_species, mode=None):

        options = self.parse_options(input_fasta, label, busco_db, augustus_species, mode=mode)

        self.execute(options=options)

    def assess_multiple_genomes(self, input_fasta_list, busco_db, augustus_species, label_list=None,
                                output_dir="./", mode=None):

        outdir = os.path.abspath(output_dir)

        for fasta, label in zip(input_fasta_list, label_list if label_list else ["A%i" % i for i in range(1,
                                                                                                          len(input_fasta_list) + 1)]):
            genome_fasta = os.path.abspath(fasta)
            options = self.parse_options(genome_fasta, label, busco_db, augustus_species, mode=mode)
            busco_dir = "%s/%s/" % (outdir, label)
            self.safe_mkdir(busco_dir)

            self.execute(options=options, directory=busco_dir)

    def compare_multiple_genome_results(self, busco_file_list, output_prefix, label_list=None,
                                        black_scaffold_list=(), white_scaffold_list=()):

        busco_table_dict = OrderedDict()
        gene_id_dict = OrderedDict()
        counts_dict = OrderedDict()

        output_path_list = self.split_filename(output_prefix)

        pairwise_overlaps_dir = "%s/pairwise_overlaps/" % (output_path_list[0] if output_path_list[0] else ".")
        pairwise_overlap_counts_dir = "%s/pairwise_overlap_counts/" % (output_path_list[0] if output_path_list[0] else ".")
        self.safe_mkdir(pairwise_overlaps_dir)
        self.safe_mkdir(pairwise_overlap_counts_dir)

        lllabels_list = label_list if label_list else ["A%i" % i for i in range(1, len(busco_file_list) + 1)]

        for busco_table, label in zip(busco_file_list, lllabels_list):
            busco_table_dict[label] = BUSCOtable(in_file=busco_table, black_list=black_scaffold_list,
                                                 white_list=white_scaffold_list)

            gene_id_dict[label] = OrderedDict()
            counts_dict[label] = OrderedDict()

            gene_id_dict[label], counts_dict[label] = busco_table_dict[label].count_statuses()

        # TODO: draw piecharts


        # TODO: count overlaps

        pairwise_overlap_dict = OrderedDict()
        count_pairwise_overlap_dict = OrderedDict()
        for label1 in lllabels_list:
            for label2 in lllabels_list:
                if label1 == label2:
                    continue
                overlap_id = "%s_vs_%s" % (label1, label2)
                pairwise_overlap_dict[overlap_id] = TwoLvlDict()
                count_pairwise_overlap_dict[overlap_id] = TwoLvlDict()
                for status1 in self.status_list:
                    pairwise_overlap_dict[overlap_id]["%s@%s" % (label1, status1)] = OrderedDict()
                    count_pairwise_overlap_dict[overlap_id]["%s@%s" % (label1, status1)] = OrderedDict()
                    for status2 in self.status_list:
                        pairwise_overlap_dict[overlap_id]["%s@%s" % (label1, status1)]["%s@%s" % (label2, status2)] = IdSet(gene_id_dict[label1][status1] & gene_id_dict[label2][status2])
                        count_pairwise_overlap_dict[overlap_id]["%s@%s" % (label1, status1)]["%s@%s" % (label2, status2)] = len(pairwise_overlap_dict[overlap_id]["%s@%s" % (label1, status1)]["%s@%s" % (label2, status2)])
                        pairwise_overlap_dict[overlap_id]["%s@%s" % (label1, status1)]["%s@%s" % (label2, status2)].write("%s/%s.%s_vs_%s.ids" % (pairwise_overlaps_dir, output_prefix, "%s@%s" % (label1, status1), "%s@%s" % (label2, status2)))

                count_pairwise_overlap_dict[overlap_id].write("%s/%s.overlap.%s.tsv" % (pairwise_overlap_counts_dir, output_prefix, overlap_id))

        if 2 <= len(busco_file_list) <= 3:
            fig, subplot_list = plt.subplots(2, 2, figsize=(6, 6))
            plt.suptitle("Overlaps for BUSCO categories between assemblies/genomes")
            print subplot_list
            for status, index in zip(self.status_list, range(0, 4)):

                plt.sca(subplot_list[index // 2][index % 2])
                plt.title(status)
                MatplotlibRoutines.venn_diagram_from_sets(gene_id_dict[lllabels_list[0]][status],
                                                          gene_id_dict[lllabels_list[1]][status],
                                                          set3=gene_id_dict[lllabels_list[2]][status] if len(lllabels_list) > 2 else None,
                                                          set_labels=lllabels_list, set_colors=["red", "yellow", "green"],
                                                          output_prefix=None, extensions=("png",), title=None)

            plt.savefig("%s.venn.png" % output_prefix)

            plt.close()

if __name__ == "__main__":
    pass

__author__ = 'mahajrod'
import os
import re
import sys
import math
import pickle
import pandas as pd

if sys.version_info[0] == 2:
    from string import maketrans

from RouToolPa.Collections.General import TwoLvlDict, SynDict, IdList, IdSet
from RouToolPa.Routines.SequenceCluster import SequenceClusterRoutines


class HaplotypeRoutines(SequenceClusterRoutines):

    def __init__(self):
        SequenceClusterRoutines.__init__(self)

    @staticmethod
    def find_indistinguishable_haplotypes(seq_file,  output_prefix, haplotype_syn_file=None,
                                          threads=1, cdhit_dir=None):
        from RouToolPa.Routines import SequenceClusterRoutines
        from RouToolPa.Tools.Clustering import CDHit

        sequence_fam = "%s.fam" % output_prefix
        haplotypes_fam = "%s.haplotypes.fam" % output_prefix
        indistinguishable_haplotypes_fam = "%s.haplotypes.indistinguishable.fam" % output_prefix

        CDHit.threads = threads
        CDHit.path = cdhit_dir
        CDHit.find_duplicates(seq_file, output_prefix)
        if haplotype_syn_file:
            SequenceClusterRoutines.rename_elements_in_clusters(sequence_fam,
                                                                haplotype_syn_file,
                                                                haplotypes_fam,
                                                                keep_only_unique_elements=True)
        SequenceClusterRoutines.extract_clusters_by_size_from_file(haplotypes_fam if haplotype_syn_file else sequence_fam,
                                                                   min_cluster_size=2,
                                                                   max_cluster_size=None,
                                                                   white_list_ids=None,
                                                                   out_file=indistinguishable_haplotypes_fam)

    @staticmethod
    def prepare_template_for_popart(alignment_file, output_file,  haplotype_fam_file=None, traits_file=None,
                                    whitelist_file=None):
        from RouToolPa.Parsers.Sequence import CollectionSequence
        sequence_collection = CollectionSequence(in_file=alignment_file, parsing_mode="parse")
        sequence_collection.get_stats_and_features(count_gaps=False, sort=False)
        whitelist = IdSet(filename=whitelist_file)
        alignment_len = sequence_collection.seq_lengths["length"].unique()
        if len(alignment_len) > 1:
            raise ValueError("ERROR!!! Sequences in alignment have different lengths!")
        alignment_len = alignment_len[0]

        haplotype_selected_sequence_dict = SynDict()
        haplotypes_without_sequences_ids = IdList()

        traits_df = pd.read_csv(traits_file, sep="\t", index_col=0) if traits_file else pd.DataFrame()

        if haplotype_fam_file:
            haplotype_dict = SynDict(filename=haplotype_fam_file, split_values=True)
            for haplotype_id in haplotype_dict:
                for sequence_id in haplotype_dict[haplotype_id]:
                    if sequence_id in sequence_collection.records:
                        haplotype_selected_sequence_dict[haplotype_id] = sequence_id
                        break
                else:
                    haplotypes_without_sequences_ids.append(haplotype_id)
        else:
            haplotype_dict = dict([(entry, [entry]) for entry in sequence_collection.scaffolds])
            haplotype_selected_sequence_dict = dict([(entry, entry) for entry in sequence_collection.scaffolds])

        final_haplotype_set = (set(haplotype_selected_sequence_dict.keys()) & whitelist) if whitelist else set(haplotype_selected_sequence_dict.keys())

        with open(output_file, "w") as out_fd:
            #out_fd.write("#NEXUS\nBEGIN TAXA;\nDIMENSIONS\nNTAX = %i;\nTAXLABELS\n%s\n;\nEND;\n\n" % (len(haplotype_selected_sequence_dict),
            #                                                                                          "\n".join(haplotype_selected_sequence_dict.keys())))
            out_fd.write("#NEXUS\n\n")
            out_fd.write("BEGIN DATA;\n\tDIMENSIONS NTAX=%i NCHAR=%i;\n\tFORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. ;\n" % (len(final_haplotype_set),
                                                                                                                                 alignment_len))
            out_fd.write("\tMATRIX\n")

            for haplotype_id in final_haplotype_set:
                out_fd.write("\t\t%s %s\n" % (haplotype_id, sequence_collection.records[haplotype_selected_sequence_dict[haplotype_id]]))
            out_fd.write("\t;\nEND;\n\n")

            if not traits_df.empty:
                traits_number = len(traits_df.columns)
                out_fd.write("BEGIN TRAITS;\n\tDimensions NTRAITS={0};\n\tFormat labels=yes missing=? separator=Comma;\n".format(traits_number))
                out_fd.write("\tTraitLabels {0};\n".format(" ".join(traits_df.columns)))
                out_fd.write("\tMATRIX\n")
                for haplotype_id in final_haplotype_set:
                    out_fd.write("\t\t%s %s\n" % (haplotype_id,
                                                  ",".join(map(str, traits_df.loc[haplotype_id])) if haplotype_id in traits_df.index else ("0," * traits_number)[:-1]))
            else:
                out_fd.write("BEGIN TRAITS;\n\tDimensions NTRAITS=1;\n\tFormat labels=yes missing=? separator=Comma;\n")
                out_fd.write("\tTraitLabels Area;\n")
                out_fd.write("\tMATRIX\n")
                for haplotype_id in final_haplotype_set:
                    out_fd.write("\t\t%s %i\n" % (haplotype_id, len(haplotype_dict[haplotype_id])))
            out_fd.write("\t;\nEND;\n\n")

    @staticmethod
    def count_haplotype_traits(haplotype_traits, output_file):
        traits_df = pd.read_csv(haplotype_traits, sep="\t", index_col=0, header=0)
        count_series = traits_df.groupby(by=["haplotype", "location"]).size()
        haplotypes = count_series.index.get_level_values(level=0).unique()
        locations = count_series.index.get_level_values(level=1).unique()
        count_table_df = pd.DataFrame(0, index=haplotypes,
                                      columns=locations)
        for hap in haplotypes:
            for loc in locations:
                if (hap, loc) in count_series.index:
                    count_table_df.loc[hap, loc] = count_series.loc[(hap,loc)]
        #print(count_table_df)
        count_table_df.to_csv(output_file, sep="\t", header=True, index=True)


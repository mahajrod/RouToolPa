#!/usr/bin/env python
"""
Sequence Parser Module
"""
__author__ = 'Sergei F. Kliver'

import re
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import pandas as pd

from RouToolPa.Routines.Sequence import FileRoutines
from RouToolPa.Parsers.GFF import CollectionGFF


class CollectionSequence(FileRoutines):

    def __init__(self, in_file=None, records=None, description=None,
                 external_description=None,
                 format="fasta",
                 parsing_mode="parse", black_list=(), white_list=(),
                 masking=None, masking_file=None, masking_filetype="gff",
                 verbose=False, seq_expression=None, get_stats=False):
        self.formats = ["fasta"]
        self.format = format
        self.parsing_mode = parsing_mode
        self.seq_file = in_file
        self.seq_file_format = format
        self.white_list = white_list
        self.black_list = black_list
        self.seq_file = in_file
        self.description = OrderedDict()

        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode,
                      black_list=black_list,  white_list=white_list, verbose=verbose,
                      seq_expression=seq_expression)

        elif records is None:
            self.records = OrderedDict()
            self.description = description if description else OrderedDict()
        else:
            self.records = records
            self.description = description if description else OrderedDict()

        self.external_description = external_description if external_description else OrderedDict()

        if masking_file:
            print("Parsing masking...")
            self.masking = CollectionGFF(masking_file, format=masking_filetype, parsing_mode="only_coordinates",
                                         black_list=black_list, white_list=white_list)
        else:
            self.masking = masking

        self.seq_lengths = None   # None or pandas dataframe with seq_id as index
        self.length = 0
        self.scaffolds = None
        self.gaps = None          # None or pandas dataframe with seq_id as index
        self.stats = None         # None or pandas dataframe with seq_id as index

        if get_stats:
            self.get_stats_and_features(count_gaps=False, sort=False, min_gap_length=1)

    def sequence_generator(self, sequence_file, format="fasta", black_list=(), white_list=(), verbose=False):

        if format == "fasta":
            with self.metaopen(sequence_file, "r") as seq_fd:
                seq_id = None
                description = None
                seq = ""
                for line in seq_fd:
                    if line[0] == ">":
                        if seq_id and (seq_id not in black_list):
                            if (not white_list) or (seq_id in white_list):
                                if verbose:
                                    print("Parsing %s" % seq_id)
                                yield seq_id, description, seq
                        tmp = line[1:].strip().split(None, 1)
                        seq_id, description = tmp if len(tmp) == 2 else (tmp[0], "")

                        seq = ""
                    else:
                        seq += line[:-1]
                else:
                    if seq_id and (seq_id not in black_list):
                        if (not white_list) or (seq_id in white_list):
                            if verbose:
                                print("Parsing %s" % seq_id)
                            yield seq_id, description, seq

    def sequence_generator_with_expression(self, sequence_file, seq_expression, format="fasta",
                                           black_list=(), white_list=(), verbose=False,):
        if format == "fasta":
            with self.metaopen(sequence_file, "r") as seq_fd:
                seq_id = None
                description = None
                seq = ""
                for line in seq_fd:
                    if line[0] == ">":
                        if seq_id and (seq_id not in black_list):
                            if (not white_list) or (seq_id in white_list):
                                if verbose:
                                    print("Parsing %s" % seq_id)
                                yield seq_id, description, seq
                        tmp = line[1:].strip().split(None, 1)
                        seq_id, description = tmp if len(tmp) == 2 else (tmp[0], "")

                        seq = ""
                    else:
                        seq += line[:-1]
                else:
                    if seq_id and (seq_id not in black_list):
                        if (not white_list) or (seq_id in white_list):
                            if verbose:
                                print("Parsing %s" % seq_id)
                            yield seq_id, description, seq_expression(seq)

    def sequence_tuple_generator(self):
        if self.parsing_mode == "parse":
            for scaffold_id in self.scaffolds:
                yield scaffold_id, self.records[scaffold_id]

    def reset_seq_generator(self):
        self.records = self.sequence_generator(self.seq_file, format=self.seq_file_format,
                                               black_list=self.black_list,  white_list=self.white_list)

    def read(self, seq_file, format="fasta", parsing_mode="generator", black_list=(), white_list=(),
             verbose=False, seq_expression=None):
        if format not in self.formats:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        if parsing_mode == "generator":
            print("Creating sequence generator...")
            if seq_expression:
                self.records = self.sequence_generator_with_expression(seq_file,
                                                                       seq_expression,
                                                                       format=format,
                                                                       black_list=black_list,
                                                                       white_list=white_list,
                                                                       verbose=verbose)
            else:
                self.records = self.sequence_generator(seq_file,
                                                       format=format,
                                                       black_list=black_list,
                                                       white_list=white_list,
                                                       verbose=verbose)
        elif parsing_mode == "parse":
            print("Parsing sequences...")
            self.records = OrderedDict()
            for seq_id, description, seq in self.sequence_generator(seq_file, format=format,
                                                                    black_list=black_list,
                                                                    white_list=white_list):
                self.records[seq_id] = seq
                self.description[seq_id] = description

            return self.records

    def __len__(self):
        """
        :return: length of genome or None if genome was not parsed yet
        """
        return self.length

    def scaffold_number(self):
        """
        :return: number of regions/scaffolds/chromosomes in genome
        """
        return len(self.scaffolds)

    def get_stats_and_features(self, thresholds_list=(0, 500, 1000),
                               count_gaps=True, sort=True, min_gap_length=1):
        length_list = []
        gaps_list = []
        if self.parsing_mode == "generator":
            for seq_id, description, seq in self.records:
                length_list.append([seq_id, len(seq)])
                if count_gaps:
                    gaps_list.append(self.find_gaps_in_seq(seq, seq_id, min_gap_length=min_gap_length))

            self.reset_seq_generator()
        else:
            for seq_id in self.records:
                length_list.append([seq_id, len(self.records[seq_id])])
                if count_gaps:
                    gaps_list.append(self.find_gaps_in_seq(self.records[seq_id], seq_id, min_gap_length=min_gap_length))

        if count_gaps:
            self.gaps = CollectionGFF(records=pd.concat(gaps_list), format="gff", parsing_mode="only_coordinates",
                                      black_list=self.black_list, white_list=self.white_list
                                      )

            self.gaps.records.sort_values(by=["scaffold", "start", "end"])

        self.seq_lengths = pd.DataFrame.from_records(length_list, columns=("scaffold", "length"), index="scaffold").sort_values(by=["length", "scaffold"], ascending=(False, True))
        self.length = np.sum(self.seq_lengths["length"])
        self.scaffolds = self.seq_lengths.index.values

        self.get_length_stats(thresholds_list=thresholds_list, count_gaps=count_gaps)

        return self.stats

    def get_length_stats(self, thresholds_list=(0, 500, 1000), count_gaps=True):
        stats = OrderedDict()
        for threshold in thresholds_list:
            stats[threshold] = OrderedDict()
            lengths_df = self.seq_lengths[self.seq_lengths["length"] >= threshold].copy()
            lengths_df["cum_length"] = lengths_df["length"].cumsum()
            half_length = float(lengths_df["cum_length"][-1]) / 2
            lengths_df["cumlen_longer"] = lengths_df["cum_length"] >= half_length
            middle_element_index = lengths_df["cumlen_longer"].idxmax()
            L50 = lengths_df.index.get_loc(middle_element_index) + 1
            N50 = lengths_df["length"][middle_element_index]

            stats[threshold]["Total length"] = lengths_df["length"].sum()
            stats[threshold]["Total scaffolds"] = len(lengths_df["length"])
            stats[threshold]["Longest scaffold"] = lengths_df["length"][0]
            stats[threshold]["L50"] = L50
            stats[threshold]["N50"] = N50
            if count_gaps:
                gaps_df = self.gaps.records[self.gaps.records.index.isin(lengths_df.index)]
                stats[threshold]["Ns"] = sum(gaps_df["end"] - gaps_df["start"])

        stats = pd.DataFrame.from_dict(stats)
        self.stats = stats
        return self.stats

    @staticmethod
    def find_gaps_in_seq(sequence, seq_id=None, min_gap_length=1):
        """
        Finds gaps (N) in sequence
        :return: None
        """
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        gaps_list = []
        gaps = gap_reg_exp.finditer(sequence)  # iterator with
        for match in gaps:
            if (match.end() - match.start()) >= min_gap_length:
                gaps_list.append([match.start(), match.end()])
        gaps_list = pd.DataFrame(gaps_list, columns=("start", "end"))
        if seq_id:
            # add seq id as index
            gaps_list.index = pd.Index([seq_id] * len(gaps_list), name="scaffold")
        return gaps_list

    def remove_small_gaps(self, min_gap_length):
        self.gaps = self.gaps[(self.gaps['end'] - self.gaps['start']) >= min_gap_length]

    def get_merged_gaps_and_masking(self):
        if self.masking and self.gaps:
            merged = self.masking + self.gaps
            merged.collapse_records(sort=False, verbose=True)  # sorting is called during addition
            return merged
        elif self.masking:
            return self.masking
        elif self.gaps:
            return self.gaps
        else:
            return None

    def write_splited(self, out_dir, expression=None, max_symbols_per_line=60):
        if self.parsing_mode == "parse":
            for seq_id in self.records:
                if expression:
                    if not expression(seq_id, self.records[seq_id]):
                        continue
                with self.metaopen("%s/%s.fasta" % (out_dir, seq_id), "w") as out_fd:
                    out_fd.write(">%s\n" % seq_id if seq_id not in self.description else ">%s %s\n" % (seq_id, self.description[seq_id]))
                    length = self.seq_lengths[seq_id][0] if self.seq_lengths else len(self.records[seq_id])
                    line_number = length // max_symbols_per_line
                    index = 0
                    while index < line_number:
                        out_fd.write(self.records[seq_id][
                                     index * max_symbols_per_line:(index + 1) * max_symbols_per_line] + "\n")
                        index += 1
                    if line_number * max_symbols_per_line < length:
                        out_fd.write(self.records[seq_id][index * max_symbols_per_line:] + "\n")
        else:
            raise ValueError("ERROR!!! Writing was implemented only for parsing mode yet!")

    def write(self, outfile, expression=None, max_symbols_per_line=60):
        if self.parsing_mode == "parse":
            with self.metaopen(outfile, "w") as out_fd:
                for seq_id in self.records:
                    if expression:
                        if not expression(seq_id, self.records[seq_id]):
                            continue
                    out_fd.write(">%s\n" % seq_id if seq_id not in self.description else ">%s %s\n" % (seq_id, self.description[seq_id]))
                    length = self.seq_lengths[seq_id][0] if self.seq_lengths else len(self.records[seq_id])
                    line_number = length // max_symbols_per_line
                    index = 0
                    while index < line_number:
                        out_fd.write(self.records[seq_id][index*max_symbols_per_line:(index+1)*max_symbols_per_line] + "\n")
                        index += 1
                    if line_number * max_symbols_per_line < length:
                        out_fd.write(self.records[seq_id][index*max_symbols_per_line:] + "\n")
        if self.parsing_mode == "generator":
            with self.metaopen(outfile, "w") as out_fd:
                for seq_id, description, seq in self.records:
                    if expression:
                        if not expression(seq_id, seq):
                            continue
                    out_fd.write(">%s\n" % seq_id if not description else ">%s %s\n" % (seq_id, description))
                    length = len(seq)
                    line_number = length // max_symbols_per_line
                    index = 0
                    while index < line_number:
                        out_fd.write(seq[index*max_symbols_per_line:(index+1)*max_symbols_per_line] + "\n")
                        index += 1
                    if line_number * max_symbols_per_line < length:
                        out_fd.write(seq[index*max_symbols_per_line:] + "\n")
        else:
            raise ValueError("ERROR!!! Writing was implemented only for parsing mode yet!")

    def write_by_syn(self, outfile, syn_dict, max_symbols_per_line=60, absent_symbol="."):
        """
        :param outfile: output fasta file
        :param syn_dict: dictionary with synonyms to be used for sequence ids.
               Note: sequence form collection will be written as many times as it appears in values of syn_dict
        :param max_symbols_per_line: Maximum symbols per line in output fasta
        :param absent_symbol:
        :return:
        """
        if self.parsing_mode == "parse":
            with self.metaopen(outfile, "w") as out_fd:
                for seq_id in syn_dict:
                    if syn_dict[seq_id] == absent_symbol:
                        continue

                    out_fd.write(">%s\n" % seq_id if syn_dict[seq_id] not in self.description else ">%s %s\n" % (seq_id, self.description[syn_dict[seq_id]]))

                    length = self.seq_lengths[syn_dict[seq_id]][0] if self.seq_lengths else len(self.records[syn_dict[seq_id]])
                    line_number = length // max_symbols_per_line
                    index = 0
                    while index < line_number:
                        out_fd.write(self.records[syn_dict[seq_id]][index*max_symbols_per_line:(index+1)*max_symbols_per_line] + "\n")
                        index += 1
                    if line_number * max_symbols_per_line < length:
                        out_fd.write(self.records[syn_dict[seq_id]][index*max_symbols_per_line:] + "\n")

        else:
            raise ValueError("ERROR!!! Writing was implemented only for parsing mode yet!")

    def convert_rna_to_dna(self):
        if self.parsing_mode == "parse":
            for seq_id in self.scaffolds:
                self.records[seq_id] = self.records[seq_id].replace("U", "T").replace("u", "t")
        else:
            raise ValueError("ERROR!!! Convertation from RNA to DNA was implemented only for parsing mode yet!")

    @staticmethod
    def count_window_number_in_scaffold(scaffold_length, window_size, window_step):
        if scaffold_length < window_size:
            return 0
        return int((scaffold_length - window_size)/window_step) + 1

    def count_window_number(self, window_size, window_step):

        def count(scaffold_length):
            return self.count_window_number_in_scaffold(scaffold_length, window_size, window_step)

        return self.seq_lengths.apply(count)

    @staticmethod
    def expression_unmask(s):
        return s.upper()

    def unmask(self, in_place=False, verbose=None):

        if in_place:
            collection = self
        else:
            collection = deepcopy(self)

        if collection.parsing_mode == "generator":
            print("Creating unmasked sequence generator...")

            collection.records = collection.sequence_generator_with_expression(collection.seq_file,
                                                                               collection.expression_unmask,
                                                                               format=collection.format,
                                                                               black_list=collection.black_list,
                                                                               white_list=collection.white_list,
                                                                               verbose=verbose)
        elif collection.parsing_mode == "parse":
            print("Unmasking sequences...")
            for seq_id in collection.records:
                collection.records[seq_id] = collection.records[seq_id].upper()

        if not in_place:
            return collection

    def calculate_p_distance(self):
        from RouToolPa.Routines import SequenceRoutines
        if self.parsing_mode == "generator":
            raise ValueError("ERROR!!! P distance calculation was not implemented for generator mode!")

        if self.seq_lengths is None:
            self.get_stats_and_features(count_gaps=False, sort=False)
        seq_len = self.seq_lengths["length"].unique()
        if len(seq_len) > 1:
            raise ValueError("ERROR!!! Some sequences have different length!")
        else:
            seq_len = seq_len[0]

        distance_df = pd.DataFrame(0, index=self.scaffolds, columns=self.scaffolds)
        for record_id_a in self.scaffolds:
            for record_id_b in self.scaffolds:
                if record_id_a == record_id_b:
                    continue

                distance_df.loc[record_id_a, record_id_b] = distance_df.loc[record_id_b, record_id_a] = SequenceRoutines.p_distance(self.records[record_id_a], self.records[record_id_b], seq_len)

        return distance_df




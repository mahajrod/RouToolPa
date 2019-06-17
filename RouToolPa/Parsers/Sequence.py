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

from RouToolPa.Routines.File import FileRoutines
from RouToolPa.Parsers.GFF import CollectionGFF


class CollectionSequence(FileRoutines):

    def __init__(self, in_file=None, records=None, format="fasta",
                 parsing_mode="parse", black_list=(), white_list=(),
                 masking=None, masking_file=None, masking_filetype="gff",
                 verbose=False):
        self.formats = ["fasta"]
        self.parsing_mode = parsing_mode
        self.seq_file = in_file
        self.seq_file_format = format
        self.white_list = white_list
        self.black_list = black_list

        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode,
                      black_list=black_list,  white_list=white_list, verbose=verbose)
        elif records is None:
            self.records = OrderedDict()
        else:
            self.records = records

        if masking_file:
            print("Parsing masking...")
            self.masking = CollectionGFF(masking_file, format=masking_filetype, parsing_mode="only_coordinates",
                                         black_list=black_list, white_list=white_list)
        else:
            self.masking = masking

        self.seq_lengths = None # None or pandas dataframe with seq_id as index
        self.length = 0
        self.scaffolds = None
        self.gaps = None          # None or pandas dataframe with seq_id as index

    def sequence_generator(self, sequence_file, format="fasta", black_list=(), white_list=(), verbose=False):

        if format == "fasta":
            with self.metaopen(sequence_file, "r") as seq_fd:
                seq_id = None
                seq = ""
                for line in seq_fd:
                    if line[0] == ">":
                        if seq_id and (seq_id not in black_list):
                            if (not white_list) or (seq_id in white_list):
                                if verbose:
                                    print("Parsing %s" % seq_id)
                                yield seq_id, seq
                        seq_id = line[1:].split()[0]
                        seq = ""
                    else:
                        seq += line[:-1]
                else:
                    if seq_id and (seq_id not in black_list):
                        if (not white_list) or (seq_id in white_list):
                            if verbose:
                                print("Parsing %s" % seq_id)
                            yield seq_id, seq

    def reset_seq_generator(self):
        self.records = self.sequence_generator(self.seq_file, format=self.seq_file_format,
                                               black_list=self.black_list,  white_list=self.white_list)

    def read(self, seq_file, format="fasta", parsing_mode="generator", black_list=(), white_list=(),
             verbose=False):
        if format not in self.formats:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        if parsing_mode == "generator":
            print("Creating sequence generator...")
            self.records = self.sequence_generator(seq_file, format=format,
                                                   black_list=black_list,  white_list=white_list,
                                                   verbose=verbose)
        elif parsing_mode == "parse":
            print("Parsing sequences...")
            self.records = OrderedDict()
            for seq_id, seq in self.sequence_generator(seq_file, format=format,
                                                       black_list=black_list,
                                                       white_list=white_list):
                self.records[seq_id] = seq

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

    def get_stats_and_features(self, count_gaps=True, sort="True", min_gap_length=1):
        length_list = []
        gaps_list = []
        if self.parsing_mode == "generator":
            for seq_id, seq in self.records:
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

        self.seq_lengths = pd.DataFrame.from_records(length_list, columns=("scaffold", "length"), index="scaffold")
        if sort:
            self.seq_lengths.sort_values(by=["length", "scaffold"])
        self.length = np.sum(self.seq_lengths["length"])
        self.scaffolds = self.seq_lengths.index.values

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
            merged.collapse_records(sort=False, verbose=True) # sorting is called during addition
            return merged
        elif self.masking:
            return self.masking
        elif self.gaps:
            return self.gaps
        else:
            return None

    def write(self, outfile, expression=None, max_symbols_per_line=60):
        if self.parsing_mode == "parse":
            with self.metaopen(outfile, "w") as out_fd:
                for seq_id in self.records:
                    if expression:
                        if not expression(seq_id, self.records[seq_id]):
                            continue
                    out_fd.write(">%s\n" % seq_id)
                    length = self.seq_lengths[seq_id][0] if self.seq_lengths else len(self.records[seq_id])
                    line_number = length // max_symbols_per_line
                    index = 0
                    while index < line_number:
                        out_fd.write(self.records[seq_id][index*max_symbols_per_line:(index+1)*max_symbols_per_line] + "\n")
                        index += 1
                    if line_number * max_symbols_per_line < length:
                        out_fd.write(self.records[seq_id][index*max_symbols_per_line:-1] + "\n")

        else:
            raise ValueError("ERROR!!! Writing was implemented only for parsing mode yet!")

    def write_by_syn(self, outfile, syn_dict, max_symbols_per_line=60, absent_symbol="."):
        """
        :param outfile: output fasta file
        :param syn_dict: dictionary with synonyms to be used for sequence ids.
               Note: sequence form collection will be written as many times as it appears in values of syn_dict
        :param max_symbols_per_line: Maximum symbols per line in output fasta
        :return:
        """
        if self.parsing_mode == "parse":
            with self.metaopen(outfile, "w") as out_fd:
                for seq_id in syn_dict:
                    if syn_dict[seq_id] == ".":
                        continue
                    out_fd.write(">%s\n" % seq_id)
                    length = self.seq_lengths[syn_dict[seq_id]][0] if self.seq_lengths else len(self.records[syn_dict[seq_id]])
                    line_number = length // max_symbols_per_line
                    index = 0
                    while index < line_number:
                        out_fd.write(self.records[syn_dict[seq_id]][index*max_symbols_per_line:(index+1)*max_symbols_per_line] + "\n")
                        index += 1
                    if line_number * max_symbols_per_line < length:
                        out_fd.write(self.records[syn_dict[seq_id]][index*max_symbols_per_line:-1] + "\n")

        else:
            raise ValueError("ERROR!!! Writing was implemented only for parsing mode yet!")


    @staticmethod
    def count_window_number_in_scaffold(scaffold_length, window_size, window_step):
        if scaffold_length < window_size:
            return 0
        return int((scaffold_length - window_size)/window_step) + 1

    def count_window_number(self, window_size, window_step):

        def count(scaffold_length):
            return self.count_window_number_in_scaffold(scaffold_length, window_size, window_step)

        return self.length.apply(count)

#!/usr/bin/env python
"""
Last output parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'
import datetime
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import pandas as pd

import RouToolPa.Formats.AlignmentFormats as AlignmentFormats


class CollectionBLAST:

    def __init__(self, in_file=None, records=None, format="tab6", parsing_mode="complete",
                 target_black_list=(), target_white_list=(),
                 query_black_list=(), query_white_list=(),
                 target_syn_dict=None, query_syn_dict=None,
                 min_target_hit_len=None, min_query_hit_len=None,
                 min_target_len=None, min_query_len=None):

        self.formats = ["tab6", "busco_tblastn"]
        self.parsing_parameters = {"tab6": {"complete": {
                                                   "col_names": ["query_id",
                                                                 "target_id",
                                                                 "identity,%%",
                                                                 "length",
                                                                 "mismatch",
                                                                 "gapopen",
                                                                 "query_start",
                                                                 "query_end",
                                                                 "target_start",
                                                                 "target_end",
                                                                 "evalue",
                                                                 "bitscore"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                                                   "index_cols": ["qseqid", "sseqid"],
                                                   "converters": {
                                                                  "qseqid":     str,
                                                                  "sseqid":     str,
                                                                  "pident":     np.float32,
                                                                  "length":     np.int64,
                                                                  "mismatch":   np.int32,
                                                                  "gapopen":    np.int32,
                                                                  "qstart":     np.int64,
                                                                  "qend":       np.int64,
                                                                  "sstart":     np.int64,
                                                                  "send":       np.int64,
                                                                  "evalue":     np.float64,
                                                                  "bitscore":   np.float32
                                                                 },
                                                   "col_name_indexes": {
                                                                        "query_id":         0,
                                                                        "targetid":         1,
                                                                        "identity,%%":      2,
                                                                        "length":           3,
                                                                        "mismatch":         4,
                                                                        "gapopen":          5,
                                                                        "query_start":      6,
                                                                        "query_end":        7,
                                                                        "target_start":     8,
                                                                        "target_end":       9,
                                                                        "evalue":           10,
                                                                        "bitscore":         11
                                                                        },
                                                   "original_col_names": ["qseqid",
                                                                          "sseqid",
                                                                          "pident",
                                                                          "length",
                                                                          "mismatch",
                                                                          "gapopen",
                                                                          "qstart",
                                                                          "qend",
                                                                          "sstart",
                                                                          "send",
                                                                          "evalue",
                                                                          "bitscore"],
                                                   "original_col_name_indexes": {
                                                                                 "qseqid":     0,
                                                                                 "sseqid":     1,
                                                                                 "pident":     2,
                                                                                 "length":     3,
                                                                                 "mismatch":   4,
                                                                                 "gapopen":    5,
                                                                                 "qstart":     6,
                                                                                 "qend":       7,
                                                                                 "sstart":     8,
                                                                                 "send":       9,
                                                                                 "evalue":     10,
                                                                                 "bitscore":   11
                                                                                 },
                                                   },

                                           },
                                 "busco_tblastn": {"complete": {
                                                               "col_names": ["query_id",
                                                                             "target_id",
                                                                             "identity,%%",
                                                                             "length",
                                                                             "mismatch",
                                                                             "gapopen",
                                                                             "query_start",
                                                                             "query_end",
                                                                             "target_start",
                                                                             "target_end",
                                                                             "evalue",
                                                                             "bitscore"],
                                                               "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                                                               "index_cols": ["qseqid", "sseqid"],
                                                               "converters": {
                                                                              "qseqid":     str,
                                                                              "sseqid":     str,
                                                                              "pident":     np.float32,
                                                                              "length":     np.int64,
                                                                              "mismatch":   np.int32,
                                                                              "gapopen":    np.int32,
                                                                              "qstart":     np.int64,
                                                                              "qend":       np.int64,
                                                                              "sstart":     np.int64,
                                                                              "send":       np.int64,
                                                                              "evalue":     np.float64,
                                                                              "bitscore":   np.float32
                                                                             },
                                                               "col_name_indexes": {
                                                                                    "query_id":         0,
                                                                                    "targetid":         1,
                                                                                    "identity,%%":      2,
                                                                                    "length":           3,
                                                                                    "mismatch":         4,
                                                                                    "gapopen":          5,
                                                                                    "query_start":      6,
                                                                                    "query_end":        7,
                                                                                    "target_start":     8,
                                                                                    "target_end":       9,
                                                                                    "evalue":           10,
                                                                                    "bitscore":         11
                                                                                    },
                                                               "original_col_names": ["qseqid",
                                                                                      "sseqid",
                                                                                      "pident",
                                                                                      "length",
                                                                                      "mismatch",
                                                                                      "gapopen",
                                                                                      "qstart",
                                                                                      "qend",
                                                                                      "sstart",
                                                                                      "send",
                                                                                      "evalue",
                                                                                      "bitscore"],
                                                               "original_col_name_indexes": {
                                                                                             "qseqid":     0,
                                                                                             "sseqid":     1,
                                                                                             "pident":     2,
                                                                                             "length":     3,
                                                                                             "mismatch":   4,
                                                                                             "gapopen":    5,
                                                                                             "qstart":     6,
                                                                                             "qend":       7,
                                                                                             "sstart":     8,
                                                                                             "send":       9,
                                                                                             "evalue":     10,
                                                                                             "bitscore":   11
                                                                                             },
                                                                },

                                           },
                                    }
        self.parsing_mode = parsing_mode

        self.format = format
        self.target_black_list = target_black_list
        self.target_white_list = target_white_list
        self.target_syn_dict = target_syn_dict

        self.query_black_list = query_black_list
        self.query_white_list = query_white_list
        self.query_syn_dict = query_syn_dict
        # attributes type conversion parameters

        self.converters = OrderedDict()
        self.pandas_int_type_correspondence = OrderedDict({
                                                           "Int8":  np.float16,
                                                           "Int16": np.float16,
                                                           "Int32": np.float32,
                                                           "Int64": np.float64,
                                                           })

        # init aliases
        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]
        self.target_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_id"]
        self.target_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_start"]
        #self.target_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_hit_len"]
        #self.target_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_strand"]

        self.query_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_id"]
        self.query_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_start"]
        #self.query_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_hit_len"]
        #self.query_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_strand"]

        #
        self.target_scaffold_list = None
        self.query_scaffold_list = None
        self.target_scaffold_lengths = None
        self.query_scaffold_lengths = None

        if in_file:
            self.read(in_file,
                      format=format,
                      parsing_mode=parsing_mode,
                      target_black_list=target_black_list,
                      target_white_list=target_white_list,
                      query_black_list=query_black_list,
                      query_white_list=query_white_list,
                      min_target_hit_len=min_target_hit_len,
                      min_query_hit_len=min_query_hit_len,
                      min_target_len=min_target_len,
                      min_query_len=min_query_len)

        else:
            self.records = records

    def read(self, in_file,
             format="tab6", parsing_mode="all",
             target_black_list=(), target_white_list=(),
             query_black_list=(), query_white_list=(),
             min_target_hit_len=None, min_query_hit_len=None,
             min_target_len=None, min_query_len=None):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        print("%s\tReading input..." % str(datetime.datetime.now()))
        self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                   comment="#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])

        print("%s\tReading input finished..." % str(datetime.datetime.now()))
        print("%s\tFiltering..." % str(datetime.datetime.now()))
        if target_white_list or target_black_list:
            target_scaffolds_to_keep = self.get_filtered_entry_list(self.records["target_id"].tolist(),
                                                                    entry_black_list=target_black_list,
                                                                    entry_white_list=target_white_list)
            #print target_scaffolds_to_keep
            self.records = self.records[self.records["target_id"].isin(target_scaffolds_to_keep)]

        if query_white_list or query_black_list:
            query_scaffolds_to_keep = self.get_filtered_entry_list(self.records["query_id"].tolist(),
                                                                   entry_black_list=query_black_list,
                                                                   entry_white_list=query_white_list)
            #print query_scaffolds_to_keep
            self.records = self.records[self.records["query_id"].isin(query_scaffolds_to_keep)]
        if self.target_syn_dict:
            self.records["target_id"].replace(self.target_syn_dict, inplace=True)
        if self.query_syn_dict:
            self.records["query_id"].replace(self.query_syn_dict, inplace=True)
        # retain only automatically generated index by row number
        #self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
        #                                               names=("scaffold", "row"))
        print("%s\tFiltering finished..." % str(datetime.datetime.now()))

        if min_target_len and min_query_len:
            self.records = self.records[(self.records["target_len"] >= min_target_len) & (self.records["query_len"] >= min_query_len)]
        elif min_target_len:
            self.records = self.records[self.records["target_len"] >= min_target_len]
        elif min_query_len:
            self.records = self.records[self.records["query_len"] >= min_query_len]
            
        if min_target_hit_len and min_query_hit_len:
            self.records = self.records[(self.records["target_hit_len"] >= min_target_hit_len) & (self.records["query_hit_len"] >= min_query_hit_len)]
        elif min_target_hit_len:
            self.records = self.records[self.records["target_hit_len"] >= min_target_hit_len]
        elif min_query_hit_len:
            self.records = self.records[self.records["query_hit_len"] >= min_query_hit_len]

        if "target_len" in self.parsing_parameters[self.format][self.parsing_mode]["col_names"]:
            self.target_scaffold_lengths = self.records[["target_id", "target_len"]].drop_duplicates()
            self.target_scaffold_lengths.columns = ["id", "length"]
            self.target_scaffold_lengths.set_index("id", inplace=True)
            self.target_scaffold_list = self.target_scaffold_lengths.index.tolist()
        if "query_len" in self.parsing_parameters[self.format][self.parsing_mode]["col_names"]:
            self.query_scaffold_lengths = self.records[["query_id", "query_len"]].drop_duplicates()
            self.query_scaffold_lengths.columns = ["id", "length"]
            self.query_scaffold_lengths.set_index("id", inplace=True)
            self.query_scaffold_list = self.query_scaffold_lengths.index.tolist()

        retained_columns = deepcopy(self.parsing_parameters[self.format][self.parsing_mode]["col_names"])
        for entry in "target_len", "query_len":
            if entry in retained_columns:
                retained_columns.remove(entry)

        self.records = self.records[retained_columns]
    """
    def sort(self, inplace=False,
             sorting_order=("target_id", "query_id", "target_start", "target_hit_len", "query_start", "query_hit_len")):
        return self.records.sort_values(by=sorting_order, inplace=inplace)
    """

    @staticmethod
    def get_filtered_entry_list(entry_list,
                                entry_black_list=[],
                                sort_entries=False,
                                entry_ordered_list=None,
                                entry_white_list=[]):
        white_set = set(entry_white_list)
        black_set = set(entry_black_list)
        entry_set = set(entry_list)

        if white_set:
            entry_set = entry_set & white_set
        if black_set:
            entry_set = entry_set - black_set

        filtered_entry_list = list(entry_set)
        if sort_entries:
            filtered_entry_list.sort()

        final_entry_list = []

        if entry_ordered_list:
            for entry in entry_ordered_list:
                if entry in filtered_entry_list:
                    final_entry_list.append(entry)
                    filtered_entry_list.remove(entry)
                else:
                    print("WARNING!!!Entry(%s) from order list is absent in list of entries!" % entry)
            return final_entry_list + filtered_entry_list
        else:
            return filtered_entry_list

    def write(self, output, separator="\t", header=False):
        self.records.to_csv(output,
                            sep=separator,
                            header=header,
                            index=False)

    def rename_target_ids(self, syn_dict):
        self.records[["target_id"]].replace(syn_dict, inplace=True)

    def rename_query_ids(self, syn_dict):
        self.records[["query_id"]].replace(syn_dict, inplace=True)


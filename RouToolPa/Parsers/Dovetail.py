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


class CollectionLinks:

    def __init__(self, in_file=None, records=None, format="tab", parsing_mode="complete",
                 target_black_list=(), target_white_list=(),
                 query_black_list=(), query_white_list=(),
                 target_syn_dict=None, query_syn_dict=None,
                 min_target_hit_len=None, min_query_hit_len=None,
                 min_target_len=None, min_query_len=None):

        self.formats = ["tab", "busco_tblastn"]
        self.parsing_parameters = {"tab": {"complete": {
                                                   "col_names": ["hic_scaffold_id",
                                                                 "raw_scaffold_id",
                                                                 "raw_scaffold_start",
                                                                 "raw_scaffold_end",
                                                                 "raw_scaffold_strand",
                                                                 "hic_scaffold_start",
                                                                 "hic_scaffold_end"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6],
                                                   "index_cols": ["hic_scaffold_id"],
                                                   "converters": {
                                                                  "hic_scaffold_id":        str,
                                                                  "raw_scaffold_id":        str,
                                                                  "raw_scaffold_start":     np.int64,
                                                                  "raw_scaffold_end":       np.int64,
                                                                  "raw_scaffold_strand":    str,
                                                                  "hic_scaffold_start":     np.int64,
                                                                  "hic_scaffold_end":       np.int64
                                                                 },
                                                   "col_name_indexes": {
                                                                        "hic_scaffold_id":      0,
                                                                        "raw_scaffold_id":      1,
                                                                        "raw_scaffold_start":   2,
                                                                        "raw_scaffold_end":     3,
                                                                        "raw_scaffold_strand":  4,
                                                                        "hic_scaffold_start":   5,
                                                                        "hic_scaffold_end":     6
                                                                        },
                                                   "original_col_names": ["hic_scaffold_id",
                                                                          "raw_scaffold_id",
                                                                          "raw_scaffold_start",
                                                                          "raw_scaffold_end",
                                                                          "raw_scaffold_strand",
                                                                          "hic_scaffold_start",
                                                                          "hic_scaffold_end"],
                                                   "original_col_name_indexes": {
                                                                                 "hic_scaffold_id":      0,
                                                                                 "raw_scaffold_id":      1,
                                                                                 "raw_scaffold_start":   2,
                                                                                 "raw_scaffold_end":     3,
                                                                                 "raw_scaffold_strand":  4,
                                                                                 "hic_scaffold_start":   5,
                                                                                 "hic_scaffold_end":     6
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
        #self.target_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_id"]
        #self.target_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_start"]
        #self.target_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_hit_len"]
        #self.target_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_strand"]

        #self.query_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_id"]
        #self.query_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_start"]
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
             format="tab", parsing_mode="all",
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

    def get_raw_coordinate(self, hic_scaffold, start, end):
        hic_scaffold_df = self.records.loc[hic_scaffold]
        hic_scaffold_df = hic_scaffold_df[(hic_scaffold_df["hic_scaffold_start"] <= start) & (hic_scaffold_df["hic_scaffold_end"] >= end)]

        if len(hic_scaffold_df["raw_scaffold_id"]) > 1:
            print("Ambigious coordinates")
            return None
        elif hic_scaffold_df.empty:
            return None

        if hic_scaffold_df["raw_scaffold_strand"][0] == "+":
            raw_start = start - hic_scaffold_df["hic_scaffold_start"][0]
            raw_end = end - hic_scaffold_df["hic_scaffold_start"][0]
        else:
            raw_start = hic_scaffold_df["raw_scaffold_end"][0] - start + hic_scaffold_df["hic_scaffold_start"][0]
            raw_end = hic_scaffold_df["raw_scaffold_end"][0] - end + hic_scaffold_df["hic_scaffold_start"][0]

        return hic_scaffold_df["raw_scaffold_id"][0], raw_start, raw_end







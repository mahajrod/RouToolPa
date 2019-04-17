#!/usr/bin/env python
"""
Last output parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'
import datetime
from copy import deepcopy
from collections import OrderedDict, Iterable

import numpy as np
import pandas as pd


class CollectionLast:

    def __init__(self, in_file=None, records=None, format="tab", parsing_mode="only_coordinates",
                 black_list=(), white_list=(), min_target_hit_len=None, min_query_hit_len=None,
                 min_target_len=None, min_query_len=None):

        self.formats = ["tab"]
        self.TAB_COLS = OrderedDict({
                                     "score": 0,
                                     "target_id":   1,
                                     "target_start": 2,
                                     "target_hit_len": 3,
                                     "target_strand": 4,
                                     "target_len": 5,
                                     "query_id": 6,
                                     "query_start": 7,
                                     "query_hit_len": 8,
                                     "query_strand": 9,
                                     "query_len": 10,
                                     "alignment": 11,
                                     "EG2": 12,
                                     "E": 13,
                                     })
        self.parsing_parameters = {"tab": {
                                           "all": {
                                                   "col_names": ["score", "target_id", "target_start", "target_hit_len",
                                                                 "target_strand", "target_len", "query_id",
                                                                 "query_start", "query_hit_len", "query_strand",
                                                                 "query_len", "alignment", "EG", "E"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                                   "index_cols": None,
                                                   "converters": {
                                                                  "score":          np.int64,
                                                                  "target_id":      str,
                                                                  "target_start":   np.int64,
                                                                  "target_hit_len": np.int64,
                                                                  "target_strand":  str,
                                                                  "target_len":     np.int64,
                                                                  "query_id":       str,
                                                                  "query_start":    np.int64,
                                                                  "query_hit_len":  np.int64,
                                                                  "query_strand":   str,
                                                                  "query_len":      np.int64,
                                                                  "alignment":      str,
                                                                  "EG":             str,
                                                                  "E":              str,
                                                                  },
                                                   "col_name_indexes": {
                                                                        "score": 0,
                                                                        "target_id":   1,
                                                                        "target_start": 2,
                                                                        "target_hit_len": 3,
                                                                        "target_strand": 4,
                                                                        "target_len": 5,
                                                                        "query_id": 6,
                                                                        "query_start": 7,
                                                                        "query_hit_len": 8,
                                                                        "query_strand": 9,
                                                                        "query_len": 10,
                                                                        "alignment": 11,
                                                                        "EG": 12,
                                                                        "E": 13,
                                                                        },
                                                   },
                                           "complete": {
                                                        "col_names": ["score", "target_id", "target_start", "target_hit_len",
                                                                      "target_strand", "target_len", "query_id",
                                                                      "query_start", "query_hit_len", "query_strand",
                                                                      "query_len", "alignment", "EG", "E"],
                                                        "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                                        "index_cols": None,
                                                        "converters": {
                                                                       "score":          np.int64,
                                                                       "target_id":      str,
                                                                       "target_start":   np.int64,
                                                                       "target_hit_len": np.int64,
                                                                       "target_strand":  str,
                                                                       "target_len":     np.int64,
                                                                       "query_id":       str,
                                                                       "query_start":    np.int64,
                                                                       "query_hit_len":  np.int64,
                                                                       "query_strand":   str,
                                                                       "query_len":      np.int64,
                                                                       "alignment":      str,
                                                                       "EG":             str,
                                                                       "E":              str,
                                                                       },
                                                        "col_name_indexes": {
                                                                             "score": 0,
                                                                             "target_id":   1,
                                                                             "target_start": 2,
                                                                             "target_hit_len": 3,
                                                                             "target_strand": 4,
                                                                             "target_len": 5,
                                                                             "query_id": 6,
                                                                             "query_start": 7,
                                                                             "query_hit_len": 8,
                                                                             "query_strand": 9,
                                                                             "query_len": 10,
                                                                             "alignment": 11,
                                                                             "EG": 12,
                                                                             "E": 13,
                                                                             },
                                                   },
                                           },

                                   }
        self.parsing_mode = parsing_mode
        self.alignment_parsing_modes = ["all", "complete"]
        #self.attributes_parsing_modes = ["complete", "coord_and_attr"]
        self.format = format
        self.black_list = black_list
        self.white_list = white_list

        # attributes type conversion parameters

        self.converters = OrderedDict()
        self.pandas_int_type_correspondence = OrderedDict({
                                                           "Int8":  np.float16,
                                                           "Int16": np.float16,
                                                           "Int32": np.float32,
                                                           "Int64": np.float64,
                                                           })

        # init aliases
        self.record_id_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["scaffold"]
        self.record_start_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["start"]
        self.record_end_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["end"]

        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]

        #
        self.target_scaffold_list = None
        self.query_scaffold_list = None
        self.target_scaffold_lengths = None
        self.query_scaffold_lengths = None

        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode, black_list=black_list, white_list=white_list,
                      min_target_hit_len=min_target_hit_len, min_query_hit_len=min_query_hit_len, 
                      min_target_len=min_target_len, min_query_len=min_query_len)

        else:
            self.records = records

    def read(self, in_file, format="tab", parsing_mode="only_coordinates",
             black_list=(), white_list=(), min_target_hit_len=None, min_query_hit_len=None,
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

        if white_list or black_list:
            scaffolds_to_keep = self.get_filtered_entry_list(self.records.index, entry_black_list=black_list,
                                                             entry_white_list=white_list)
            self.records = self.records[self.records.index.get_level_values('scaffold').isin(scaffolds_to_keep)]

        # retain only automatically generated index by row number
        #self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
        #                                               names=("scaffold", "row"))
        print("%s\tReading input finished..." % str(datetime.datetime.now()))

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

        self.target_scaffold_lengths = self.records[["target_id", "target_len"]].drop_duplicates()
        self.target_scaffold_lengths.set_index("target_id")
        self.query_scaffold_lengths = self.records[["query_id", "query_len"]].drop_duplicates()
        self.query_scaffold_lengths.set_index("query_id")
        self.target_scaffold_list = self.target_scaffold_lengths.index.to_list()
        self.query_scaffold_list = self.query_scaffold_lengths.index.to_list()

        retained_columns = deepcopy(self.parsing_parameters[self.format][self.parsing_mode]["col_names"])
        for entry in "target_len", "query_len":
            retained_columns.remove(entry)

        self.records = self.records[retained_columns]

        if parsing_mode == "complete":
            self.records["EG2"] = map(lambda s: np.float32(s.split("=")[1]), list(self.records["EG2"]))
            self.records["E"] = map(lambda s: np.float32(s.split("=")[1]), list(self.records["E"]))

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

    @staticmethod
    def last_alignment_to_(last_aln_str):

        splited_string = map(lambda s: s.split[":"], last_aln_str.split(","))

    def get_coordinates_from_(self):
        pass


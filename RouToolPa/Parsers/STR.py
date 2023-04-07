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


class CollectionSTR:

    def __init__(self, in_file=None, records=None, format="str", parsing_mode="all",
                 black_list=(), white_list=(),
                 syn_dict=None):

        self.formats = ["str", "filtered_str"]
        self.TAB_COLS = AlignmentFormats.ALN_FMT_COLS["tab"]
        self.parsing_parameters = {
                                   "str": {
                                           "coordinates_only": {
                                                   "col_names": ["scaffold_id", "start", "end"],
                                                   "cols":      [1, 18, 19],
                                                   "index_cols": None,
                                                   "converters": {
                                                                  "scaffold_id":      str,
                                                                  "start":            np.int64,
                                                                  "end":              np.int64,
                                                                  },
                                                   "col_name_indexes": {
                                                                        "scaffold_id":   0,
                                                                        "start":         1,
                                                                        "end":           2,
                                                                        },
                                                                 },

                                           "all": {
                                                   "col_names": ["amplicon_len",
                                                                 "scaffold_id",
                                                                 "forward_primer_id",
                                                                 "forward_primer_seq",
                                                                 "forward_primer_degeneracies",
                                                                 "forward_primer_mismatches",
                                                                 "reverse_primer_id",
                                                                 "reverse_primer_seq",
                                                                 "reverse_primer_degeneracies",
                                                                 "reverse_primer_mismatches",
                                                                 "rev_comp_reverse_primer",
                                                                 "probe_id",
                                                                 "probe_seq",
                                                                 "revcomp_probe",
                                                                 "probe_degeneracies",
                                                                 "probe_mismatches",
                                                                 "probe_start_on_amplicon",
                                                                 "probe_strand",
                                                                 "start",
                                                                 "end",
                                                                 "full_hit_id",
                                                                 "scaffold_len",
                                                                 "amplicon_seq",
                                                                 "primary_tag",
                                                                 "gene",
                                                                 "product",
                                                                 "protein_id",
                                                                 "note"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                                 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
                                                   "index_cols": None,
                                                   "converters": {
                                                                  "amplicon_len":                 np.int32,
                                                                  "scaffold_id":                  str,
                                                                  "forward_primer_id":            str,
                                                                  "forward_primer_seq":           str,
                                                                  "forward_primer_degeneracies":  np.int32,
                                                                  "forward_primer_mismatches":    np.int32,
                                                                  "reverse_primer_id":            str,
                                                                  "reverse_primer_seq":           str,
                                                                  "reverse_primer_degeneracies":  np.int32,
                                                                  "reverse_primer_mismatches":    np.int32,
                                                                  "rev_comp_reverse_primer":      str,
                                                                  "probe_id":                     str,
                                                                  "probe_seq":                    str,
                                                                  "revcomp_probe":                str,
                                                                  "probe_degeneracies":           np.int32,
                                                                  "probe_mismatches":             np.int32,
                                                                  "probe_start_on_amplicon":      np.int32,
                                                                  "probe_strand":                 str,
                                                                  "start":                        np.int64,
                                                                  "end":                          np.int64,
                                                                  "full_hit_id":                  str,
                                                                  "scaffold_len":                 np.int64,
                                                                  "amplicon_seq":                 str,
                                                                  "primary_tag":                  str,
                                                                  "gene":                         str,
                                                                  "product":                      str,
                                                                  "protein_id":                   str,
                                                                  "note":                         str
                                                                  },
                                                   "col_name_indexes": {
                                                                        "amplicon_len":                 0,
                                                                        "scaffold_id":                  1,
                                                                        "forward_primer_id":            2,
                                                                        "forward_primer_seq":           3,
                                                                        "forward_primer_degeneracies":  4,
                                                                        "forward_primer_mismatches":    5,
                                                                        "reverse_primer_id":            6,
                                                                        "reverse_primer_seq":           7,
                                                                        "reverse_primer_degeneracies":  8,
                                                                        "reverse_primer_mismatches":    9,
                                                                        "rev_comp_reverse_primer":      10,
                                                                        "probe_id":                     11,
                                                                        "probe_seq":                    12,
                                                                        "revcomp_probe":                13,
                                                                        "probe_degeneracies":           14,
                                                                        "probe_mismatches":             15,
                                                                        "probe_start_on_amplicon":      16,
                                                                        "probe_strand":                 17,
                                                                        "start":                        18,
                                                                        "end":                          19,
                                                                        "full_hit_id":                  20,
                                                                        "scaffold_len":                 21,
                                                                        "amplicon_seq":                 22,
                                                                        "primary_tag":                  23,
                                                                        "gene":                         24,
                                                                        "product":                      25,
                                                                        "protein_id":                   26,
                                                                        "note":                         27
                                                                        },
                                                   },
                                           },
                                "filtered_str": {
                                                "coordinates_only": {
                                                                     "col_names": ["scaffold_id", "start", "end"],
                                                                     "cols":      [1, 7, 8],
                                                                     "index_cols": None,
                                                                     "converters": {
                                                                                    "scaffold_id":      str,
                                                                                    "start":            np.int64,
                                                                                    "end":              np.int64,
                                                                                   },
                                                                     "col_name_indexes": {
                                                                                          "scaffold_id":   0,
                                                                                          "start":         1,
                                                                                          "end":           2,
                                                                                         },
                                                                    },

                                                "all": {
                                                    "col_names": ["primer_pair",
                                                                  "amplicon_len",
                                                                  "scaffold_id",
                                                                  "forward_primer_id",
                                                                  "forward_primer_mismatches",
                                                                  "reverse_primer_id",
                                                                  "reverse_primer_mismatches",
                                                                  "start",
                                                                  "end",
                                                                  "amplicon_seq",
                                                                  "max_mismatches",
                                                                  "total_mismatches",
                                                                  "max_mis_min_dist",
                                                                  "tot_mis_min_dist",
                                                                  "min_len",
                                                                  "max_len",
                                                                  "len_in_expected_interval",
                                                                  "color",
                                                                  ],
                                                    "cols": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                             11, 12, 13, 14, 15, 16, 17],
                                                    "index_cols": None,
                                                    "converters": {
                                                                   "primer_pair":                 str,
                                                                   "amplicon_len":                np.int32,
                                                                   "scaffold_id":                 str,
                                                                   "forward_primer_id":           str,
                                                                   "forward_primer_mismatches":   np.int32,
                                                                   "reverse_primer_id":           str,
                                                                   "reverse_primer_mismatches":   np.int32,
                                                                   "start":                       np.int64,
                                                                   "end":                         np.int64,
                                                                   "amplicon_seq":                str,
                                                                   "max_mismatches":              np.int32,
                                                                   "total_mismatches":            np.int32,
                                                                   "max_mis_min_dist":            str,
                                                                   "tot_mis_min_dist":            str,
                                                                   "min_len":                     np.int32,
                                                                   "max_len":                     np.int32,
                                                                   "len_in_expected_interval":    np.str,
                                                                   "color":                       np.str,
                                                    },
                                                    "col_name_indexes": {
                                                                         "primer_pair":                 0,
                                                                         "amplicon_len":                1,
                                                                         "scaffold_id":                 2,
                                                                         "forward_primer_id":           3,
                                                                         "forward_primer_mismatches":   4,
                                                                         "reverse_primer_id":           5,
                                                                         "reverse_primer_mismatches":   6,
                                                                         "start":                       7,
                                                                         "end":                         8,
                                                                         "amplicon_seq":                9,
                                                                         "max_mismatches":              10,
                                                                         "total_mismatches":            11,
                                                                         "max_mis_min_dist":            12,
                                                                         "tot_mis_min_dist":            13,
                                                                         "min_len":                     14,
                                                                         "max_len":                     15,
                                                                         "len_in_expected_interval":    16,
                                                                         "color":                       17,
                                                                         },
                                                         },
                                                },
                                   }

        self.parsing_mode = parsing_mode

        self.format = format
        self.black_list = black_list
        self.white_list = white_list
        self.syn_dict = syn_dict
        self.scaffold_lengths = None

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
        """
        self.target_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_id"]
        self.target_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_start"]
        self.target_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_hit_len"]
        self.target_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_strand"]

        self.query_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_id"]
        self.query_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_start"]
        self.query_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_hit_len"]
        self.query_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_strand"]

        self.target_id_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_id"]
        self.query_id_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_id"]
        """

        self.scaffold_list = None

        if in_file:
            self.read(in_file,
                      format=format,
                      parsing_mode=parsing_mode,
                      black_list=black_list,
                      white_list=white_list
                      )

        else:
            self.records = records

    def read(self, in_file,
             format="str", parsing_mode="all",
             black_list=(), white_list=(),
             keep_seq_length_in_df=False):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        print("%s\tReading input..." % str(datetime.datetime.now()))
        self.records = pd.read_csv(in_file, sep='\t', header=0, na_values=".",
                                   comment="#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])
        self.records.index.name = "row"
        print("%s\tReading input finished..." % str(datetime.datetime.now()))
        print("%s\tFiltering..." % str(datetime.datetime.now()))
        if white_list or black_list:
            scaffolds_to_keep = self.get_filtered_entry_list(self.records["scaffold_id"].tolist(),
                                                             entry_black_list=black_list,
                                                             entry_white_list=white_list)
            #print target_scaffolds_to_keep
            self.records = self.records[self.records["scaffold_id"].isin(scaffolds_to_keep)]

        if self.syn_dict:
            self.records["scaffold_id"].replace(self.syn_dict, inplace=True)

        #                                               names=("scaffold", "row"))
        print("%s\tFiltering finished..." % str(datetime.datetime.now()))
        """
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
        """
        if "scaffold_len" in self.records.columns:
            self.scaffold_lengths = self.records[["scaffold_id", "scaffold_len"]].drop_duplicates()
            self.scaffold_lengths.columns = ["id", "length"]
            self.scaffold_lengths.set_index("id", inplace=True)
            self.scaffold_list = self.scaffold_lengths.index.tolist()
        else:
            self.scaffold_list = list(self.records["scaffold_id"].unique())

        retained_columns = deepcopy(self.parsing_parameters[self.format][self.parsing_mode]["col_names"])
        if "scaffold_len" in self.records.columns:
            if not keep_seq_length_in_df:
                for entry in ("scaffold_len",):
                    retained_columns.remove(entry)

        self.records = self.records[retained_columns]

        self.records.index.name = "row"
        self.records.set_index(keys=self.records["forward_primer_id"].apply(lambda x: x.split("|")[0]), drop=False,
                               inplace=True, append=True)
        self.records.index.set_names(("row", "primer_pair"), inplace=True)
        self.records.index = self.records.index.swaplevel(0, 1)

    def sort(self, inplace=False,
             sorting_order=None):
        if sorting_order is None:
            if self.parsing_mode == "str":
                sort_order = ["scaffold_id", "forward_primer_mismatches", "reverse_primer_mismatches", "amplicon_len"]
            elif self.parsing_mode == "filtered_str":
                sort_order = ["primer_pair", "max_mismatches", "total_mismatches", "amplicon_len"]
        else:
            sort_order = sorting_order

        return self.records.sort_values(by=sort_order, inplace=inplace)

    #def invert_coordinates(self, scaffold_id_df, scaffold_len_df):
    #    for scaffold in scaffold_id_list:
    #        self.records.loc[self.records["scaffold_id"] == scaffold, "start"] = scaff

    def filter_str(self):

        tmp_df = deepcopy(self.records)

        tmp_df["max_mismatches"] = df_dict[entry][["forward_primer_mismatches", "reverse_primer_mismatches"]].max(axis=1)
        tmp_df["total_mismatches"] = df_dict[entry][["forward_primer_mismatches", "reverse_primer_mismatches"]].sum(axis=1)
        tmp_df.sort_values(by=["primer_pair", "max_mismatches", "total_mismatches", "amplicon_len"], inplace=True)
        tmp_df[["max_mis_min_dist", "tot_mis_min_dist"]] = - tmp_df[["max_mismatches", "total_mismatches"]].diff(periods=-1)
        tmp_df.loc[tmp_df.groupby('primer_pair').tail(1).index, ["max_mis_min_dist", "tot_mis_min_dist"]] = np.nan

        return tmp_df

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

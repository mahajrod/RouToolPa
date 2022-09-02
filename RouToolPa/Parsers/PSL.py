# !/usr/bin/env python
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


class CollectionPSL:

    def __init__(self, in_file=None, records=None, format="psl", parsing_mode="all",
                 target_black_list=None, target_white_list=None,
                 query_black_list=None, query_white_list=None,
                 target_syn_dict=None, query_syn_dict=None,
                 min_target_hit_len=None, min_query_hit_len=None,
                 min_target_len=None, min_query_len=None, keep_seq_length_in_df=False,
                 invert_coordinates_for_target_negative_strand=False):

        self.formats = ["psl"]
        self.PSL_COLS = AlignmentFormats.ALN_FMT_COLS["psl"]
        self.invert_coordinates_for_target_negative_strand = invert_coordinates_for_target_negative_strand
        self.parsing_parameters = {
            "psl": {
                "coordinates_only": {
                    "col_names": ["strand",
                                  "qName",
                                  "qSize",
                                  "qStart",
                                  "qEnd",
                                  "tName",
                                  "tSize",
                                  "tStart",
                                  "tEnd"
                                  ],
                    "cols": [8, 9, 10, 11, 12, 13, 14, 15, 16, ],
                    "index_cols": None,
                    "converters": {
                        "strand": str,
                        "qName": str,
                        "qSize": np.int64,
                        "qStart": np.int64,
                        "qEnd": np.int64,
                        "tName": str,
                        "tSize": np.int64,
                        "tStart": np.int64,
                        "tEnd": np.int64,
                    },
                    "col_name_indexes": {
                        "strand": 0,
                        "qName": 1,
                        "qSize": 2,
                        "qStart": 3,
                        "qEnd": 4,
                        "tName": 5,
                        "tSize": 6,
                        "tStart": 7,
                        "tEnd": 8,
                    },
                },
                "all": {
                    "col_names": list(self.PSL_COLS.keys()),
                    "cols": list(self.PSL_COLS.values()),
                    "index_cols": None,
                    "converters": {
                        "matches": np.int64,
                        "misMatches": np.int64,
                        "repMatches": np.int64,
                        "nCount": np.int64,
                        "qNumInsert": np.int64,
                        "qBaseInsert": np.int64,
                        "tNumInsert": np.int64,
                        "tBaseInsert": np.int64,
                        "strand": str,
                        "qName": str,
                        "qSize": np.int64,
                        "qStart": np.int64,
                        "qEnd": np.int64,
                        "tName": str,
                        "tSize": np.int64,
                        "tStart": np.int64,
                        "tEnd": np.int64,
                        "blockCount": np.int64,
                        "blockSizes": str,
                        "qStarts": str,
                        "tStarts": str,
                    },
                    "col_name_indexes": self.PSL_COLS
                },
            },

        }
        self.parsing_mode = parsing_mode
        self.alignment_parsing_modes = ["all", "complete"]
        # self.attributes_parsing_modes = ["complete", "coord_and_attr"]
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
            "Int8": np.float16,
            "Int16": np.float16,
            "Int32": np.float32,
            "Int64": np.float64,
        })

        # init aliases
        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]

        self.target_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["tName"]
        self.target_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"][
            "tStart"]

        self.query_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["qName"]
        self.query_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"][
            "qStart"]
        self.query_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"][
            "strand"]

        self.target_id_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_id"]
        self.query_id_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_id"]

        self.target_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_len"]
        self.query_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_len"]

        self.target_start_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_start"]
        self.query_start_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_start"]

        self.target_end_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_end"]
        self.query_end_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_end"]

        self.target_hit_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_hit_len"]
        self.query_hit_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_hit_len"]

        self.target_strand_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_strand"]
        self.query_strand_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_strand"]

        self.target_hit_len_index = None
        self.target_strand_index = None
        self.query_hit_len_index = None

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
                      min_query_len=min_query_len,
                      keep_seq_length_in_df=keep_seq_length_in_df)

        else:
            self.records = records

    def read(self, in_file,
             format="psl", parsing_mode="only_coordinates",
             target_black_list=None, target_white_list=None,
             query_black_list=None, query_white_list=None,
             min_target_hit_len=None, min_query_hit_len=None,
             min_target_len=None, min_query_len=None,
             keep_seq_length_in_df=False):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (
            format, parsing_mode))

        print("%s\tReading input..." % str(datetime.datetime.now()))
        self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                   comment="#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])
        self.records.index.name = "row"
        print("%s\tReading input finished..." % str(datetime.datetime.now()))
        print("%s\tFiltering..." % str(datetime.datetime.now()))
        if (target_white_list is not None) or (target_black_list is not None):
            target_scaffolds_to_keep = self.get_filtered_entry_list(self.records[self.target_id_syn].tolist(),
                                                                    entry_black_list=target_black_list,
                                                                    entry_white_list=target_white_list)
            # print target_scaffolds_to_keep
            self.records = self.records[self.records[self.target_id_syn].isin(target_scaffolds_to_keep)]

        if (query_white_list is not None) or (query_black_list is not None):
            query_scaffolds_to_keep = self.get_filtered_entry_list(self.records[self.query_id_syn].tolist(),
                                                                   entry_black_list=query_black_list,
                                                                   entry_white_list=query_white_list)
            # print query_scaffolds_to_keep
            self.records = self.records[self.records[self.query_id_syn].isin(query_scaffolds_to_keep)]
        if self.target_syn_dict is not None:
            self.records[self.target_id_syn].replace(self.target_syn_dict, inplace=True)
        if self.query_syn_dict is not None:
            self.records[self.query_id_syn].replace(self.query_syn_dict, inplace=True)
        # retain only automatically generated index by row number
        # self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
        #                                               names=("scaffold", "row"))
        print("%s\tFiltering finished..." % str(datetime.datetime.now()))
        self.records["strand"].replace({"++": "+", "+-": "-"}, inplace=True)
        # convert target coordinates to forward strand, where necessary
        if self.invert_coordinates_for_target_negative_strand:
            self.records.loc[self.records["strand"] == "-", "tStart"], \
             self.records.loc[self.records["strand"] == "-", "tEnd"] = self.records.loc[self.records["strand"] == "-", "tSize"] - self.records.loc[self.records["strand"] == "-", "tEnd"], \
                                                                       self.records.loc[self.records["strand"] == "-", "tSize"] - self.records.loc[self.records["strand"] == "-", "tStart"]

        if min_target_len and min_query_len:
            self.records = self.records[
                (self.records[self.target_len_syn] >= min_target_len) & (self.records[self.query_len_syn] >= min_query_len)]
        elif min_target_len:
            self.records = self.records[self.records[self.target_len_syn] >= min_target_len]
        elif min_query_len:
            self.records = self.records[self.records[self.query_len_syn] >= min_query_len]

        self.records[self.target_hit_len_syn] = self.records[self.target_end_syn] - self.records[self.target_start_syn]
        self.records[self.query_hit_len_syn] = self.records[self.query_end_syn] - self.records[self.query_start_syn]
        # print(self.records.columns)
        if min_target_hit_len and min_query_hit_len:
            self.records = self.records[(self.target_hit_len_syn >= min_target_hit_len) & (
                        self.query_hit_len_syn >= min_query_hit_len)]
        elif min_target_hit_len:
            self.records = self.records[self.target_hit_len_syn >= min_target_hit_len]
        elif min_query_hit_len:
            self.records = self.records[self.query_hit_len_syn >= min_query_hit_len]

        self.target_scaffold_lengths = self.records[[self.target_id_syn, self.target_len_syn]].drop_duplicates()
        self.target_scaffold_lengths.columns = ["id", "length"]
        self.target_scaffold_lengths.set_index("id", inplace=True)
        self.query_scaffold_lengths = self.records[[self.query_id_syn,
                                                    self.query_len_syn]].drop_duplicates()
        self.query_scaffold_lengths.columns = ["id", "length"]
        self.query_scaffold_lengths.set_index("id", inplace=True)
        self.target_scaffold_list = self.target_scaffold_lengths.index.tolist()
        self.query_scaffold_list = self.query_scaffold_lengths.index.tolist()

        retained_columns = deepcopy(self.parsing_parameters[self.format][self.parsing_mode]["col_names"] + [self.target_hit_len_syn, self.query_hit_len_syn])
        if not keep_seq_length_in_df:
            for entry in self.target_len_syn, self.query_len_syn:
                retained_columns.remove(entry)

        self.records = self.records[retained_columns]

    def sort(self, inplace=False,
             sorting_order=(AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_id"],
                            AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_id"],
                            AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_start"],
                            AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["target_hit_len"],
                            AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_start"],
                            AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["psl"]["query_hit_len"],
                            )):
        return self.records.sort_values(by=sorting_order, inplace=inplace)

    @staticmethod
    def get_filtered_entry_list(entry_list,
                                entry_black_list=None,
                                sort_entries=False,
                                entry_ordered_list=None,
                                entry_white_list=None):
        white_set = set(entry_white_list) if entry_white_list is not None else set()
        black_set = set(entry_black_list) if entry_black_list is not None else set()
        entry_set = set(entry_list)

        if white_set:
            entry_set = entry_set & white_set
        if black_set:
            entry_set = entry_set - black_set

        filtered_entry_list = list(entry_set)
        if sort_entries:
            filtered_entry_list.sort()

        final_entry_list = []

        if entry_ordered_list is not None:
            for entry in entry_ordered_list:
                if entry in filtered_entry_list:
                    final_entry_list.append(entry)
                    filtered_entry_list.remove(entry)
                else:
                    print("WARNING!!!Entry(%s) from order list is absent in list of entries!" % entry)
            return final_entry_list + filtered_entry_list
        else:
            return filtered_entry_list

    def rename_target_ids(self, syn_dict):
        self.records[self.target_id_syn].replace(syn_dict, inplace=True)

    def rename_query_ids(self, syn_dict):
        self.records[self.query_id_syn].replace(syn_dict, inplace=True)


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

import RouToolPa.Formats.ScaffoldingFormats as ScaffoldingFormats

AGP_PART_TYPES_INT_COORDINATES = ["W"]


class CollectionAGP:

    def __init__(self, in_file=None, records=None, format="agp", parsing_mode="all",
                 black_list=(), white_list=()):

        self.formats = ["agp"]
        self.TAB_COLS = ScaffoldingFormats.SCAF_FMT_COLS["agp"]
        self.parsing_parameters = {
                                "agp": {
                                           "all": {
                                                   "col_names": ["scaffold", "start", "end", "part_number", "part_type", "part_id/gap_length", "part_start/gap_type", "part_end/linkage", "orientation/evidence"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8],
                                                   "index_cols": None,
                                                   "converters": {
                                                                   "scaffold": str,
                                                                   "start": np.int64,
                                                                   "end": np.int64,
                                                                   "part_number": np.int32,
                                                                   "part_type": str,
                                                                   "part_id/gap_length": str,
                                                                   "part_start/gap_type": str,
                                                                   "part_end/linkage": str,
                                                                   "orientation/evidence": str
                                                                  },
                                                   "col_name_indexes": ScaffoldingFormats.SCAF_FMT_COLS["agp"]
                                                   },
                                    },

                                }
        self.format = format
        self.parsing_mode = parsing_mode

        self.pandas_int_type_correspondence = OrderedDict({
                                                           "Int8":  np.float16,
                                                           "Int16": np.float16,
                                                           "Int32": np.float32,
                                                           "Int64": np.float64,
                                                           })

        self.target_scaffold_list = None
        self.query_scaffold_list = None
        self.target_scaffold_lengths = None
        self.query_scaffold_lengths = None
        self.length_df = None

        self.target_black_list = black_list
        self.target_white_list = white_list

        if in_file:
            self.read(in_file,
                      format=format,
                      parsing_mode=parsing_mode,
                      black_list=black_list,
                      white_list=white_list)
                      #target_black_list=target_black_list,
                      #target_white_list=target_white_list,
                      #query_black_list=query_black_list,
                      #query_white_list=query_white_list,
                      #min_target_hit_len=min_target_hit_len,
                      #min_query_hit_len=min_query_hit_len,
                      #min_target_len=min_target_len,
                      #min_query_len=min_query_len,
                      #keep_seq_length_in_df=keep_seq_length_in_df)

        else:
            self.records = records

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

    def read(self, in_file, format="agp", parsing_mode="all",
             black_list=None, white_list=None):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % format)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        print("%s\tReading input..." % str(datetime.datetime.now()))
        self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                   comment="#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])

        if (white_list is not None) or (black_list is not None):
            scaffolds_to_keep = self.get_filtered_entry_list(self.records["scaffold"].tolist(),
                                                             entry_black_list=black_list,
                                                             entry_white_list=white_list)
            self.records = self.records[self.records["scaffold"].isin(scaffolds_to_keep)]

        if self.parsing_parameters[format][parsing_mode]["index_cols"] is None:
            self.records["segment_id"] = [f"SEG{i}" for i in range(1, len(self.records) + 1)]
            self.records.set_index("segment_id", inplace=True)

        #self.records.index.name = "row"
        self.records["start"] = self.records["start"] - 1
        self.records.loc[self.records["part_type"].isin(AGP_PART_TYPES_INT_COORDINATES), "part_start/gap_type"] = self.records.loc[self.records["part_type"].isin(AGP_PART_TYPES_INT_COORDINATES), "part_start/gap_type"].astype(int) - 1
        self.records.loc[self.records["part_type"].isin(AGP_PART_TYPES_INT_COORDINATES), "part_end/linkage"] = self.records.loc[self.records["part_type"].isin(AGP_PART_TYPES_INT_COORDINATES), "part_end/linkage"].astype(int)
        self.get_seq_len()
        print("%s\tReading input finished..." % str(datetime.datetime.now()))
        print("%s\tFiltering..." % str(datetime.datetime.now()))
        """
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

        self.target_scaffold_lengths = self.records[["target_id", "target_len"]].drop_duplicates()
        self.target_scaffold_lengths.columns = ["id", "length"]
        self.target_scaffold_lengths.set_index("id", inplace=True)
        self.query_scaffold_lengths = self.records[["query_id", "query_len"]].drop_duplicates()
        self.query_scaffold_lengths.columns = ["id", "length"]
        self.query_scaffold_lengths.set_index("id", inplace=True)
        self.target_scaffold_list = self.target_scaffold_lengths.index.tolist()
        self.query_scaffold_list = self.query_scaffold_lengths.index.tolist()

        retained_columns = deepcopy(self.parsing_parameters[self.format][self.parsing_mode]["col_names"])
        if not keep_seq_length_in_df:
            for entry in "target_len", "query_len":
                retained_columns.remove(entry)

        self.records = self.records[retained_columns]

        if (self.format == "tab") and (parsing_mode == "complete"):
            self.records["EG2"] = map(lambda s: np.float32(s.split("=")[1]), list(self.records["EG2"]))
            self.records["E"] = map(lambda s: np.float32(s.split("=")[1]), list(self.records["E"]))
        elif (self.format == "tab_mismap") and (parsing_mode == "complete"):
            self.records["mismap"] = map(lambda s: np.float32(s.split("=")[1]), list(self.records["mismap"]))
        """

    def filter_records(self, white_list=None, black_list=None):
        if (white_list is not None) or (black_list is not None):
            scaffolds_to_keep = self.get_filtered_entry_list(self.records["scaffold"].tolist(),
                                                             entry_black_list=black_list,
                                                             entry_white_list=white_list)
            self.records = self.records[self.records["scaffold"].isin(scaffolds_to_keep)]

    def get_seq_len(self):
        self.length_df = self.records.groupby(by="scaffold").apply(lambda df: df[["end"]].iloc[-1])
        self.length_df.columns = pd.Index(["length"])
        self.length_df.sort_values(by="length", ascending=False, inplace=True)
        return self.length_df

    def get_parts_df(self, part_type_whitelist=("W",), part_type_blacklist=(), keep_index=True, sort=False):
        if part_type_whitelist and part_type_whitelist:
            filtered_df = self.records[self.records["part_type"].isin(part_type_whitelist) & (~self.records["part_type"].isin(part_type_blacklist))]
        elif part_type_whitelist:
            filtered_df = self.records[self.records["part_type"].isin(part_type_whitelist)]
        elif part_type_blacklist:
            filtered_df = self.records[~self.records["part_type"].isin(part_type_blacklist)]
        else:
            filtered_df = self.records

        agp_parts_df = filtered_df[["part_id/gap_length",
                                    "part_start/gap_type",
                                    "part_end/linkage"]]
        agp_parts_df.loc[:, "part_start/gap_type"] = agp_parts_df["part_start/gap_type"].astype(int)
        agp_parts_df.loc[:, "part_end/linkage"] = agp_parts_df["part_end/linkage"].astype(int)

        if sort:
            agp_parts_df = agp_parts_df.sort_values(by=["part_id/gap_length", "part_start/gap_type", "part_end/linkage"])
        if not keep_index:
            agp_parts_df.set_index("part_id/gap_length", inplace=True)

        return agp_parts_df

    def get_nongap_records(self): # coordinates of segments are unchanged!!!!
        return self.records[self.records["part_type"] != "U"]

    def write(self, output_agp):
        # convert internal 0-based coordinate system to 1-based used by agp
        tmp = deepcopy(self.records)
        tmp.loc[:, "start"] += 1
        tmp.loc[self.records["part_type"].isin(AGP_PART_TYPES_INT_COORDINATES), "part_start/gap_type"] += 1
        tmp.to_csv(output_agp, sep="\t", index=False, header=False)

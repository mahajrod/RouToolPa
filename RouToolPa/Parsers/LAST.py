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


class CollectionLast:

    def __init__(self, in_file=None, records=None, format="tab", parsing_mode="all",
                 target_black_list=(), target_white_list=(),
                 query_black_list=(), query_white_list=(),
                 target_syn_dict=None, query_syn_dict=None,
                 min_target_hit_len=None, min_query_hit_len=None,
                 min_target_len=None, min_query_len=None, keep_seq_length_in_df=False):

        self.formats = ["tab", "tab_mismap"]
        self.TAB_COLS = AlignmentFormats.ALN_FMT_COLS["tab"]
        self.parsing_parameters = {
                                "tab": {
                                           "coordinates_only": {
                                                   "col_names": ["target_id", "target_start", "target_hit_len",
                                                                 "target_strand", "target_len", "query_id",
                                                                 "query_start", "query_hit_len", "query_strand",
                                                                 "query_len"],
                                                   "cols":      [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                   "index_cols": None,
                                                   "converters": {
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
                                                                  },
                                                   "col_name_indexes": {
                                                                        "target_id":   0,
                                                                        "target_start": 1,
                                                                        "target_hit_len": 2,
                                                                        "target_strand": 3,
                                                                        "target_len": 4,
                                                                        "query_id": 5,
                                                                        "query_start": 6,
                                                                        "query_hit_len": 7,
                                                                        "query_strand": 8,
                                                                        "query_len": 9,
                                                                        },
                                                                 },

                                           "all": {
                                                   "col_names": ["score", "target_id", "target_start", "target_hit_len",
                                                                 "target_strand", "target_len", "query_id",
                                                                 "query_start", "query_hit_len", "query_strand",
                                                                 "query_len", "alignment", "EG2", "E"],
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
                                                                  "EG2":             str,
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
                                                                        "EG2": 12,
                                                                        "E": 13,
                                                                        },
                                                   },
                                           "complete": {
                                                        "col_names": ["score", "target_id", "target_start", "target_hit_len",
                                                                      "target_strand", "target_len", "query_id",
                                                                      "query_start", "query_hit_len", "query_strand",
                                                                      "query_len", "alignment", "EG2", "E"],
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
                                                                       "EG2":             str,
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
                                                                             "EG2": 12,
                                                                             "E": 13,
                                                                             },
                                                   },
                                           },

                                "tab_mismap": {
                                                "coordinates_only": {
                                                    "col_names": ["target_id", "target_start", "target_hit_len",
                                                                  "target_strand", "target_len", "query_id",
                                                                  "query_start", "query_hit_len", "query_strand",
                                                                  "query_len"],
                                                    "cols": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                    "index_cols": None,
                                                    "converters": {
                                                        "target_id": str,
                                                        "target_start": np.int64,
                                                        "target_hit_len": np.int64,
                                                        "target_strand": str,
                                                        "target_len": np.int64,
                                                        "query_id": str,
                                                        "query_start": np.int64,
                                                        "query_hit_len": np.int64,
                                                        "query_strand": str,
                                                        "query_len": np.int64,
                                                    },
                                                    "col_name_indexes": {
                                                        "target_id": 0,
                                                        "target_start": 1,
                                                        "target_hit_len": 2,
                                                        "target_strand": 3,
                                                        "target_len": 4,
                                                        "query_id": 5,
                                                        "query_start": 6,
                                                        "query_hit_len": 7,
                                                        "query_strand": 8,
                                                        "query_len": 9,
                                                    },
                                                },

                                                "all": {
                                                    "col_names": ["score", "target_id", "target_start", "target_hit_len",
                                                                  "target_strand", "target_len", "query_id",
                                                                  "query_start", "query_hit_len", "query_strand",
                                                                  "query_len", "alignment", "mismap"],
                                                    "cols": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                                    "index_cols": None,
                                                    "converters": {
                                                        "score": np.int64,
                                                        "target_id": str,
                                                        "target_start": np.int64,
                                                        "target_hit_len": np.int64,
                                                        "target_strand": str,
                                                        "target_len": np.int64,
                                                        "query_id": str,
                                                        "query_start": np.int64,
                                                        "query_hit_len": np.int64,
                                                        "query_strand": str,
                                                        "query_len": np.int64,
                                                        "alignment": str,
                                                        "mismap": str,
                                                    },
                                                    "col_name_indexes": {
                                                        "score": 0,
                                                        "target_id": 1,
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
                                                        "mismap": 12,
                                                    },
                                                },
                                                "complete": {
                                                    "col_names": ["score", "target_id", "target_start", "target_hit_len",
                                                                  "target_strand", "target_len", "query_id",
                                                                  "query_start", "query_hit_len", "query_strand",
                                                                  "query_len", "alignment", "mismap"],
                                                    "cols": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                                    "index_cols": None,
                                                    "converters": {
                                                        "score": np.int64,
                                                        "target_id": str,
                                                        "target_start": np.int64,
                                                        "target_hit_len": np.int64,
                                                        "target_strand": str,
                                                        "target_len": np.int64,
                                                        "query_id": str,
                                                        "query_start": np.int64,
                                                        "query_hit_len": np.int64,
                                                        "query_strand": str,
                                                        "query_len": np.int64,
                                                        "alignment": str,
                                                        "mismap": str
                                                    },
                                                    "col_name_indexes": {
                                                        "score": 0,
                                                        "target_id": 1,
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
                                                        "mismap": 12,
                                                    },
                                                },
                                },
                                   }
        self.parsing_mode = parsing_mode
        self.alignment_parsing_modes = ["all", "complete"]
        #self.attributes_parsing_modes = ["complete", "coord_and_attr"]
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
        self.target_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_hit_len"]
        self.target_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["target_strand"]

        self.query_id_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_id"]
        self.query_start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_start"]
        self.query_hit_len_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_hit_len"]
        self.query_strand_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["query_strand"]

        self.target_id_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_id"]
        self.query_id_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_id"]

        self.target_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_len"]
        self.query_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_len"]

        self.target_start_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_start"]
        self.query_start_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_start"]

        self.target_end_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_end"]
        self.query_end_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_end"]

        self.target_hit_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_hit_len"]
        self.query_hit_len_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_hit_len"]

        self.target_strand_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["target_strand"]
        self.query_strand_syn = AlignmentFormats.ALN_FMT_COLUMN_NAMES_SYN["tab"]["query_strand"]

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
                      min_query_len=min_query_len,
                      keep_seq_length_in_df=keep_seq_length_in_df)

        else:
            self.records = records

    def read(self, in_file,
             format="tab", parsing_mode="only_coordinates",
             target_black_list=(), target_white_list=(),
             query_black_list=(), query_white_list=(),
             min_target_hit_len=None, min_query_hit_len=None,
             min_target_len=None, min_query_len=None,
             keep_seq_length_in_df=False):
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
        self.records.index.name = "row"
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

    def sort(self, inplace=False,
             sorting_order=("target_id", "query_id", "target_start", "target_hit_len", "query_start", "query_hit_len")):
        return self.records.sort_values(by=sorting_order, inplace=inplace)

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

    def write(self, output, separator="\t", header=False):
        self.records.to_csv(output,
                            sep=separator,
                            header=header,
                            index=False)

    def merge_adjucent_blocks(self, max_dist_to_collapse=None):
        """
        ONLY COORDINATES ARE MERGED. OTHER FIELDS HAVE NO SENSE AFTER THIS PROCEDUREAND ARE RETAINED ONLY FOR PARSING COMPATIBILITY!
        :return:
        """
        row_iterator = self.records.itertuples(index=False)
        merged_row_list = []

        for row in row_iterator:
            curr_row = list(row)
            for index in self.target_id_index, self.target_strand_index, self.query_id_index, self.query_strand_index:

                if prev_row[index] != curr_row[index]:
                    merged_row_list.append(prev_row)
                    prev_row = curr_row
                    break
            else:
                if max_dist_to_collapse:
                    if (curr_row[self.target_start_index] - prev_row[self.target_start_index] -
                                prev_row[self.target_hit_len_index] + 1 <= max_dist_to_collapse) and \
                            (curr_row[self.query_start_index] - prev_row[self.query_start_index] -
                                 prev_row[self.query_hit_len_index] <= max_dist_to_collapse):
                        prev_row[self.target_hit_len_index] += curr_row[self.target_hit_len_index]
                        prev_row[self.query_hit_len_index] += curr_row[self.query_hit_len_index]
                    else:
                        merged_row_list.append(prev_row)
                        prev_row = curr_row
                else:
                    prev_row[self.target_hit_len_index] += curr_row[self.target_hit_len_index]
                    prev_row[self.query_hit_len_index] += curr_row[self.query_hit_len_index]

        merged_row_list.append(prev_row)

    def extract_coordinates(self):

        return self.records[["target_id", "target_start", "target_hit_len", "target_strand",
                             "query_id", "query_start", "query_hit_len", "query_strand"]]

    def extract_coordinates_to_file(self, output, separator="\t", header=False):
        self.records.to_csv(output,
                            columns=["target_id", "target_start", "target_hit_len", "target_strand",
                                     "query_id", "query_start", "query_hit_len", "query_strand"],
                            sep=separator,
                            header=header,
                            index=False)

    def extract_and_convert_coordinates(self, output=None, separator="\t", header=False):
        coords = deepcopy(self.records[["target_id", "target_start", "target_hit_len", "target_strand",
                                        "query_id", "query_start", "query_hit_len", "query_strand"]])
        coords.columns = ["target_id", "target_start", "target_end", "target_strand",
                        "query_id", "query_start", "query_end", "query_strand"]
        coords["target_end"] = coords["target_start"] + coords["target_end"]
        coords["query_end"] = coords["query_start"] + coords["query_end"]
        if output:
            coords.to_csv(output, sep=separator, header=header,
                          index=False)

        return coords

    def collapse_coordinates(self, output=None, max_dist_to_collapse=1000, separator="\t", header=False):

        coords = self.extract_and_convert_coordinates()
        row_iterator = coords.itertuples(index=False)

        prev_row = list(row_iterator.next())
        collapsed_row_list = []

        for row in row_iterator:
            curr_row = list(row)
            for index in 0, 3, 4, 7:
                # check for target and query scaffold and id
                if prev_row[index] != curr_row[index]:
                    collapsed_row_list.append(prev_row)
                    prev_row = curr_row
                    break
            else:
                if max_dist_to_collapse:
                    if (curr_row[1] - prev_row[2] <= max_dist_to_collapse) and (curr_row[5] - prev_row[6] <= max_dist_to_collapse):
                        prev_row[2] = curr_row[2]
                        prev_row[6] = curr_row[6]
                    else:
                        collapsed_row_list.append(prev_row)
                        prev_row = curr_row
                else:
                    prev_row[2] = curr_row[2]
                    prev_row[6] = curr_row[6]

        collapsed_row_list.append(prev_row)

        collapsed_row_list = pd.DataFrame.from_records(collapsed_row_list, columns=coords.columns)

        if output:
            collapsed_row_list.to_csv(output, sep=separator, header=header, index=False)

        return collapsed_row_list

    def rename_target_ids(self, syn_dict):
        self.records[["target_id"]].replace(syn_dict, inplace=True)

    def rename_query_ids(self, syn_dict):
        self.records[["query_id"]].replace(syn_dict, inplace=True)

    @staticmethod
    def parse_alignment_string(alignment_string):

        str_list = map(lambda s: s.split(":"), alignment_string.split(","))

        target_list = []
        query_list = []

        for entry in str_list:
            target_list.append(int(entry[0]))
            query_list.append(int(entry[0]) if len(entry) == 1 else int(entry[1]))
        
        return np.array(target_list), np.array(query_list)

    def transfer_coordinates(self, scaffold, start, stop, source="target"):

        if source == "target":
            pass
        elif source == "query":
            pass
        else:
            raise ValueError("ERROR!!! Unrecognized source! Only 'target' or 'query' are allowed.")

        self.records.loc[scaffold]
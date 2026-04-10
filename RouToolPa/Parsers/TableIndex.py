# !/usr/bin/env python
"""
Last output parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'
import sys
import datetime
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import pandas as pd

import RouToolPa.Formats.AnnotationFormats as AnnotationFormats


class CollectionTableIndex:

    def __init__(self, in_file=None, records=None,  records_columns=None, format="fai", parsing_mode="complete",
                 scaffold_black_list=None, scaffold_white_list=None, scaffold_syn_dict=None, header_in_file=False,
                 rename_dict=None, white_list_dict={}, black_list_dict={},
                 scaffold_column_name=None, length_column_name=None,
                 ):

        self.formats = ["fai", "fastq_fai"]
        self.parsing_mode = parsing_mode
        self.parsing_parameters = {
            "fai": {
                "length_only": {
                    "col_names": ["scaffold",
                                  "length",
                                  ],
                    "cols": [0, 1, ],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "length": np.int64,
                        },
                    "col_name_indexes": {
                                        "scaffold": 0,
                                        "length": 1,
                                        },
                    "original_col_names": {},
                    },
                "complete": {
                    "col_names": ["scaffold",
                                  "length",
                                  "offset",
                                  "linebases",
                                  "linewidth"
                                  ],
                    "cols": [0, 1, 2, 3, 4],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "length": np.int64,
                        "offset": np.int64,
                        "linebases": np.float64,
                        "linewidth": np.float64,
                        },
                    "col_name_indexes": {
                        "scaffold":   0,
                        "length":     1,
                        "offset":     2,
                        "linebases":  3,
                        "linewidth":  4
                        },
                    "original_col_names": {},
                    },

                },
            "fastq_fai": {
                "length_only": {
                    "col_names": ["scaffold",
                                  "length",
                                  ],
                    "cols": [0, 1, ],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "length": np.int64,
                    },
                    "col_name_indexes": {
                        "scaffold": 0,
                        "length": 1,
                    },
                    "original_col_names": {},
                },
                "complete": {
                    "col_names": ["scaffold",
                                  "length",
                                  "offset",
                                  "linebases",
                                  "linewidth",
                                  "qualoffset"
                                  ],
                    "cols": [0, 1, 2, 3, 4],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "length": np.int64,
                        "offset": np.int64,
                        "linebases": np.float64,
                        "linewidth": np.float64,
                        "qualoffset": np.float64
                    },
                    "col_name_indexes": {
                        "scaffold":    0,
                        "length":      1,
                        "offset":      2,
                        "linebases":   3,
                        "linewidth":   4,
                        "qualoffset":  5
                    },
                    "original_col_names": {},
                },

            },
                                  }
        self.parsing_mode = parsing_mode
        self.format = format
        self.scaffold_black_list = scaffold_black_list
        self.scaffold_white_list = scaffold_white_list
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

        self.scaffold_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["scaffold"]
        self.length_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["length"]

        self.scaffold_syn = "scaffold"
        self.length_syn = "length"

        self.scaffold_syn_dict = scaffold_syn_dict
        self.scaffold_list = []
        if in_file:
            self.read(in_file,
                      format=format,
                      parsing_mode=parsing_mode,
                      scaffold_black_list=scaffold_black_list,
                      scaffold_white_list=scaffold_white_list,
                      scaffold_syn_dict=scaffold_syn_dict,
                      header_in_file=header_in_file,
                      records_columns=records_columns,
                      rename_dict=rename_dict,
                      white_list_dict=white_list_dict,
                      black_list_dict=black_list_dict,
                      scaffold_column_name=scaffold_column_name,
                      length_column_name=length_column_name,
                      )

        else:
            self.records = records

            if records_columns:
                self.records.columns = pd.Index(records_columns)

    def read(self, in_file,
             format="bed", parsing_mode="only_coordinates",
             scaffold_black_list=None, scaffold_white_list=None,
             scaffold_syn_dict=None,
             header_in_file=False,
             records_columns=None,
             rename_dict=None,
             white_list_dict={},
             black_list_dict={},
             scaffold_column_name=None,
             length_column_name=None):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        sys.stderr.write(f"{str(datetime.datetime.now())}\tReading input: format={format}, parsing mode={parsing_mode}...\n")

        self.records = pd.read_csv(in_file, sep='\t', header=None if header_in_file is False else 0, na_values=".",
                                   comment="#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])

        if records_columns is not None:
            self.records.columns = records_columns if isinstance(pd.Index, records_columns) else pd.Index(records_columns)
        sys.stderr.write("%s\tReading input finished...\n" % str(datetime.datetime.now()))
        sys.stderr.write("%s\tFiltering...\n" % str(datetime.datetime.now()))
        if (scaffold_white_list is not None) or (scaffold_black_list is not None):
            scaffolds_to_keep = self.get_filtered_entry_list(self.records[self.scaffold_syn].tolist(),
                                                             entry_black_list=scaffold_black_list,
                                                             entry_white_list=scaffold_white_list)
            self.records = self.records[self.records[self.scaffold_syn].isin(scaffolds_to_keep)]

        columns_to_filter_set = set(white_list_dict.keys()) | set(black_list_dict.keys())
        for column in columns_to_filter_set:
            black_list = black_list_dict[column] if column in black_list_dict else None

            values_to_keep = self.get_filtered_entry_list(self.records[column].tolist(),
                                                          entry_black_list=black_list_dict[column] if column in black_list_dict else None,
                                                          entry_white_list=white_list_dict[column] if column in white_list_dict else None)

            self.records = self.records[self.records[column].isin(values_to_keep)]

        if self.scaffold_syn_dict is not None:
            self.records.rename(index=scaffold_syn_dict, inplace=True)

        if rename_dict is not None:
            for column in rename_dict:
                self.records[column].replace(rename_dict[column], inplace=True)

        sys.stderr.write("%s\tFiltering finished...\n" % str(datetime.datetime.now()))

        self.scaffold_list = self.records.index.tolist()

    def sort(self, sorting_order, inplace=False):
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
                    sys.stderr.write("WARNING!!!Entry(%s) from order list is absent in list of entries!\n" % entry)
            return final_entry_list + filtered_entry_list
        else:
            return filtered_entry_list

    def rename_scaffolds(self, syn_dict):
        self.records[self.scaffold_syn].replace(syn_dict, inplace=True)

    def rename_column_values(self, column_id, syn_dict):
        self.records[column_id].replace(syn_dict, inplace=True)

    def write(self, out_file, write_header=False):
        self.records.to_csv(out_file, sep="\t", index=True, header=write_header,)

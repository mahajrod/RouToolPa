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

import RouToolPa.Formats.AnnotationFormats as AnnotationFormats


class CollectionBED:

    def __init__(self, in_file=None, records=None,  records_columns=None, format="bed", parsing_mode="all",
                 scaffold_black_list=None, scaffold_white_list=None, scaffold_syn_dict=None, header_in_file=False,
                 rename_dict=None, white_list_dict={}, black_list_dict={}
                 ):

        self.formats = ["bed"]
        self.parsing_parameters = {
            "bed": {
                "coordinates_only": {
                    "col_names": ["scaffold",
                                  "start",
                                  "end",
                                  ],
                    "cols": [0, 1, 2, ],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "start": np.int64,
                        "end": np.int64,
                    },
                    "col_name_indexes": {
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                    },
                },
                "all": {
                    "col_names": None,
                    "cols": None,
                    "index_cols": 0,
                    "converters": {},
                    "col_name_indexes": {
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                    },
                },
            },
            "bed_track": {
                "coordinates_only": {
                    "col_names": ["scaffold",
                                  "start",
                                  "end",
                                  ],
                    "cols": [0, 1, 2, ],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "start": np.int64,
                        "end": np.int64,
                    },
                    "col_name_indexes": {
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                    },
                },
                "all": {
                    "col_names": ["scaffold",
                                  "start",
                                  "end",
                                  "query",
                                  "color"
                                  ],
                    "cols": [0, 1, 2, 3, 4],
                    "index_cols": 0,
                    "converters": {
                        "scaffold": str,
                        "start": np.int64,
                        "end": np.int64,
                        "query": str,
                        "color": str
                    },
                    "col_name_indexes": {
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                        "query": 3,
                        "color": 4
                    },
                },
                            },
            "bed_synteny_track": {
                "coordinates_only": {
                    "col_names": ["scaffold",
                                  "start",
                                  "end",
                                  ],
                    "cols": [0, 1, 2, ],
                    "index_cols": "scaffold",
                    "converters": {
                        "scaffold": str,
                        "start": np.int64,
                        "end": np.int64,
                    },
                    "col_name_indexes": {
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                    },
                },
                "all": {
                    "col_names": ["scaffold",
                                  "start",
                                  "end",
                                  "query",
                                  "query_start",
                                  "query_end",
                                  "strand",
                                  "color"
                                  ],
                    "cols": [0, 1, 2, 3, 4, 5, 6, 7],
                    "index_cols": 0,
                    "converters": {
                        "scaffold": str,
                        "start": np.int64,
                        "end": np.int64,
                        "query": str,
                        "query_start": np.int64,
                        "query_end": np.int64,
                        "strand": str,
                        "color": str

                    },
                    "col_name_indexes": {
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                        "query": 3,
                        "query_start": 4,
                        "query_end": 5,
                        "strand": 6,
                        "color": 7
                    },
                },
            },
                                  }
        self.parsing_mode = parsing_mode
        self.alignment_parsing_modes = ["all", "complete"]
        # self.attributes_parsing_modes = ["complete", "coord_and_attr"]
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
        self.start_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["start"]
        self.end_index = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["end"]

        self.scaffold_syn = "scaffold"
        self.start_syn = "start"
        self.end_syn = "end"

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
                      black_list_dict=black_list_dict
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
             black_list_dict={}):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        print("%s\tReading input..." % str(datetime.datetime.now()))
        self.records = pd.read_csv(in_file, sep='\t', header=None if header_in_file is None else 0, na_values=".",
                                   comment=None if (self.format =="bed_track") or (self.format =="bed_synteny_track") else "#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])
        if records_columns is not None:
            self.records.columns = records_columns if isinstance(pd.Index, records_columns) else pd.Index(records_columns)
        print("%s\tReading input finished..." % str(datetime.datetime.now()))
        print("%s\tFiltering..." % str(datetime.datetime.now()))
        if (scaffold_white_list is not None) or (scaffold_black_list is not None):
            scaffolds_to_keep = self.get_filtered_entry_list(self.records[self.scaffold_syn].tolist(),
                                                             entry_black_list=scaffold_black_list,
                                                             entry_white_list=scaffold_white_list)
            # print target_scaffolds_to_keep
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

        print("%s\tFiltering finished..." % str(datetime.datetime.now()))

        self.scaffold_list = self.records.index.tolist()

    def sort(self, inplace=False,
             sorting_order=(AnnotationFormats.BED_COLS["scaffold"],
                            AnnotationFormats.BED_COLS["start"],
                            AnnotationFormats.BED_COLS["end"],
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

    def rename_scaffolds(self, syn_dict):
        self.records[self.scaffold_syn].replace(syn_dict, inplace=True)

    def rename_column_values(self, column_id, syn_dict):
        self.records[column_id].replace(syn_dict, inplace=True)


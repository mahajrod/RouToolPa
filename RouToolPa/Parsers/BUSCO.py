#!/usr/bin/env python
"""
GFF Parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'
import datetime
from copy import deepcopy
from collections import OrderedDict, Iterable
import numpy as np
import pandas as pd

from RouToolPa.Parsers.Abstract import Parser
import RouToolPa.Formats.AnnotationFormats as AnnotationFormats


class BUSCOtable(Parser):

    def __init__(self, in_file=None, records=None, busco_version=None, dataset=None, dataset_size=None,
                 format="busco_table", parsing_mode="complete",
                 black_list=(), white_list=(), add_row_index=True):

        self.formats = ["busco_table"]

        self.status_list = ["Complete",
                            "Duplicated",
                            "Fragmented",
                            "Missing"]

        self.parsing_parameters = {"busco_table": {

                                           "complete": {
                                                   "col_names": ["id", "status", "scaffold", "start", "end", "strand",
                                                                 "score", "length", "url", "description"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                                                   "index_cols": ["status"],
                                                   "converters": {
                                                                  "id":             str,
                                                                  "status":         str,
                                                                  "scaffold":       str,
                                                                  "start":          str, #lambda x: np.int32(x) - 1,
                                                                  "end":            str, #np.int32,
                                                                  "strand":         str,
                                                                  "score":          str, #float,
                                                                  "length":         str, #np.int32,
                                                                  "url":            str, #np.int32,
                                                                  "description":    str, #np.int32,
                                                                  },
                                                   "col_name_indexes": {
                                                                        "id":           0,
                                                                        "status":       1,
                                                                        "scaffold":     2,
                                                                        "start":        3,
                                                                        "end":          4,
                                                                        "strand":       5,
                                                                        "score":        6,
                                                                        "length":       7,
                                                                        "url":          8,
                                                                        "description":  9,
                                                                        },
                                                   },
                                           },
                                    }

        self.parsing_mode = parsing_mode
        self.format = format
        self.busco_version = None
        self.dataset = None
        self.dataset_size = None
        self.original_header = None

        self.black_list = black_list
        self.white_list = white_list

        # init aliases
        self.record_id_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["id"]
        self.record_start_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["start"]
        self.record_scaffold_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["scaffold"]
        self.record_end_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["end"]
        self.record_status_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["status"]

        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]

        # load records
        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode, black_list=black_list, white_list=white_list,
                      add_row_index=add_row_index)
            self.read_metadata(in_file)

        else:
            self.records = records
            self.busco_version = None
            self.dataset = None
            self.dataset_size = None

    def read_metadata(self, in_file, format="busco_table"):
        version_parsed = False
        dataset_parsed = False
        with open(in_file, "r") as in_fd:
            if format == "busco_table":
                for line in in_fd:
                    if line[0] == "# Busco id":
                        self.original_header = line
                        break
                    if "BUSCO version is:" in line:
                        if version_parsed:
                            raise ValueError("ERROR!!! Two or more lines of metadata contains BUSCO version!")
                        self.busco_version = line.strip().split()[-1]
                        version_parsed = True
                    if "The lineage dataset is:" in line:
                        if dataset_parsed:
                            raise ValueError("ERROR!!! Two or more lines of metadata contains datase info!")
                        line_list = line.split("The lineage dataset is:")[-1].strip().split()
                        self.dataset = line_list[0]
                        self.dataset_size = int(line_list[-1][:-1])
                        dataset_parsed = True

    def read(self, in_file, format="busco_table", parsing_mode="only_coordinates", sort=False,
             black_list=(), white_list=(), add_row_index=True):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        print("%s\tReading input..." % str(datetime.datetime.now()))
        if format == "busco_table":
            self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                       comment="#",
                                       usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                       converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                       names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                       index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])

            if white_list or black_list:
                scaffolds_to_keep = self.get_filtered_entry_list(self.records.index, entry_black_list=black_list,
                                                                 entry_white_list=white_list)
                self.records = self.records[self.records.index.get_level_values('status').isin(scaffolds_to_keep)]

            if add_row_index:
                self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
                                                               names=("status", "row"))
            print("%s\tReading input finished..." % str(datetime.datetime.now()))

            if sort:
                self.records.sort_values(by=["scaffold", "start", "end"])

    def count_statuses(self):
        status_dict = OrderedDict()
        count_dict = OrderedDict()
        for status in self.status_list:
            status_dict[status] = set(self.records.loc[status]["id"].unique())
            count_dict[status] = len(status_dict[status])

        return status_dict, count_dict

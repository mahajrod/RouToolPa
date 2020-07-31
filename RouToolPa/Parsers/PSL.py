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


class PSL(Parser):

    def __init__(self, in_file=None, records=None, format="psl", parsing_mode="complete",
                 black_list=(), white_list=(), add_row_index=True):

        self.formats = ["psl"]
        self.parsing_parameters = {"psl": {

                                           "all": {
                                                   "col_names": ["matches", # - Number of bases that match that aren't repeats
                                                                 "misMatches", # - Number of bases that don't match
                                                                 "repMatches", # - Number of bases that match but are part of repeats
                                                                 "nCount", # - Number of "N" bases
                                                                 "qNumInsert", # - Number of inserts in query
                                                                 "qBaseInsert", # - Number of bases inserted in query
                                                                 "tNumInsert", # - Number of inserts in target
                                                                 "tBaseInsert", # - Number of bases inserted in target
                                                                 "strand", # - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
                                                                 "qName", # - Query sequence name
                                                                 "qSize", # - Query sequence size.
                                                                 "qStart", # - Alignment start position in query
                                                                 "qEnd", # - Alignment end position in query
                                                                 "tName", #- Target sequence name
                                                                 "tSize", # - Target sequence size
                                                                 "tStart", # - Alignment start position in target
                                                                 "tEnd", # - Alignment end position in target
                                                                 "blockCount", # - Number of blocks in the alignment (a block contains no gaps)
                                                                 "blockSizes", # - Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
                                                                 "qStarts", # - Comma-separated list of starting positions of each block in query
                                                                 "tStarts" #  - Comma-separated list of starting positions of each block in target
                                                                  ],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                 10, 11, 12, 13, 14, 15, 16,
                                                                 17, 18, 19, 20, 21],
                                                   "index_cols": ["qName", "tName"],
                                                   "converters": {
                                                                  "matches":        int,
                                                                  "misMatches":     int,
                                                                  "repMatches":     int,
                                                                  "nCount":         int, # - Number of "N" bases
                                                                  "qNumInsert":     int, # - Number of inserts in query
                                                                  "qBaseInsert":    int, # - Number of bases inserted in query
                                                                  "tNumInsert":     int, # - Number of inserts in target
                                                                  "tBaseInsert":    int, # - Number of bases inserted in target
                                                                  "strand":         str, # - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
                                                                  "qName":          str, # - Query sequence name
                                                                  "qSize":          int, # - Query sequence size.
                                                                  "qStart":         int, # - Alignment start position in query
                                                                  "qEnd":           int, # - Alignment end position in query
                                                                  "tName":          str, #- Target sequence name
                                                                  "tSize":          int, # - Target sequence size
                                                                  "tStart":         int, # - Alignment start position in target
                                                                  "tEnd":           int, # - Alignment end position in target
                                                                  "blockCount":     int, # - Number of blocks in the alignment (a block contains no gaps)
                                                                  "blockSizes":     str, # - Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
                                                                  "qStarts":        str, # - Comma-separated list of starting positions of each block in query
                                                                  "tStarts":        str #  - Comma-separated list of starting positions of each block in target
                                                                  },
                                                   },
                                           "coordinates": {
                                                            "col_names": ["matches",  # - Number of bases that match that aren't repeats
                                                                          "misMatches",  # - Number of bases that don't match
                                                                          "repMatches",  # - Number of bases that match but are part of repeats
                                                                          "nCount",  # - Number of "N" bases
                                                                          "qNumInsert",  # - Number of inserts in query
                                                                          "qBaseInsert",  # - Number of bases inserted in query
                                                                          "tNumInsert",  # - Number of inserts in target
                                                                          "tBaseInsert",  # - Number of bases inserted in target
                                                                          "strand",
                                                                          # - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
                                                                          "qName",  # - Query sequence name
                                                                          "qSize",  # - Query sequence size.
                                                                          "qStart",  # - Alignment start position in query
                                                                          "qEnd",  # - Alignment end position in query
                                                                          "tName",  # - Target sequence name
                                                                          "tSize",  # - Target sequence size
                                                                          "tStart",  # - Alignment start position in target
                                                                          "tEnd",  # - Alignment end position in target
                                                                          "blockCount",  # - Number of blocks in the alignment (a block contains no gaps)

                                                                          ],
                                                            "cols": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                     10, 11, 12, 13, 14, 15, 16,
                                                                     17, 18, ],
                                                            "index_cols": ["qName", "tName"],
                                                            "converters": {
                                                                "matches": int,
                                                                "misMatches": int,
                                                                "repMatches": int,
                                                                "nCount": int,  # - Number of "N" bases
                                                                "qNumInsert": int,  # - Number of inserts in query
                                                                "qBaseInsert": int,  # - Number of bases inserted in query
                                                                "tNumInsert": int,  # - Number of inserts in target
                                                                "tBaseInsert": int,  # - Number of bases inserted in target
                                                                "strand": str, # - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
                                                                "qName": str,  # - Query sequence name
                                                                "qSize": int,  # - Query sequence size.
                                                                "qStart": int,  # - Alignment start position in query
                                                                "qEnd": int,  # - Alignment end position in query
                                                                "tName": str,  # - Target sequence name
                                                                "tSize": int,  # - Target sequence size
                                                                "tStart": int,  # - Alignment start position in target
                                                                "tEnd": int,  # - Alignment end position in target
                                                                "blockCount": int, # - Number of blocks in the alignment (a block contains no gaps)
                                                                          },
                                                           },

                                           },
                                    }

        self.parsing_mode = parsing_mode
        self.format = format

        self.black_list = black_list
        self.white_list = white_list

        # init aliases

        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]

        # load records
        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode, black_list=black_list, white_list=white_list,
                      add_row_index=add_row_index)

        else:
            self.records = records

    def read(self, in_file, format="busco_table", parsing_mode="only_coordinates", sort=False,
             black_list=(), white_list=(), add_row_index=True):
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

        if add_row_index:
            self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
                                                           names=("scaffold", "row"))
        print("%s\tReading input finished..." % str(datetime.datetime.now()))

        #if sort:
        #    self.records.sort_values(by=["scaffold", "start", "end"])


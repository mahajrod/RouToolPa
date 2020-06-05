#!/usr/bin/env python
"""
VCF Parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'

import os
import re
import sys
import datetime

from math import sqrt
from copy import deepcopy

from collections import OrderedDict, Iterable

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from RouToolPa.Collections.General import IdList, IdSet, SynDict
from RouToolPa.Routines import DrawingRoutines
from RouToolPa.Routines.File import FileRoutines
import RouToolPa.Formats.VariantFormats as VariantFormats

ref_alt_variants = {"deaminases": [("C", ["T"]), ("G", ["A"])]
                    }


class MetadataVCF(OrderedDict):
    """
    MetadataVCF class
    """
    def __init__(self, metadata=(), from_file=False, in_file=None,
                 parsing_mode="all", sparse=False):
        OrderedDict.__init__(self)
        self.sparse = sparse
        if from_file:
            self.metadata = []
            self.read(in_file)
        else:
            self.metadata = metadata
        self.converters = OrderedDict()
        self.info_flag_list = []
        self.info_nonflag_list = []
        self.format_flag_list = []
        self.format_nonflag_list = []
        if (metadata or from_file) and (parsing_mode in ("all", "complete")):
            self.create_converters(parsing_mode=parsing_mode, sparse=sparse)
        self.parameter_separator_reg_exp_dict = OrderedDict({
                                                             "GT": re.compile("[\/\|]")
                                                            })
        self.parameter_separator_dict = OrderedDict({
                                                     "GT": "/"
                                                     })
        self.parameter_replace_dict = OrderedDict({
                                                   "GT": {
                                                          ".": np.nan,
                                                          }
                                                   })
        self.default_replace_dict = OrderedDict({
                                                 ".": np.nan
                                                 })
        self.pandas_int_type_correspondence = OrderedDict({
                                                           "Int8": np.float16,
                                                           "Int16": np.float16,
                                                           "Int32": np.float32,
                                                           "Int64": np.float64,
                                                           })

    def create_converters(self, parsing_mode="all", sparse=False):
        if parsing_mode in ("genotypes", "coordinates_and_genotypes"):
            self.converters["FORMAT"] = OrderedDict()
            if sparse:
                self.converters["FORMAT"]["GT"] = pd.SparseDtype(np.int8, fill_value=np.nan)
            else:
                self.converters["FORMAT"]["GT"] = "Int8"
        elif parsing_mode == "pos_gt_dp":
            self.converters["FORMAT"] = OrderedDict()
            if sparse:
                self.converters["FORMAT"]["GT"] = pd.SparseDtype(np.int8, fill_value=np.nan)
                self.converters["FORMAT"]["DP"] = pd.SparseDtype(np.int32, fill_value=np.nan)
            else:
                self.converters["FORMAT"]["GT"] = "Int8"
                self.converters["FORMAT"]["DP"] = "Int32"
        elif parsing_mode in ("all", "complete"):
            for field in "INFO", "FORMAT":
                self.converters[field] = OrderedDict()
                for entry in self[field]:
                    try:
                        a = int(self[field][entry]["Number"])
                    except:
                        a = 2
                    if self[field][entry]["Type"] == "Flag":
                        a = 1
                        self.info_flag_list.append(entry) if field == "INFO" else self.format_flag_list.append(entry)
                    else:
                        self.info_nonflag_list.append(entry) if field == "INFO" else self.format_nonflag_list.append(entry)

                    if entry == "GT":
                        if sparse:
                            self.converters[field][entry] = pd.SparseDtype(np.int8, fill_value=np.nan) if parsing_mode == "complete" else str
                        else:
                            self.converters[field][entry] = "Int8" if parsing_mode == "complete" else str
                    else:
                        if a == 1:
                            if self[field][entry]["Type"] == "Integer":
                                if sparse:
                                    self.converters[field][entry] = pd.SparseDtype(np.int32, fill_value=np.nan)
                                else:
                                    self.converters[field][entry] = "Int32"
                            elif self[field][entry]["Type"] == "Float":
                                if sparse:
                                    self.converters[field][entry] = pd.SparseDtype(np.float32, fill_value=np.nan)
                                else:
                                    self.converters[field][entry] = np.float32
                            elif self[field][entry]["Type"] == "String":
                                self.converters[field][entry] = str
                            elif self[field][entry]["Type"] == "Flag":
                                self.converters[field][entry] = lambda s: True
                            else:
                                raise ValueError("ERROR!!! Unknown value type in metadata: %s, %s, %s" % (field,
                                                                                                          entry,
                                                                                                          self[field][entry]))
                        else:
                            if self[field][entry]["Type"] == "Integer":
                                if sparse:
                                    self.converters[field][entry] = pd.SparseDtype(np.int32, fill_value=np.nan) if parsing_mode == "complete" else str
                                else:
                                    self.converters[field][entry] = "Int32" if parsing_mode == "complete" else str
                            elif self[field][entry]["Type"] == "Float":
                                if sparse:
                                    self.converters[field][entry] = pd.SparseDtype(np.float32, fill_value=np.nan) if parsing_mode == "complete" else str
                                else:
                                    self.converters[field][entry] = np.float32 if parsing_mode == "complete" else str
                            elif self[field][entry]["Type"] == "String":
                                self.converters[field][entry] = str
                            elif self[field][entry]["Type"] == "Flag":
                                self.converters[field][entry] = lambda s: True
                            else:
                                raise ValueError("ERROR!!! Unknown value type in metadata: %s, %s, %s" % (field,
                                                                                                          entry,
                                                                                                          self[field][entry]))

    def read(self, in_file):
        """
        while cycle and readline method are used for compatibility with parsing the data segment of file by pandas

        :param in_file:
        :return:
        """
        #print type(in_file)
        if isinstance(in_file, str):
            with FileRoutines.metaopen(in_file, "r") as fd:
                while True:
                    line = fd.readline()
                    if line[:2] != "##":
                        # self.header = HeaderVCF(line[1:].strip().split("\t"))   # line[1:].strip().split("\t")
                        # self.samples = self.header[9:]
                        if "contig" in self:
                            self["contig"] = pd.DataFrame.from_dict(self["contig"], orient="index")
                            self["contig"].columns = ["length"]
                        return line
                    self.add_metadata(line)
        else:
            while True:
                line = in_file.readline()
                if line[:2] != "##":
                    # self.header = HeaderVCF(line[1:].strip().split("\t"))   # line[1:].strip().split("\t")
                    # self.samples = self.header[9:]
                    if "contig" in self:
                        self["contig"] = pd.DataFrame.from_dict(self["contig"], orient="index")
                        self["contig"].columns = ["length"]
                    return line
                self.add_metadata(line)

    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            # if "=" is not present in string (case of flag type in INFO field)
            return string, None
        return string[:index], string[index+1:]

    @staticmethod
    def _split_by_comma_sign(string):
        index_list = [-1]
        i = 1
        while i < len(string):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == ",":
                index_list.append(i)
            i += 1
        index_list.append(len(string))
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]

    def add_metadata(self, line):
        """
        Adds vcf-like metadata from line
        :param line: string containing metadata info
        :return: None
        """

        key, value = self._split_by_equal_sign(line[2:].strip())
        if key == "contig":
            if "contig" not in self:
                self["contig"] = OrderedDict({})
            value = map(self._split_by_equal_sign, self._split_by_comma_sign(value[1:-1]))
            for entry_tuple in value:
                if entry_tuple[0] == "length":
                    self["contig"][value[0][1]] = int(entry_tuple[1])
                    break

        else:
            if value[0] == "<" and value[-1] == ">":
                value = self._split_by_comma_sign(value[1:-1])
                value_id = self._split_by_equal_sign(value[0])[1]
                value = OrderedDict(self._split_by_equal_sign(entry) for entry in value[1:])
                if key not in self:
                    self[key] = OrderedDict({})
                self[key][value_id] = value
            else:
                self[key] = value

    def add_metadata_from_values(self, name, number, ntype, description):
        """
        Adds vcf-like metadata from values
        :param name: name of parameter
        :param number: number of values in parameter
        :param ntype: type of values in parameter
        :param description: description of parameter
        :return: None
        """
        self[name] = OrderedDict({})
        self[name]["Number"] = number
        self[name]["Type"] = ntype
        self[name]["Description"] = description

    def __str__(self):
        """
        :return: vcf-like string representation of metadata
        """
        metadata_string = ""
        for key in self:
            if key == "contig":
                for contig_tuple in self["contig"].itertuples():
                    metadata_string += "##contig=<ID=%s,length=%i>\n" % (contig_tuple[0], contig_tuple[1])
            elif not isinstance(self[key], dict):
                metadata_string += "##%s=%s\n" % (key, self[key])
            else:

                prefix = "##%s=<" % key
                suffix = ">\n"
                try:
                    for att_id in self[key]:
                        middle = "ID=%s," % att_id + ",".join(["%s=%s" % (param, self[key][att_id][param])
                                                               for param in self[key][att_id]])
                        metadata_string += prefix + middle + suffix
                except:
                    raise ValueError()

        return metadata_string[:-1]


class HeaderVCF(list):
    """
    HeaderVCF class
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SB6536
    """

    def __str__(self):
        """
        :return: vcf-like string representation of header
        """
        return "#" + "\t".join(self)

    def get_samples_list(self):
        return self[9:]


class CollectionVCF:
    """
    CollectionVCF class

    """

    def __init__(self, in_file=None, metadata=None, records=None, header=None, samples=None,
                 external_metadata=None, threads=1, parsing_mode="all", sparse=False,
                 scaffold_black_list=(), scaffold_white_list=(),
                 scaffold_syn_dict=None,):
        """
        Initializes collection. If from_file is True collection will be read from file (arguments other then in_file, external_metadata and threads are ignored)
        Otherwise collection will be initialize from meta, records_dict, header, samples

        IMPORTANT: coordinates are converted to 0-based
        """
        # vcf file columns
        self.formats = ["vcf"]
        self.VCF_COLS = VariantFormats.VCF_COLS

        self.sparse = sparse
        self.parsing_parameters = {
                                   "only_coordinates": {
                                                        "col_names": ["CHROM", "POS"],
                                                        "cols":   [0, 1],
                                                        "index_cols": "CHROM",
                                                        "converters": {
                                                                       "CHROM":  str,
                                                                       "POS":    np.int32,
                                                                       },
                                                        },
                                   "coordinates_and_genotypes": {
                                                                 "col_names": ["CHROM", "POS", "FORMAT"],
                                                                 "cols": [0, 1, 8],
                                                                 "index_cols": "CHROM",
                                                                 "converters": {
                                                                                "CHROM":  str,
                                                                                "POS":    np.int32,
                                                                                "FORMAT": str
                                                                                 },
                                                                 },
                                   "pos_gt_dp": {
                                                 "col_names": ["CHROM", "POS", "FORMAT"],
                                                 "cols": [0, 1, 8],
                                                 "index_cols": "CHROM",
                                                 "converters": {
                                                                "CHROM":  str,
                                                                "POS":    np.int32,
                                                                "FORMAT": str
                                                                },
                                                 },
                                   "except_data":      {
                                                        "col_names": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
                                                        "cols":   [0, 1, 2, 3, 4, 5, 6],
                                                        "index_cols": "CHROM",
                                                        "converters": {
                                                                       "CHROM":  str,
                                                                       "POS":    np.int32,
                                                                       "ID":     str,
                                                                       "REF":    str,
                                                                       "ALT":    str, #lambda s: s.split(","),
                                                                       "QUAL":   self.to_numeric,
                                                                       "FILTER": lambda s: s.split(","),
                                                                       },
                                                        },
                                   "genotypes":        {
                                                        "col_names": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"],
                                                        "cols": [0, 1, 2, 3, 4, 5, 6, 8],
                                                        "index_cols": "CHROM",
                                                        "converters": {
                                                                       "CHROM":  str,
                                                                       "POS":    np.int32,
                                                                       "ID":     str,
                                                                       "REF":    str,
                                                                       "ALT":    str, #lambda s: s.split(","),
                                                                       "QUAL":   self.to_numeric,
                                                                       "FILTER": lambda s: s.split(","),
                                                                       "FORMAT": str #lambda s: s.split(":")
                                                                       },
                                                        },
                                   "pos_ref_alt_id_gt_ad":         {
                                                                    "col_names": ["CHROM", "POS", "ID", "REF", "ALT", "FORMAT"],
                                                                    "cols": [0, 1, 2, 3, 4, 8],
                                                                    "index_cols": "CHROM",
                                                                    "converters": {
                                                                        "CHROM": str,
                                                                        "POS": np.int32,
                                                                        "ID": str,
                                                                        "REF": str,
                                                                        "ALT": str,  # lambda s: s.split(","),
                                                                        "FORMAT": str  # lambda s: s.split(":")
                                                                                    },
                                                                                 },
                                   "read": {
                                        "col_names": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
                                        "cols": None,
                                        "index_cols": "CHROM",
                                        "converters": {
                                                       "CHROM":  str,
                                                       "POS":    np.int32,
                                                       "ID":     str,
                                                       "REF":    str,
                                                       "ALT":    str,
                                                       "QUAL":   str,
                                                       "FILTER": str,
                                                       "INFO":   str,  # self.parse_info_field,
                                                       "FORMAT": str  # lambda s: s.split(":")
                                                                       }
                                    },
                                   "all":              {
                                                        "col_names": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
                                                        "cols": None,
                                                        "index_cols": "CHROM",
                                                        "converters": {
                                                                       "CHROM":  str,
                                                                       "POS":    np.int32,
                                                                       "ID":     str,
                                                                       "REF":    str,
                                                                       "ALT":    str, #lambda s: s.split(","),
                                                                       "QUAL":   self.to_numeric,
                                                                       "FILTER": lambda s: s.split(","),
                                                                       "INFO":   str,   # self.parse_info_field,
                                                                       "FORMAT": str   # lambda s: s.split(":")
                                                                       },
                                                        },

                                   "complete":         {
                                                        "col_names": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
                                                        "cols": None,
                                                        "index_cols": "CHROM",
                                                        "converters": {
                                                                       "CHROM":  str,
                                                                       "POS":    np.int32,
                                                                       "ID":     str,
                                                                       "REF":    str,
                                                                       "ALT":    str, #lambda s: s.split(","),
                                                                       "QUAL":   self.to_numeric,
                                                                       "FILTER": lambda s: s.split(","),
                                                                       "INFO":   str,   # self.parse_info_field,
                                                                       "FORMAT": str   # lambda s: s.split(":")
                                                                       },
                                                        },
                                   }

        self.parsing_modes_with_genotypes = ["complete", "genotypes", "coordinates_and_genotypes", 
                                             "pos_gt_dp"]
        self.parsing_modes_with_sample_coverage = ["complete", "pos_gt_dp"]
        self.parsing_modes_with_allel_coverage = ["pos_ref_alt_id_gt_ad"]
        self.parsing_modes_with_alt_allells = []

        for p_mode in self.parsing_parameters:
            if p_mode != "read":
                if "ALT" in self.parsing_parameters[p_mode]["col_names"]:
                    self.parsing_modes_with_alt_allells.append(p_mode)

        self.linkage_dict = None
        self.parsing_mode = parsing_mode

        self.scaffold_black_list = scaffold_black_list
        self.scaffold_white_list = scaffold_white_list
        self.scaffold_syn_dict = scaffold_syn_dict
        
        if in_file:
            self.read(in_file, external_metadata=external_metadata,
                      parsing_mode=self.parsing_mode, sparse=sparse)
            if "contig" in self.metadata:
                self.scaffold_length = self.metadata["contig"]
        else:
            self.metadata = metadata
            self.records = None if records is None else records
            self.header = header
            self.samples = samples
            self.scaffold_length = None
        self.record_number = len(self.records)
        self.sample_number = len(self.samples)
        self.per_scaffold_record_number = self.records.groupby(self.records.index).size() # pandas Series with scaffold ids as index
        self.scaffold_list = self.records.index.get_level_values('CHROM').unique().to_list()

        self.number_of_scaffolds = len(self.scaffold_list)
        self.threads = threads

    # ========================================= Parsing section=========================================================
    @staticmethod
    def to_numeric(value):
        return pd.to_numeric(value, errors="coerce")

    def read(self, in_file, external_metadata=None, parsing_mode=None, sparse=False):

        """
        Reads collection from vcf file
        :param in_file: path to file
        :param external_metadata: external(not from input file) metadata that could be used to parse records
        :param parsing_mode:
        :param sparse:
        :return: None
        """
        if parsing_mode is not None:
            self.parsing_mode = parsing_mode

        self.metadata = MetadataVCF(sparse=sparse)
        self.records = None

        fd = FileRoutines.metaopen(in_file, "r")

        header_line = self.metadata.read(in_file=fd)
        self.header = HeaderVCF(header_line[1:].strip().split("\t"))
        self.samples = self.header[9:]

        if self.parsing_mode in ("all", "read", "complete", "genotypes", "coordinates_and_genotypes", "pos_gt_dp", "pos_ref_alt_id_gt_ad"):
            self.metadata.create_converters(parsing_mode=self.parsing_mode, sparse=sparse)
            self.parsing_parameters[self.parsing_mode]["col_names"] = self.header
            for sample_col in range(9, 9 + len(self.samples)):
                self.parsing_parameters[self.parsing_mode]["converters"][self.header[sample_col]] = str  # self.parse_sample_field_simple
            if self.parsing_mode in ("genotypes", "coordinates_and_genotypes", "pos_gt_dp", "pos_ref_alt_id_gt_ad"):
                self.parsing_parameters[self.parsing_mode]["cols"] += [i for i in range(9, 9 + len(self.samples))]

        print("%s\tReading file..." % str(datetime.datetime.now()))
        self.records = pd.read_csv(fd, sep='\t', header=None, na_values=".",
                                   usecols=self.parsing_parameters[self.parsing_mode]["cols"],
                                   converters=self.parsing_parameters[self.parsing_mode]["converters"],
                                   names=self.parsing_parameters[self.parsing_mode]["col_names"],
                                   index_col=self.VCF_COLS["CHROM"])
        fd.close()
        print("%s\tReading file finished..." % str(datetime.datetime.now()))

        # convert to 0-based representation

        self.records['POS'] -= 1

        if self.scaffold_white_list and self.scaffold_black_list:
            self.records = self.records[self.records.index.isin(self.scaffold_white_list) & (~self.records.index.isin(self.scaffold_black_list))]
        elif self.scaffold_white_list:
            self.records = self.records[self.records.index.isin(self.scaffold_white_list)]
        elif self.scaffold_black_list:
            self.records = self.records[~self.records.index.isin(self.scaffold_black_list)]

        if self.scaffold_syn_dict:
            self.records.rename(index=self.scaffold_syn_dict, inplace=True)

        self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
                                                       names=("CHROM", "ROW"))
        if parsing_mode in self.parsing_modes_with_alt_allells:
            #print parsing_mode
            #print self.parsing_modes_with_alt_allells
            alt = pd.DataFrame.from_records(map(lambda s: s.split(","), self.records["ALT"].to_list()),
                                            index=self.records.index).astype("category")
            alt_colomn_number = np.shape(alt)[1] if len(np.shape(alt)) > 1 else 1
            alt.columns = pd.MultiIndex.from_arrays([["ALT"] * alt_colomn_number,
                                                     [i for i in range(0, alt_colomn_number)]])
        if "REF" in self.records:
            self.records["REF"] = self.records["REF"].astype("category")

        if self.parsing_mode in ("genotypes", "coordinates_and_genotypes"):
            sample_genotypes = self.parse_samples(["GT"])
            self.records.columns = pd.MultiIndex.from_arrays([
                                                              self.records.columns,
                                                              self.records.columns,
                                                              self.records.columns
                                                              ])
            if self.parsing_mode == "coordinates_and_genotypes":
                self.records = self.records[["POS"]]
                self.records.columns = pd.MultiIndex.from_arrays([["POS"], ["POS"], ["POS"]])
                self.records = pd.concat([self.records] + sample_genotypes, axis=1)
            else:
                self.records = self.records[["POS", "ID", "REF", "QUAL", "FILTER"]],
                self.records.columns = pd.MultiIndex.from_arrays([self.records.columns,
                                                                  self.records.columns,
                                                                  self.records.columns
                                                                  ])
                alt.columns = pd.MultiIndex.from_arrays([["ALT"] * alt_colomn_number,
                                                         ["ALT"] * alt_colomn_number,
                                                         [i for i in range(0, alt_colomn_number)]])

                self.records = pd.concat([self.records, alt] + sample_genotypes, axis=1)

        elif parsing_mode == "pos_gt_dp":
            sample_data = self.parse_samples(["GT", "DP"])
            self.records.columns = pd.MultiIndex.from_arrays([
                                                              self.records.columns,
                                                              self.records.columns,
                                                              self.records.columns
                                                              ])
            self.records = pd.concat([self.records[["POS"]]] + sample_data, axis=1)
            
        elif parsing_mode == "pos_ref_alt_id_gt_ad":
            sample_data = self.parse_samples(["GT", "AD"])
            self.records.columns = pd.MultiIndex.from_arrays([
                                                              self.records.columns,
                                                              self.records.columns,
                                                              ])
            self.records = pd.concat([self.records[["POS", "REF", "ALT", "ID"]]] + sample_data, axis=1)
        elif self.parsing_mode in ("all", "complete"):
            info = self.parse_info()
            sample_list = self.parse_samples()
            self.records = self.records[["POS", "ID", "REF", "QUAL", "FILTER"]]
            if self.parsing_mode == "all":
                self.records.columns = pd.MultiIndex.from_arrays([
                                                                  self.records.columns,
                                                                  self.records.columns
                                                                  ])
            elif self.parsing_mode == "complete":
                alt.columns = pd.MultiIndex.from_arrays([["ALT"] * alt_colomn_number,
                                                         ["ALT"] * alt_colomn_number,
                                                         [i for i in range(0, alt_colomn_number)]])
                self.records.columns = pd.MultiIndex.from_arrays([
                                                                  self.records.columns,
                                                                  self.records.columns,
                                                                  self.records.columns
                                                                  ])
            self.records = pd.concat([self.records, alt,  info] + sample_list, axis=1)

    def parse_column(self, column, param, param_group):
        if (self.parsing_mode == "all") or (self.parsing_mode == "pos_ref_alt_id_gt_ad"):
            if self.metadata.converters[param_group][param] == str:
                return column
            col = column.replace(self.metadata.default_replace_dict)
            if self.sparse:
                col = col.astype(self.metadata.converters[param_group][param])
            else:
                if self.metadata.converters[param_group][param] in self.metadata.pandas_int_type_correspondence:
                    col = col.apply(self.metadata.pandas_int_type_correspondence[self.metadata.converters[param_group][param]]).astype(self.metadata.converters[param_group][param])
                else:
                    col = col.apply(self.metadata.converters[param_group][param])
        elif self.parsing_mode in ("complete", "genotypes", "coordinates_and_genotypes", "pos_gt_dp"):
            index = column.index
            if param in self.metadata.parameter_separator_reg_exp_dict:
                col = pd.DataFrame(map(self.metadata.parameter_separator_reg_exp_dict[param].split, column.to_list()))
            else:
                #print param
                #print column
                #print map(lambda s: s.split(self.metadata.parameter_separator_dict[param] if param in self.metadata.parameter_separator_dict else ",") if type(s) == str else [s], column.to_list())
                col = pd.DataFrame.from_records(map(lambda s: s.split(self.metadata.parameter_separator_dict[param] if param in self.metadata.parameter_separator_dict else ",") if type(s) == str else [s], column.to_list()))
                #col = column.str.split(self.metadata.parameter_separator_dict[param] if param in self.metadata.parameter_separator_dict else ",", expand=True)
                #print col
            col.index = index
            if self.metadata[param_group][param]["Type"] != "Flag":
                col.replace(self.metadata.default_replace_dict, inplace=True)
            if self.metadata.converters[param_group][param] == str:
                return col
            if self.sparse:
                col = col.astype(self.metadata.converters[param_group][param])
            else:
                if self.metadata.converters[param_group][param] in self.metadata.pandas_int_type_correspondence:
                    #print self.metadata.converters[param_group][param]
                    #print col
                    col = col.apply(self.metadata.pandas_int_type_correspondence[self.metadata.converters[param_group][param]]).astype(self.metadata.converters[param_group][param])

                else:
                    col = col.apply(self.metadata.converters[param_group][param])
        else:
            raise ValueError("ERROR!!! Unknown parsing mode: %s" % self.parsing_mode)

        return col
    
    def parse_info(self):
        print("%s\tParsing info field..." % str(datetime.datetime.now()))

        def split_info_entry(info_entry):
            tmp = info_entry.split("=")
            if len(tmp) == 2:
                return tmp
            else:
                return tmp[0], True

        tmp_info = pd.DataFrame(map(lambda s: OrderedDict(map(split_info_entry, s.split(";"))), list(self.records["INFO"])))
        tmp_info.index = self.records.index

        #print tmp_info
        info_df_list = []
        for param in self.metadata.info_flag_list + self.metadata.info_nonflag_list:
            #print self.metadata.info_flag_list
            #print self.metadata.info_nonflag_list
            if param in tmp_info:
                print("%s\tParsing info parameter %s..." % (str(datetime.datetime.now()), param))
                if param in self.metadata.info_nonflag_list:
                    column_df = self.parse_column(tmp_info[param], param, "INFO")
                else:
                    column_df = deepcopy(tmp_info[[param]])
                shape = np.shape(column_df)
                column_number = 1 if len(shape) == 1 else shape[1]
                # print column_number
                if self.parsing_mode == "all":
                    column_df.columns = pd.MultiIndex.from_arrays([
                                                  ["INFO"] * column_number,
                                                  [param] * column_number
                                                  ])
                if self.parsing_mode == "complete":

                    column_df.columns = pd.MultiIndex.from_arrays([
                                                     ["INFO"] * column_number,
                                                     [param] * column_number,
                                                     np.arange(0, column_number)
                                                     ])

                info_df_list.append(column_df)
        #print info_df_list
        info = pd.concat(info_df_list, axis=1)
        info.sort_index(level=1, inplace=True)

        del info_df_list
        #print info

        print("%s\tParsing info field finished..." % str(datetime.datetime.now()))
        return info

    def parse_samples(self, parameter_list=()):
        print("%s\tParsing samples..." % str(datetime.datetime.now()))
        uniq_format_set = self.records['FORMAT'].drop_duplicates()
        uniq_format_dict = OrderedDict([(format_entry, format_entry.split(":")) for format_entry in uniq_format_set])
        sample_data_dict = {}

        present_parameter_dict = OrderedDict()

        for format_entry in uniq_format_dict:
            present_parameter_dict[format_entry] = []
            if parameter_list:
                for parameter in parameter_list:
                    if parameter in uniq_format_dict[format_entry]:
                        #print parameter, uniq_format_dict[format_entry]
                        present_parameter_dict[format_entry].append(parameter)
            else:
                present_parameter_dict[format_entry] = uniq_format_dict[format_entry]

        for sample in self.samples:
            sample_data_dict[sample] = OrderedDict()
            for format_entry in uniq_format_dict:
                sample_data_dict[sample][format_entry] = list()
                print(self.records[self.records['FORMAT'] == format_entry])
                tmp = self.records[self.records['FORMAT'] == format_entry][sample]
                tmp_index = deepcopy(tmp.index)
                tmp = pd.DataFrame(map(lambda s: s.split(":"), list(tmp)))
                tmp.index = tmp_index
                tmp.columns = uniq_format_dict[format_entry]
                sample_data_dict[sample][format_entry] = []

                for parameter in present_parameter_dict[format_entry] if parameter_list else uniq_format_dict[format_entry]:
                    print(tmp[parameter])
                    parameter_col = self.parse_column(tmp[parameter], parameter, "FORMAT")
                    sample_data_dict[sample][format_entry].append(parameter_col)

                for i in range(0, len(present_parameter_dict[format_entry])) if parameter_list else range(0, len(uniq_format_dict[format_entry])):
                    shape = np.shape(sample_data_dict[sample][format_entry][i])
                    column_number = 1 if len(shape) == 1 else shape[1]
                    if (self.parsing_mode == "all") or (self.parsing_mode == "pos_ref_alt_id_gt_ad"):
                        column_index = pd.MultiIndex.from_arrays([
                                                                  [sample] * column_number,
                                                                  [uniq_format_dict[format_entry][i]] * column_number
                                                                  ],)
                    elif self.parsing_mode in self.parsing_modes_with_genotypes:

                        column_index = pd.MultiIndex.from_arrays([
                                                                  [sample] * column_number,
                                                                  [present_parameter_dict[format_entry][i]] * column_number,
                                                                  np.arange(0, column_number)
                                                                  ],)
                    """
                    elif self.parsing_mode in ("genotypes", "coordinates_and_genotypes", "pos_gt_dp"):
                        column_index = pd.MultiIndex.from_arrays([
                                                                  [sample] * column_number,
                                                                  [present_parameter_dict[format_entry][i]] * column_number,
                                                                  np.arange(0, column_number)
                                                                  ],)
                    """
                    sample_data_dict[sample][format_entry][i].columns = column_index
                if sample_data_dict[sample][format_entry]:
                    #print sample_data_dict[sample][format_entry]
                    sample_data_dict[sample][format_entry] = pd.concat(sample_data_dict[sample][format_entry],
                                                                       axis=1)
            if sample_data_dict[sample]:
                sample_data_dict[sample] = pd.concat(sample_data_dict[sample].values(),
                                                     axis=0)
                if (self.parsing_mode == "all") or (self.parsing_mode == "pos_ref_alt_id_gt_ad"):
                    column_index = pd.MultiIndex.from_arrays([
                                                              [sample] * len(sample_data_dict[sample].columns),
                                                              sample_data_dict[sample].columns
                                                              ],)
                    sample_data_dict[sample].columns = column_index
                # sort by row number
                sample_data_dict[sample].sort_index(level=1, inplace=True)
            else:
                sample_data_dict.pop(sample, None)
        print("%s\tParsing sample finished..." % str(datetime.datetime.now()))
        return list(sample_data_dict.values())

    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            return string, None
        return string[:index], string[index+1:]

    def _split_by_comma_sign(self, string):
        return self._split_by_sign(string, sign=",")

    @staticmethod
    def _split_by_sign(string, sign=","):
        #IMPORTANT!!! ignores sign in "
        index_list = [-1]
        i = 1
        while (i < len(string)):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == sign:
                index_list.append(i)
            i += 1
        index_list.append(len(string))
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]
    # ========================================= Parsing section end ====================================================

    # ============================================ Writing section =====================================================

    def write(self, outfile, format='simple_bed', bed_type="0-based", samples=None, split_samples=False):
        def not_ref_variant(s):
            return (s[:3] != "0/0") and (s[:3] != "./.")

        df = self.records.reset_index(level='CHROM')
        if format == 'simple_bed':
            if bed_type == '1-based':
                df["POS"] += 1
            df[["CHROM", "POS"]].to_csv(outfile, sep='\t', index=False, header=False)

        elif format == 'bed':
            if bed_type == '1-based':
                df["POS"] += 1
            df.to_csv(outfile, sep='\t', index=False, header=False)

        elif format == "vcf":
            df["POS"] += 1
            #print split_samples
            if split_samples:
                for sample in samples if samples else self.samples:
                    with open("%s.%s.vcf" % (outfile, sample), "w") as out_fd:
                        out_fd.write(str(self.metadata))
                        out_fd.write("\n")

                        header = self.header[:9] + [sample]
                        #print sample
                        #print header
                        out_fd.write("#" + "\t".join(header) + "\n")

                        if self.parsing_mode == "read":
                            df[df[[sample]].applymap(not_ref_variant)[sample]][header].to_csv(out_fd, sep="\t", header=False, index=False, columns=header)
                        elif self.parsing_mode == "all":
                            pass
                        elif self.parsing_mode == "complete":
                            pass

                        out_fd.close()

            else:
                with open(outfile, "w") as out_fd:
                    out_fd.write(str(self.metadata))
                    out_fd.write("\n")

                    if samples:
                        if isinstance(samples, str):
                            df = df[df[[samples]].applymap(not_ref_variant)[samples]]
                            header = self.header[:9] + [samples]
                        else:
                            header = self.header[:9] + samples
                            df = df[np.logical_or.reduce([df[[sample]].applymap(not_ref_variant)[sample] for sample in samples])]
                        out_fd.write("#" + "\t".join(header))
                    else:
                        header = self.header
                        out_fd.write(str(self.header))

                    out_fd.write("\n")

                    if self.parsing_mode == "read":
                        df.to_csv(out_fd, sep="\t", header=False, index=False, columns=header)

                    out_fd.close()

    @staticmethod
    def write_df(dataframe, outfile, format='simple_bed', type="0-based"):
        if format == 'simple_bed':
            if type == "0-based":
                dataframe[["POS"]].reset_index(level='CHROM').to_csv(outfile, sep='\t', index=False, header=False)
            elif type == '1-based':
                df = dataframe[["POS"]].reset_index(level='CHROM')
                df["POS"] += 1
                df.to_csv(outfile, sep='\t', index=False, header=False)
        elif format == 'bed':
            if type == "0-based":
                dataframe.reset_index(level='CHROM').to_csv(outfile, sep='\t', index=False, header=False)
            elif type == '1-based':
                df = dataframe.reset_index(level='CHROM')
                df["POS"] += 1
                df.to_csv(outfile, sep='\t', index=False, header=False)

    # ========================================= Writing section end ====================================================

    # =========================================== Filtering section ====================================================
    def check_records_by_expression(self, expression):
        """
        Checks records in collection by expression. Expression must be a function with one argument - record entry,
        returning boolean
        :param expression:  expression to check
        :return: array of booleans with length of record number.
        """
        return self.records.apply(expression, axis=1)

    def filter_records(self, expression):
        boolean_array = self.records.apply(expression, axis=1)

        return self.records[boolean_array], self.records[~boolean_array]

    def filter(self, expression):
        """
        Splits collection based on expression. Expression must be a function with one argument - record entry
        :param expression: filtering expression
        :return: tuple of two CollectionVCF. First contains records for which expression is True, second - False.
        """
        #
        filtered_records, filtered_out_records = self.filter_records(expression)
        return CollectionVCF(metadata=self.metadata, records=filtered_records,
                             header=self.header, samples=self.samples, ),\
               CollectionVCF(metadata=self.metadata, records=filtered_out_records,
                             header=self.header, samples=self.samples, )

    # ========================================= Filtering section end ==================================================
    def count_records(self, expression):
        """
        Counts records passed expression in collection based on expression.
        Expression must be a function with one argument - record entryreturning boolean
        :param expression: filtering expression
        :return: tuple of two numbers. First is number of records for which expression is True, second - False.
        """
        return np.sum(self.check_records_by_expression(expression))

    def rainfall_plot(self, plot_name, dpi=300, figsize=(20, 20), facecolor="#D6D6D6",
                      ref_genome=None, min_masking_length=10, suptitle=None,
                      masking_color="#777777", logbase=2,
                      extension_list=("pdf", "png"),
                      scaffold_black_list=None, scaffold_white_list=None,
                      scaffold_ordered_list=None, sort_scaffolds=False,
                      color_expression=None,
                      default_point_color='blue',
                      dot_size=None,
                      label_fontsize=None, draw_masking=False):
        """

        :param plot_name:
        :param base_colors:
        :param single_fig:
        :param dpi:
        :param figsize:
        :param facecolor:
        :param ref_genome:
        :param masked_scaffolds:
        :param min_gap_length:
        :param draw_gaps:
        :param suptitle:
        :param gaps_color:
        :param masked_scaffolds_color:
        :param logbase:
        :param extension_list:
        :param scaffold_black_list:
        :param scaffold_white_list=:
        :param scaffold_order_list=None
        :return:

        """
        # TODO: add multithreading drawing if possible and multipicture drawing
        print("Drawing rainfall plot...")
        plot_dir = "rainfall_plot"

        os.system("mkdir -p %s" % plot_dir)

        fig = plt.figure(1, dpi=dpi, figsize=figsize ) #, facecolor=facecolor)
        fig.suptitle(suptitle if suptitle else "Rainfall plot", fontweight='bold', y=0.94, fontsize=label_fontsize) #
        sub_plot_dict = OrderedDict({})
        index = 1

        final_scaffold_list = DrawingRoutines.get_filtered_scaffold_list(self.scaffold_list,
                                                                         scaffold_black_list=scaffold_black_list,
                                                                         sort_scaffolds=sort_scaffolds,
                                                                         scaffold_ordered_list=scaffold_ordered_list,
                                                                         scaffold_white_list=scaffold_white_list,
                                                                         sample_level=False)
        num_of_scaffolds = len(final_scaffold_list)
        distances_dict = OrderedDict()
        height = 0

        if (ref_genome is not None) and draw_masking:
            masking_df = ref_genome.get_merged_gaps_and_masking()
            if min_masking_length > 1:
                masking_df.remove_small_records(min_masking_length)

        for scaffold in final_scaffold_list: # self.records
            print("Handling scaffold: %s ..." % scaffold)
            distances_dict[scaffold] = self.records.loc[scaffold, "POS"].diff()
            height = max(np.max(distances_dict[scaffold]), height)
            # pandas DataFrame diff methods return differences between consecutive elements in array,
            # and first distance is NaN always, so it is replaced by 0
            distances_dict[scaffold][0] = 0
            distances_dict[scaffold].name = 'DIST'
            if color_expression:
                colors = self.records.loc[scaffold].apply(color_expression, axis=1)
                colors.name = 'COLOR'
                distances_dict[scaffold] = pd.concat([self.records.loc[scaffold, "POS"],
                                                      distances_dict[scaffold],
                                                      colors],
                                                     axis=1)
                distances_dict[scaffold] = distances_dict[scaffold].set_index('COLOR')
                color_list = colors.index.values.unique().to_list()
            else:
                distances_dict[scaffold] = pd.concat([self.records.loc[scaffold, "POS"],
                                                      distances_dict[scaffold]],
                                                     axis=1)

        length = np.max(ref_genome.seq_lengths['length']) if ref_genome is not None else np.max(self.records["POS"])

        length *= 1.1
        if length // (10 ** 9) > 2:
            def tick_formater(x, pos):
                return '%1.1f Gbp' % (x*1e-9)
        elif length // (10 ** 6) > 200:
            def tick_formater(x, pos):
                return '%.0f Mbp' % (x*1e-6)
        elif length // (10 ** 6) > 2:
            def tick_formater(x, pos):
                return '%.1f Mbp' % (x*1e-6)

        formatter = FuncFormatter(tick_formater)

        for scaffold in final_scaffold_list:
            if not sub_plot_dict:
                sub_plot_dict[scaffold] = plt.subplot(num_of_scaffolds, 1, index) #, axisbg=facecolor)

            else:
                sub_plot_dict[scaffold] = plt.subplot(num_of_scaffolds, 1, index,
                                                      sharey=sub_plot_dict[final_scaffold_list[0]])
                                                      #sharex=sub_plot_dict[keys[0]],
                                                      #)
                                                      #facecolor=facecolor)
            sub_plot_dict[scaffold].xaxis.set_major_formatter(formatter)
            index += 1

            if ref_genome is not None:
                print("\tScaffold length:%i" % ref_genome.seq_lengths.loc[scaffold])
                plt.gca().add_patch(plt.Rectangle((1, 0),
                                                  ref_genome.seq_lengths.loc[scaffold],
                                                  height, facecolor=facecolor, edgecolor='none', alpha=0.5))
                if draw_masking:
                    for masked_region in masking_df.records.loc[scaffold].itertuples(index=False):
                        plt.gca().add_patch(plt.Rectangle((masked_region[0] + 1, 1),
                                                          masked_region[1] - masked_region[0],
                                                          height, facecolor=masking_color, edgecolor='none'))

            print("Drawing scaffold: %s ..." % scaffold)

            if color_expression:
                for color in color_list:
                    plt.scatter(distances_dict[scaffold].loc[color]['POS'],
                                distances_dict[scaffold]['DIST'],
                                color=color,
                                marker='.', s=dot_size)
            else:
                #print distances_dict[scaffold]
                #print distances_dict[scaffold]['POS']
                #print distances_dict[scaffold]['DIST']
                #print "UUUUUUUU"
                plt.scatter(distances_dict[scaffold]['POS'],
                            distances_dict[scaffold]['DIST'],
                            color=default_point_color,
                            marker='.', s=dot_size)

            plt.text(-0.13, 0.5, scaffold, rotation=0, fontweight="bold", transform=sub_plot_dict[scaffold].transAxes,
                     fontsize=label_fontsize,
                     horizontalalignment='center',
                     verticalalignment='center')
            plt.ylabel("Distanse")
            #plt.axhline(y=100, color="#000000")
            #plt.axhline(y=1000, color="#000000")
            #plt.axhline(y=500, color="purple")
            #plt.axhline(y=10, color="#000000")
            sub_plot_dict[scaffold].set_yscale('log', basey=logbase)
            sub_plot_dict[scaffold].get_xaxis().set_visible(False)
            sub_plot_dict[scaffold].spines['right'].set_color('none')
            sub_plot_dict[scaffold].spines['top'].set_color('none')
            sub_plot_dict[scaffold].spines['bottom'].set_color('none')
            plt.xlim(xmin=1, xmax=length)
            #plt.ylim(ymax=height)
            #plt.tight_layout()
        #sub_plot_dict[scaffold].unshare_x_axes(sub_plot_dict[first_scaffold])
        sub_plot_dict[final_scaffold_list[-1]].get_xaxis().set_visible(True)
        sub_plot_dict[scaffold].spines['bottom'].set_color('black')
        #plt.ylim(ymax=max_distance * 1.10)
        plt.subplots_adjust(left=0.175, bottom=0.05, right=0.95, top=0.90, wspace=None, hspace=None)
        for extension in extension_list:
            plt.savefig("%s/%s_log_scale.%s" % (plot_dir, plot_name, extension))
        plt.close()

    def check_variant_presence(self, outfile=None):
        if self.parsing_mode in self.parsing_modes_with_genotypes:

            variant_presence = pd.concat([((self.records[sample]["GT"][0].notna()) & (self.records[sample]["GT"][0] != 0)) | ((self.records[sample]["GT"][1].notna()) & (self.records[sample]["GT"][1] != 0)) for sample in self.samples], axis=1)
            variant_presence.columns = self.samples
            if outfile:
                variant_presence.to_csv(outfile, sep="\t", header=True, index=True)
            return variant_presence
        else:
            raise ValueError("ERROR!!! Variant presence can't be counted for this parsing mode: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete modes'" % self.parsing_mode)

    def get_uniq_variants(self, output_prefix):
        variant_presence = self.check_variant_presence(outfile="%s.variant_presence" % output_prefix)
        return variant_presence[variant_presence.apply(lambda s: True if np.sum(s) == 1 else False, axis=1)]

    def count_uniq_variants(self, output_prefix, extension_list=("png",), figsize=(5, 5), dpi=200,
                            title="Unique variants"):
        if self.parsing_mode in self.parsing_modes_with_genotypes:

            variant_presence = pd.concat([((self.records[sample]["GT"][0].notna()) & (self.records[sample]["GT"][0] != 0)) | ((self.records[sample]["GT"][1].notna()) & (self.records[sample]["GT"][1] != 0)) for sample in self.samples], axis=1)
            variant_presence.columns = self.samples
            uniq_variants = variant_presence[variant_presence.apply(lambda s: True if np.sum(s) == 1 else False, axis=1)]
            uniq_variant_counts = uniq_variants.apply(np.sum)

            if output_prefix:
                #variant_presence.to_csv("%s.variant_presence" % output_prefix, sep="\t", header=True, index=True)
                #uniq_variants.to_csv("%s.uniq_variants" % output_prefix, sep="\t", header=True, index=True)
                uniq_variant_counts.to_csv("%s.uniq_variants.counts" % output_prefix, sep="\t", header=True, index=True)

            fig = plt.figure(1, figsize=figsize, dpi=dpi)

            bar_width = 0.5
            bin_coord = np.arange(len(self.samples))

            plt.bar(bin_coord, uniq_variant_counts, width=bar_width, edgecolor='white', color='blue',)

            plt.ylabel('Variants', fontweight='bold')
            plt.xlabel('Sample', fontweight='bold')
            plt.xticks(bin_coord, self.samples, rotation=45)
            plt.title(title, fontweight='bold')

            for extension in extension_list:
                plt.savefig("%s.%s" % (output_prefix, extension), bbox_inches='tight')
            plt.close()

            return uniq_variant_counts
        else:
            raise ValueError("ERROR!!! Variant presence can't be counted for this parsing mode: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete modes'" % self.parsing_mode)

    def draw_sample_parameter_distribution(self, parameter, bin_width, output_prefix=None,
                                           extension_list=("png",), suptitle=None,
                                           xlabel=None, ylabel=None, show_median=True,
                                           show_mean=True, median_relative=False, mean_relative=False, dpi=200,
                                           subplot_size=3, xlimit=None, verbose=False, ylogbase=10):

        param = self.records.xs(parameter, axis=1, level=1, drop_level=False)
        param_mean = param.apply(np.mean)
        param_median = param.apply(np.median)
        if verbose:
            print("Median:")
            print(param_median)
            print("Mean:")
            print(param_mean)

        if median_relative:
            param = param.astype(np.float32) / param_median
            param_mean = param_mean.astype(np.float32) / param_median
            param_median = param_median.astype(np.float32) / param_median
        elif mean_relative:
            param = param.astype(np.float32) / param_mean
            param_mean = param_mean.astype(np.float32) / param_mean
            param_median = param_median.astype(np.float32) / param_mean

        param_max = param.apply(np.max)
        param_min = param.apply(np.min)

        # selection of figure size
        n = int(np.sqrt(self.sample_number))
        if n * (n + 1) >= self.sample_number:
            m = n + 1
        else:
            n = n +1
            m = n
        if median_relative or mean_relative:
            if param_max > max(param_median) * 10:
                bins = np.arange(0, max(param_median) * 10, bin_width)
                bins = np.concat(bins, [max(param_max)])
            else:
                bins = np.arange(0, max(param_max), 0.1)
        else:
            print(np.max(param_median))
            print(np.max(param_median)[0])
            print(max(param_median))
            if param_max > max(param_median) * 10:
                bins = np.arange(1, max(param_median) * 10, bin_width)
                bins = np.concat(bins, [max(param_max)])
            else:
                bins = np.arange(1, max(param_max), bin_width)
        bins = np.concatenate((bins, [bins[-1] + bin_width, bins[-1] + 2 * bin_width]))

        print("Bins:")
        print(bins)

        figure, subplot_array = plt.subplots(nrows=n, ncols=m, sharex=True, sharey=True,
                                             figsize=(m*subplot_size, n*subplot_size), dpi=dpi)
        #print subplot_array
        #print np.shape(subplot_array)
        #print n, m
        for row in range(0, n):
            for col in range(0, m):
                #print row, col
                sample_index = row * m + col
                if ylabel and col == 0:
                    subplot_array[row][col].ylabel = ylabel
                if xlabel and row == n - 1:
                    subplot_array[row][col].xlabel = xlabel

                if sample_index >= self.sample_number:
                    continue
                sample_id = self.samples[sample_index]
                #print param[sample_id]
                # TODO: adjust function to deal not only with the first column inside parameter
                subplot_array[row][col].hist(param[sample_id][parameter][0].dropna(), bins=bins, label=sample_id)
                if show_median:
                    subplot_array[row][col].axvline(x=float(param_median[sample_id]), label="median %.2f" % param_median, color="orange")
                if show_mean:
                    subplot_array[row][col].axvline(x=float(param_mean[sample_id]), label="mean %.2f" % param_mean, color="red")
                if row == 0 and col == m - 1:
                    subplot_array[row][col].legend()
                subplot_array[row][col].set_title(sample_id)
        if suptitle:
            supt = suptitle
        elif mean_relative:
            supt = "%s distribution(Mean relative)" % parameter
        elif median_relative:
            supt = "%s distribution(Median relative)" % parameter
        else:
            supt = "%s distribution" % parameter
        plt.xlim(xmin=0)
        plt.suptitle(supt)

        if output_prefix:
            for extension in extension_list:
                plt.savefig("%s.%s" % (output_prefix, extension), bbox_inches='tight')

        xlim = xlimit if xlimit else np.max(param_median)*3
        plt.xlim(xmax=xlim, xmin=0)
        if output_prefix:
            for extension in extension_list:
                plt.savefig("%s.xlim%i.%s" % (output_prefix, xlim, extension), bbox_inches='tight')
            plt.yscale('log', basey=ylogbase)
            for extension in extension_list:
                plt.savefig("%s.xlim%i.ylog.%s" % (output_prefix, xlim, extension), bbox_inches='tight')

        plt.close()

        return param

    def get_coverage_distribution(self, output_prefix, bin_width=5, dpi=200, subplot_size=3, extension_list=("png",),
                                  verbose=False):
        if self.parsing_mode in self.parsing_modes_with_sample_coverage:
            print("Drawing coverage distribution...")
            self.draw_sample_parameter_distribution("DP", bin_width, output_prefix=output_prefix,
                                                    extension_list=extension_list,
                                                    suptitle="Coverage distribution",
                                                    xlabel="Coverage", ylabel="Variants", show_median=True,
                                                    show_mean=True, median_relative=False, mean_relative=False,
                                                    dpi=dpi, subplot_size=subplot_size,
                                                    verbose=verbose)
            print("Drawing coverage distribution relative to median...")
            self.draw_sample_parameter_distribution("DP", bin_width, output_prefix="%s.median_relative" % output_prefix,
                                                    extension_list=extension_list,
                                                    suptitle="Coverage distribution(Median relative)",
                                                    xlabel="Coverage", ylabel="Variants", show_median=True,
                                                    show_mean=True, median_relative=True, mean_relative=False,
                                                    dpi=dpi, subplot_size=subplot_size)
            print("Drawing coverage distribution relative to mean...")
            self.draw_sample_parameter_distribution("DP", bin_width, output_prefix="%s.mean_relative" % output_prefix,
                                                    extension_list=extension_list,
                                                    suptitle="Coverage distribution(Mean relative)",
                                                    xlabel="Coverage", ylabel="Variants", show_median=True,
                                                    show_mean=True, median_relative=False, mean_relative=True,
                                                    dpi=dpi, subplot_size=subplot_size)
        else:
            raise ValueError("ERROR!!! Coverage distribution can't be counted for this parsing mode: %s."
                             "Use 'pos_gt_dp' or other method parsing DP column from samples fields" % self.parsing_mode)

    def calculate_masking(self, outfile, samples=None, sample_coverage=None, min_samples=1, max_coverage=2.5, min_coverage=None):
        if self.parsing_mode in self.parsing_modes_with_sample_coverage:
            samples_to_use = samples if samples else self.samples
            coverage = self.records[samples_to_use].xs("DP", axis=1, level=1, drop_level=False)
            if sample_coverage:
                sp_coverage = pd.Series(sample_coverage, dtype=np.float32)
                sp_coverage.index = pd.MultiIndex.from_arrays([samples_to_use,
                                                               ["DP"] * len(samples_to_use),
                                                               [0] * len(samples_to_use)])
            else:
                sp_coverage = coverage.apply(np.median)
            #coverage = coverage / coverage_median
            #print sp_coverage
            #print "UUUU"
            #print coverage.apply(np.median)
            boolean_array = coverage >= (max_coverage * sp_coverage)
            if min_coverage:
                boolean_array &= coverage <= (min_coverage * sp_coverage)

            outliers = boolean_array.apply(np.sum, axis=1)
            outliers = outliers[outliers >= min_samples]
            #outliers = pd.concat([self.records[self.records.index.isin(outliers.index)]["POS"], outliers], axis=1)
            outliers = self.records[self.records.index.isin(outliers.index)]["POS"]

            print("%i variants were masked" % np.shape(outliers)[0])

            self.write_df(outliers, outfile, format="simple_bed", type="1-based")

        else:
            raise ValueError("ERROR!!! Masking can't be counted for this parsing mode: %s."
                             "Use 'pos_gt_dp' or other method parsing DP column from samples fields" % self.parsing_mode)

    #########################################################################
    #                        In progress                                    #
    #########################################################################

    #########################################################################
    # methods below were not yet rewritten for compatibility with VCFpandas #
    #########################################################################

    def no_reference_allel_and_multiallel(self, record, sample_index=None, max_allels=None):
        return record.no_reference_allel_and_multiallel(sample_index=sample_index, max_allels=max_allels)

    def filter_variants_with_reference_allel_and_multiallelic(self, sample_index=None, max_allels=None):

        def expression(record):
            return self.no_reference_allel_and_multiallel(record, sample_index=sample_index, max_allels=max_allels)

        return self.filter(expression)

    @staticmethod
    def filter_zygoty_expression(record):
            for sample_dict in record.samples_list:
                zyg = sample_dict["GT"][0].split("/")
                if zyg[0] != zyg[1]:
                    return False
            return True

    def filter_by_zygoty(self):
        """
        Splits collection based on zygoty of mutation. Mutation is counted as heterozygous even if in one sample it is hetorozygous
        :return: tuple of two CollectionVCF. First contains homozygous records, second - heterozygous
        """
        """
        def filter_expression(record):
            for sample_dict in record.samples_list:
                zyg = sample_dict["GT"][0].split("/")
                if zyg[0] != zyg[1]:
                    return False
            return True
        """
        return self.filter(self.filter_zygoty_expression)

    @staticmethod
    def filter_by_filter_presence_expression(record):
        #print record.filter_list
        for filter_entry in record.filter_list:
            #print filter_entry
            if (filter_entry != "PASS") and (filter_entry != "."):
                #print "FALSE"
                return False
        #print True
        return True

    def filter_by_filter_presence(self):
        return self.filter(self.filter_by_filter_presence_expression)

    def record_coordinates(self, black_list=[], white_list=[]):
        """
        Extracts coordinates of records in collection
        :param black_list: list of scaffolds to skip
        :param white_list: list of scaffolds to consider, other are ignored
        :return: dictionary of coordinates, keys are scaffold names, values - Numpy arrays of coordinates
        """
        coord_dict = {}
        for scaffold in self.records:
            for record in self.records[scaffold]:
                if black_list and (scaffold in black_list):
                    continue
                if white_list and (scaffold not in white_list):
                    continue
                if scaffold not in coord_dict:
                    coord_dict[scaffold] = [record.pos]
                else:
                    coord_dict[scaffold].append(record.pos)
        for scaffold in coord_dict:
            coord_dict[scaffold] = np.array(coord_dict[scaffold])
        return coord_dict

    def get_positions(self):
        """
        Extracts coordinates of records in collection
        :return: dictionary of coordinates, keys are scaffold names, values - Numpy arrays of coordinates
        """
        positions_dict = OrderedDict({})
        for scaffold in self.records:
            positions_dict[scaffold] = np.array([[record.pos] for record in self.records[scaffold]])
        return positions_dict

    def check_by_ref_and_alt(self, ref_alt_list, flag, description="No description"):
        """

        :param ref_alt_list:
        :param flag:
        :return: None
        """
        self.metadata.add_metadata("##INFO=<ID=%s,Number=0,Type=Flag,Description=\"%s\">" % (flag, description))
        for record in self:
            record.check_ref_alt_list(ref_alt_list, flag)

    def filter_by_ref_and_alt(self, ref_alt_list):
        """

        :param ref_alt_list:
        :return: None
        """
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        return self.filter(lambda record: (record.ref, record.alt_list) in ref_alt_list)

    def set_filter(self, expression, filter_name):
        """
        Sets filter_name in FILTER field if expression returns True
        :param expression:
        :param filter_name:
        :return: None
        """
        for scaffold in self.records:
            for record in self.records[scaffold]:
                if expression(scaffold, record):
                    if "PASS" in record.filter_list or "." in record.filter_list:
                        record.filter_list = [filter_name]
                    else:
                        record.filter_list.append(filter_name)

    def set_filter_by_intersection_with_feature(self, annotation_dict, filter_name, mode="cross",
                                                feature_type_black_list=[]):
        """

        :param annotation_dict:
        :param filter_name:
        :param mode:
        :return:
        """
        if mode == "cross":
            def expression(scaffold, record):
                if scaffold not in annotation_dict:
                    return False
                return record.check_intersection_with_features(scaffold, annotation_dict,
                                                               feature_type_black_list=feature_type_black_list)

        elif mode == "no_cross":
            def expression(scaffold, record):
                if scaffold not in annotation_dict:
                    return True
                return not record.check_intersection_with_features(scaffold, annotation_dict,
                                                                   feature_type_black_list=feature_type_black_list)

        self.set_filter(expression, filter_name)

    def check_presence(self, chrom, position, alt_list=None):
        """
        Checks presence of variant in collection
        :param chrom:
        :param position:
        :param alt_list: optional
        :return: True if variant is present in collection(same chrom, position and optionally alts) otherwise False
        """

        if chrom not in self.scaffold_list:
            return False
        for record in self.records[chrom]:
            if record.pos > position:
                return False
            if record.pos == position:
                if alt_list:
                    if alt_list != record.alt_list:
                        return False
                return True

    def split_by_scaffolds(self):
        """

        :return:
        """
        return [CollectionVCF(metadata=self.metadata, records_dict={scaffold: self.records[scaffold]},
                              header=self.header, samples=self.samples, from_file=False)
                for scaffold in self.records]

    def get_location(self, annotation_dict, key="Loc", use_synonym=False, strand_key="strand",
                     synonym_dict=None, feature_type_black_list=[]):
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Locations of variant\">" % key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Strand\">" % strand_key)
        for scaffold in self.records:
            for record in self.records[scaffold]:
                record.get_location(scaffold, annotation_dict, key=key, use_synonym=use_synonym,
                                    synonym_dict=synonym_dict, feature_type_black_list=feature_type_black_list,
                                    strand_key=strand_key)

    @staticmethod
    def _reference(record):
        nucleotides = ["A", "C", "G", "T"]
        if record.ref in nucleotides:
            return record.ref
        return "INDEL"

    def count_heterozygous_snps(self, window_size, window_step, reference_scaffold_length_dict,
                                ignore_scaffolds_shorter_than_window=True, output_prefix=None,
                                skip_empty_windows=False, per_sample_output=False):

        def heterozygous_variant(record):
            #print record.__str__()
            #print not record.is_homozygous()
            return not record.is_homozygous()

        return self.count_variants_in_windows(window_size, window_step, reference_scaffold_length_dict,
                                              ignore_scaffolds_shorter_than_window=ignore_scaffolds_shorter_than_window,
                                              output_prefix=output_prefix,
                                              skip_empty_windows=skip_empty_windows,
                                              expression=heterozygous_variant, per_sample_output=per_sample_output)

    def draw_snps_histogram(self, window_size, window_step, output_prefix, reference_genome,
                            gaps_and_masked_positions_max_fraction=0.4,
                            expression=None, masking_gff=None, parsing_mode="parse", per_sample_output=False,
                            plot_type="concatenated",
                            xlabel="Position in genome",
                            ylabel="Number of SNPs",
                            title="SNP counts in windows",
                            suptitle="",
                            extensions=["png", ],
                            masked_or_gaped_region_mark=0,
                            figure_height_per_plot=3,
                            figure_width=12,
                            multiplier=1000):

        window_stepppp = window_size if window_step is None else window_step

        kb = multiplier / 1000
        mb = multiplier / 1000000

        normalized_ylabel = "%s per %i %s" % (ylabel, mb if mb >= 1 else kb if kb >=1 else multiplier, "Mbp" if mb >= 1 else "Kbp" if kb >= 1 else "bp")

        print("Parsing reference and...")
        reference = ReferenceGenome(reference_genome,
                                    masked_regions=None,
                                    index_file="refgen.idx",
                                    filetype="fasta",
                                    mode=parsing_mode,
                                    black_list=[],
                                    masking_gff_list=masking_gff)
        print("Merging gaps with masking...")

        gaps_and_masked_region_window_counts = reference.count_gaped_and_masked_positions_in_windows(window_size,
                                                                                                     window_stepppp,
                                                                                                     ignore_scaffolds_shorter_than_window=True,
                                                                                                     output_prefix=output_prefix,
                                                                                                     min_gap_len=1)
        print("Counting variants in windows...")
        variant_window_counts = self.count_variants_in_windows(window_size,
                                                               window_stepppp,
                                                               reference.region_length,
                                                               ignore_scaffolds_shorter_than_window=True,
                                                               output_prefix=output_prefix,
                                                               skip_empty_windows=False,
                                                               expression=expression,
                                                               per_sample_output=per_sample_output)

        normalized_variant_window_counts = SynDict()
        filtered_normalized_variant_window_counts = SynDict()

        # normalization
        if per_sample_output:
            for sample in variant_window_counts:
                normalized_variant_window_counts[sample] = SynDict()
                filtered_normalized_variant_window_counts[sample] = SynDict()
                for scaffold_id in variant_window_counts[sample]:
                        #print sample
                        #print scaffold_id
                        #print variant_window_counts[sample][scaffold_id]
                    normalized_variant_window_counts[sample][scaffold_id] = np.divide(variant_window_counts[sample][scaffold_id].astype(float), window_stepppp - gaps_and_masked_region_window_counts[scaffold_id] + 1) * multiplier
                        #print variant_window_counts[sample][scaffold_id]
                    filtered_normalized_variant_window_counts[sample][scaffold_id] = []
        else:
            for scaffold_id in variant_window_counts:
                normalized_variant_window_counts[scaffold_id] = np.divide(variant_window_counts[scaffold_id].astype(float), window_stepppp - gaps_and_masked_region_window_counts[scaffold_id] + 1) * multiplier
                filtered_normalized_variant_window_counts[scaffold_id] = []

        # filtering
        if per_sample_output:
            for sample in variant_window_counts:
                for scaffold_id in variant_window_counts[sample]:
                    for window_index in range(0, len(variant_window_counts[sample][scaffold_id])):
                        if np.isnan(variant_window_counts[sample][scaffold_id][window_index]):
                            normalized_variant_window_counts[sample][scaffold_id][window_index] = masked_or_gaped_region_mark
                        elif float(gaps_and_masked_region_window_counts[scaffold_id][window_index])/float(window_size) > gaps_and_masked_positions_max_fraction:
                            #print variant_window_counts.keys()
                            variant_window_counts[sample][scaffold_id][window_index] = masked_or_gaped_region_mark
                            normalized_variant_window_counts[sample][scaffold_id][window_index] = masked_or_gaped_region_mark
                        else:
                            filtered_normalized_variant_window_counts[sample][scaffold_id].append(normalized_variant_window_counts[sample][scaffold_id][window_index])
        else:
            for scaffold_id in variant_window_counts:
                for window_index in range(0, len(variant_window_counts[scaffold_id])):
                    if np.isnan(variant_window_counts[scaffold_id][window_index]):
                        normalized_variant_window_counts[scaffold_id][window_index] = masked_or_gaped_region_mark
                    elif float(gaps_and_masked_region_window_counts[scaffold_id][window_index])/float(window_size) > gaps_and_masked_positions_max_fraction:
                        variant_window_counts[scaffold_id][window_index] = masked_or_gaped_region_mark #variant_window_counts[scaffold_id]
                        normalized_variant_window_counts[scaffold_id][window_index] = masked_or_gaped_region_mark
                    else:
                        filtered_normalized_variant_window_counts[scaffold_id].append(normalized_variant_window_counts[scaffold_id][window_index])
        
        if per_sample_output:
            for sample in normalized_variant_window_counts:
                normalized_variant_window_counts.write("%s.%s.normalized_variant_number.tab" % (sample, output_prefix), splited_values=True)
                filtered_normalized_variant_window_counts.write("%s.%s.filtered.normalized_variant_number.tab" % (sample, output_prefix), splited_values=True)
                
        else:
            normalized_variant_window_counts.write("%s.normalized_variant_number.tab" % output_prefix, splited_values=True)
            filtered_normalized_variant_window_counts.write("%s.filtered.normalized_variant_number.tab" % output_prefix, splited_values=True)
            
        print("Drawing...")
        if plot_type == "concatenated":
            if per_sample_output:
                data = OrderedDict()
                normalized_data = OrderedDict()
                for sample in variant_window_counts:
                    data[sample] = []
                    normalized_data[sample] = []
                    for scaffold_id in reference.region_length:
                        if scaffold_id not in variant_window_counts[sample]:
                            continue
                        len(data[sample])
                        data[sample] += list(variant_window_counts[sample][scaffold_id]) + [0, ]
                        normalized_data[sample] += list(normalized_variant_window_counts[sample][scaffold_id]) + [0, ]
                #print data
                for sample in variant_window_counts:
                    data[sample] = np.array(data[sample])
                    normalized_data[sample] = np.array(normalized_data[sample])
                    bins = np.arange(len(data[sample]))
                    #print bins
                #print data[sample]

                sample_list = list(variant_window_counts.keys())
                sample_number = len(sample_list)

                figure, subplot_list = plt.subplots(nrows=sample_number, ncols=2, sharex=True, sharey=False, figsize=(figure_width, figure_height_per_plot * sample_number))
                for row_index in range(0, sample_number):
                    if row_index > 0:
                        subplot_list[row_index][0].get_shared_x_axes().join(subplot_list[row_index][0], subplot_list[0][0])
                        subplot_list[row_index][1].get_shared_x_axes().join(subplot_list[row_index][1], subplot_list[0][1])
                    for column_index in 0, 1:

                #for subplot_index in range(0, 2 * sample_number):
                    #if subplot_index % 2 == 0:
                        if column_index == 0:
                            subplot_list[row_index][column_index].plot(bins, data[sample_list[row_index]])
                            subplot_list[row_index][column_index].set_ylabel(ylabel)
                        else:
                            subplot_list[row_index][column_index].plot(bins, normalized_data[sample_list[row_index]])
                            subplot_list[row_index][column_index].set_ylabel(normalized_ylabel)
                        subplot_list[row_index][column_index].set_xlim(xmin=0)
                        subplot_list[row_index][column_index].set_xlabel(xlabel)
                        subplot_list[row_index][column_index].set_title(self.samples[row_index])
                plt.suptitle(suptitle)
            else:
                figure, subplot_list = plt.subplots(nrows=1, ncols=2,
                                                    sharex=True, sharey=False,
                                                    figsize=(figure_width, figure_height_per_plot ))
                data = []
                normalized_data = []
                for scaffold_id in reference.region_length:
                    if scaffold_id not in variant_window_counts:
                        continue
                    data += list(variant_window_counts[scaffold_id]) + [0, ]
                    normalized_data += list(normalized_variant_window_counts[scaffold_id]) + [0, ]
                data = np.array(data)
                normalized_data = np.array(normalized_data)
                #print normalized_data
                bins = np.arange(len(data)) #* window_step
                #print data
                #print max(data)
                #print bins
                for column_index in 0, 1:
                    if column_index == 0:
                        subplot_list[column_index].plot(bins, data)
                        subplot_list[column_index].set_ylabel(ylabel)
                    else:
                        subplot_list[column_index].plot(bins, normalized_data)
                        subplot_list[column_index].set_ylabel(normalized_ylabel)

                    subplot_list[column_index].set_xlim(xmin=0)
                    subplot_list[column_index].set_xlabel(xlabel)
                    subplot_list[column_index].set_title(title)
                plt.suptitle(suptitle)
        plt.tight_layout()
        for extension in extensions:
            plt.savefig("%s.%s" % (output_prefix, extension))

    @staticmethod
    def heterozygous_variant(record):
        #print record.__str__()
        #print not record.is_homozygous()
        return not record.is_homozygous()

    def heterozygous_sample_variant(self, record, sample_index):
        #print record.__str__()
        #print sample_index, self.samples[sample_index], not record.is_homozygous()
        return not record.is_homozygous_sample(sample_index)

    def draw_heterozygous_snps_histogram(self, window_size, window_step, output_prefix, reference_genome,
                                         gaps_and_masked_positions_max_fraction=0.4,
                                         masking_gff=None, parsing_mode="parse", per_sample_output=False,
                                         plot_type="concatenated",
                                         xlabel="Position in genome",
                                         ylabel="SNPs",
                                         title="SNP counts in windows",
                                         suptitle="",
                                         extensions=["png", ],
                                         masked_or_gaped_region_mark=0,
                                         figure_height_per_plot=3,
                                         figure_width=12,
                                         multiplier=1000):

        self.draw_snps_histogram(window_size, window_step, output_prefix, reference_genome,
                                 gaps_and_masked_positions_max_fraction=gaps_and_masked_positions_max_fraction,
                                 expression=self.heterozygous_sample_variant if per_sample_output else self.heterozygous_variant,
                                 masking_gff=masking_gff,
                                 parsing_mode=parsing_mode, per_sample_output=per_sample_output,
                                 plot_type=plot_type,
                                 xlabel=xlabel,
                                 ylabel=ylabel,
                                 title=title,
                                 suptitle=suptitle,
                                 extensions=extensions,
                                 masked_or_gaped_region_mark=masked_or_gaped_region_mark,
                                 figure_height_per_plot=figure_height_per_plot,
                                 figure_width=figure_width,
                                 multiplier=multiplier)

    def draw_variant_window_densities(self, reference_fasta, output_prefix, window_size, window_step,
                                      masking_gff=None,
                                      gap_fraction_threshold=0.4,
                                      parsing_mode="index_db", min_gap_length=10,
                                      masked_region_color="grey", gap_color="grey",
                                      no_snp_color="white",
                                      ignore_scaffolds_shorter_than_window=True,
                                      skip_empty_windows=False, reference_scaffold_black_list=(),
                                      figure_extensions=("png", "svg"),
                                      suptitle="Variant density",
                                      density_multiplicator=1000,
                                      scaffold_black_list=(),
                                      sort_scaffolds=False, scaffold_ordered_list=None,
                                      scaffold_white_list=[], add_sample_name_to_labels=False,
                                      figure_width=8,
                                      figure_height_scale_factor=0.5,
                                      sample_label="SampleZZZ",
                                      dist_between_scaffolds_scaling_factor=1,
                                      colormap_tuple_list=((0.0, "#333a97"), (0.1, "#3d3795"), (0.5, "#5d3393"),
                                                           (0.75, "#813193"), (1.0, "#9d2d7f"), (1.25, "#b82861"),
                                                           (1.5, "#d33845"), (2.0, "#ea2e2e"), (2.5, "#f5ae27")),
                                      colormap=None,
                                      thresholds=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),):
        if window_step > window_size:
            raise ValueError("ERROR!!! Window step can't be larger then window size")

        print("Parsing reference...")
        reference = ReferenceGenome(reference_fasta,
                                    masked_regions=None,
                                    index_file="refgen.idx",
                                    filetype="fasta",
                                    mode=parsing_mode,
                                    black_list=reference_scaffold_black_list,
                                    masking_gff_list=masking_gff)

        print("Merging gaps with masking...")

        gaps_and_masked_region_window_count_dict = reference.count_gaped_and_masked_positions_in_windows(window_size,
                                                                                                         window_step,
                                                                                                         ignore_scaffolds_shorter_than_window=True,
                                                                                                         output_prefix=output_prefix,
                                                                                                         min_gap_len=min_gap_length)

        count_dict = {sample_label: self.count_variants_in_windows(window_size, window_step, reference.region_length,
                                                                   ignore_scaffolds_shorter_than_window=ignore_scaffolds_shorter_than_window,
                                                                   output_prefix=output_prefix,
                                                                   skip_empty_windows=skip_empty_windows)}

        DrawingRoutines.draw_variant_window_densities(count_dict, reference.region_length, window_size, window_step,
                                                      output_prefix,
                                                      masking_dict=gaps_and_masked_region_window_count_dict,
                                                      gap_fraction_threshold=gap_fraction_threshold,
                                                      colormap_tuple_list=colormap_tuple_list,
                                                      record_style=None, ext_list=figure_extensions,
                                                      label_fontsize=13, left_offset=0.2,
                                                      figure_width=figure_width,
                                                      figure_height_scale_factor=figure_height_scale_factor,
                                                      scaffold_synonym_dict=None,
                                                      id_replacement_mode="partial",
                                                      suptitle=suptitle,
                                                      density_multiplicator=density_multiplicator,
                                                      scaffold_black_list=scaffold_black_list,
                                                      sort_scaffolds=sort_scaffolds,
                                                      scaffold_ordered_list=scaffold_ordered_list,
                                                      scaffold_white_list=scaffold_white_list,
                                                      gap_color=gap_color,
                                                      masked_color=masked_region_color,
                                                      no_snp_color=no_snp_color,
                                                      add_sample_name_to_labels=add_sample_name_to_labels,
                                                      dist_between_scaffolds_scaling_factor=dist_between_scaffolds_scaling_factor,
                                                      colormap=colormap,
                                                      thresholds=thresholds,)

        DrawingRoutines.draw_window_density_distribution(count_dict, output_prefix=output_prefix,
                                                         density_multiplicator=density_multiplicator,
                                                         suptitle="SNP density distribution",
                                                         number_of_bins=None, width_of_bins=None,
                                                         max_threshold=None, min_threshold=None,
                                                         scaffold_black_list=[], scaffold_white_list=[],
                                                         sort_scaffolds=False, scaffold_ordered_list=None, subplot_size=2,
                                                         per_scaffold_histo_dir="per_scaffold_histo_dir/",
                                                         subplot_tuple=None, share_x_axis=True, share_y_axis=True,
                                                         extensions=("png",))

    def hierarchical_clustering(self, method='average', dendrogramm_max_y=2000,
                                sample_name=None, save=False, clustering_dir="clustering",
                                dendrogramm_color_threshold=1000,
                                draw_dendrogramm=True,
                                write_inconsistent=True,
                                write_correlation=True):
        """
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
        allowed methods(used to calculate distance between clusters):
        'complete'    -   Farthest Point Algorithm
        'single'      -   Nearest Point Algorithm
        'average'     -   UPGMA algorithm, distance between clusters is calculated as average from pairwise
                           distances between elements of clusters
        'weighted     -   WPGMA algorithm
        'centroid'    -   UPGMC algorithm
        'median'      -   WPGMC algorithm
        'ward'        -   incremental algorithm
        :param method:
        :param dendrogramm_max_y:
        :param sample_name:
        :param save:
        :param clustering_dir:
        :param dendrogramm_color_threshold:
        :param draw_dendrogramm:
        :param write_inconsistent:
        :param write_correlation:
        :return:
        """
        positions_dict = self.get_positions()
        correlation_dict = OrderedDict({})
        linkage_dict = OrderedDict({})
        inconsistent_dict = OrderedDict({})
        clusters_dict = OrderedDict({})
        if draw_dendrogramm or write_correlation or write_inconsistent:
            os.system("mkdir -p %s" % clustering_dir)
        for region in positions_dict:
            #print positions_dict[region]
            if len(positions_dict[region]) <= 1:
                linkage_dict[region] = None
                correlation_dict[region] = None
                continue
            else:
                distance_matrix = pdist(positions_dict[region])
            #print(distance_matrix)

            linkage_dict[region] = linkage(distance_matrix, method=method)
            if draw_dendrogramm:
                plt.figure(1, dpi=150, figsize=(50, 20))
                dendrogram(linkage_dict[region],
                           color_threshold=dendrogramm_color_threshold,
                           leaf_font_size=4,
                           distance_sort=True)
                plt.ylim(ymax=dendrogramm_max_y)
                plt.axhline(y=500, color="purple")
                plt.axhline(y=1000, color="black")
                plt.savefig("%s/clustering_%s.svg" % (clustering_dir, region))
                plt.close()

            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.cophenet.html#scipy.cluster.hierarchy.cophenet
            # calculates cophenetic correlation coefficient to estimate accuracy of clustering
            correlation_dict[region] = cophenet(linkage_dict[region], distance_matrix)[0]

            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.inconsistent.html#scipy.cluster.hierarchy.inconsistent
            # calculates inconsistent coeff

            inconsistent_dict[region] = inconsistent(linkage_dict[region])
            if write_inconsistent:
                np.savetxt("%s/inconsistent_coefficient_%s.t" % (clustering_dir, region), inconsistent_dict[region])

            # clusters_dict[region] = fcluster(linkage_dict[region], 1)
            # np.savetxt("clustering/clusters_%s.t" % region, clusters_dict[region], fmt="%i")
        if write_correlation:
            sample = sample_name
            if not sample:
                sample = self.samples[0]
            with open("%s/correlation.t" % clustering_dir, "w") as cor_fd:
                cor_fd.write("sample\t%s\n" % ("\t".join(positions_dict.keys())))
                cor_fd.write("%s\t%s\n" % (sample,
                                           "\t".join([str(correlation_dict[region]) for region in positions_dict])))

        if save:
            self.linkage_dict = linkage_dict

        return linkage_dict

    def get_clusters(self,
                     extracting_method="inconsistent",
                     threshold=0.8,
                     cluster_distance='average',
                     dendrogramm_max_y=2000,
                     sample_name=None,
                     save_clustering=False,
                     clustering_dir="clustering",
                     dendrogramm_color_threshold=1000,
                     draw_dendrogramm=True,
                     return_collection=True,
                     write_inconsistent=True,
                     write_correlation=True):

        from MACE.Parsers.CCF import RecordCCF, CollectionCCF, MetadataCCF, HeaderCCF
        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                        dendrogramm_max_y=dendrogramm_max_y,
                                                        sample_name=sample_name,
                                                        save=save_clustering,
                                                        clustering_dir=clustering_dir,
                                                        dendrogramm_color_threshold=dendrogramm_color_threshold,
                                                        draw_dendrogramm=draw_dendrogramm,
                                                        write_correlation=write_correlation,
                                                        write_inconsistent=write_inconsistent)
        mut_clusters_dict = OrderedDict({})
        clusters = OrderedDict()
        for region in linkage_dict:
            if linkage_dict[region] is None:
                clusters[region] = None
            else:
                clusters[region] = fcluster(linkage_dict[region], threshold, criterion=extracting_method)

        if return_collection:
            record_ccf_dict = OrderedDict()
            for region in self.records:
                # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
                clusters_dict = OrderedDict({})
                if clusters[region] is None:
                    continue
                for i in range(0, len(clusters[region])):
                    if clusters[region][i] not in clusters_dict:
                        clusters_dict[clusters[region][i]] = [self.records[region][i]]
                    else:
                        clusters_dict[clusters[region][i]].append(self.records[region][i])



                record_ccf_dict[region] = [RecordCCF(collection_vcf=CollectionVCF(records_dict=dict([(region, clusters_dict[cluster])]), from_file=False),
                                                     from_records=True) for cluster in clusters_dict]

            return CollectionCCF(records_dict=record_ccf_dict, metadata=MetadataCCF(self.samples,
                                                                                    vcf_metadata=self.metadata,
                                                                                    vcf_header=self.header),
                                 header=HeaderCCF("CLUSTER_ID\tCHROM\tSTART\tEND\tDESCRIPTION".split("\t")))
        else:
            return clusters

    # methods for sum of two CollectionsVCF: no check for intersections(!!!!!!!!)
    def __add__(self, other):
        new_records_dict = deepcopy(self.records)
        for scaffold in other.scaffold_list:
            if scaffold in self.scaffold_list:
                new_records_dict[scaffold] += other.records[scaffold]
            else:
                new_records_dict[scaffold] = other.records[scaffold]
        return CollectionVCF(metadata=self.metadata, records_dict=new_records_dict,
                             header=self.header, samples=self.samples, from_file=False)

    def __radd__(self, other):
        return self.__add__(other)

    def test_thresholds(self,
                        extracting_method="inconsistent",
                        threshold=None,
                        cluster_distance='average',
                        dendrogramm_max_y=2000,
                        sample_name=None,
                        save_clustering=False,
                        testing_dir="threshold_test",
                        count_singletons=True,
                        scaffold_prefix="Region",
                        extensions=("svg", "png")):
        # TODO: adjust parameters of figure
        # threshold is tuple(list) of three variables: min, max, number

        # extracting_method possible values
        #   inconsistent
        #   distance
        #   maxclust
        #   monocrit
        #   monocrit

        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                        dendrogramm_max_y=dendrogramm_max_y,
                                                        sample_name=sample_name,
                                                        save=save_clustering,
                                                        clustering_dir=testing_dir)

        num_of_regions = len(list(linkage_dict.keys()))

        side = int(sqrt(num_of_regions))
        if side*side != num_of_regions:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=150, figsize=(30, 30))
        fig.suptitle("Relashionship between number of clusters and threshold of %s" % extracting_method, fontsize=20)

        thresholds = threshold
        if extracting_method == "inconsistent":
            if not threshold:
                thresholds = (0.5, 1.5, 21)

        index = 1
        for region in linkage_dict:
            if linkage_dict[region] is None:
                continue
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
            n_clusters_list = []
            n_nonsingleton_clusters = []
            n_multiclusters = []
            n_five_plus_clusters = []
            coef_threshold_list = np.linspace(*thresholds)  # best variant 0.5, 1.5, 21
            for i in coef_threshold_list:
                clusters = fcluster(linkage_dict[region], i, criterion=extracting_method)
                n_clusters_list.append(max(clusters))

                # counting clusters with 2+ and 3+ clusters
                ua, uind = np.unique(clusters, return_inverse=True)
                counted = np.bincount(uind)
                #j = 0
                nonsingleton = 0
                multicluster = 0  # 3+
                five_plus_clusters = 0 # 5+
                for k in counted:
                    if k > 1:
                        nonsingleton += 1
                    if k > 2:
                        multicluster += 1
                    if k > 4:
                        five_plus_clusters += 1
                n_nonsingleton_clusters.append(nonsingleton)
                n_multiclusters.append(multicluster)
                n_five_plus_clusters.append(five_plus_clusters)
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            #ax = plt.gca()
            #ax.set_xticks(np.arange(0.5, 2.2, 0.1))

            plt.grid()
            if count_singletons:
                plt.plot(coef_threshold_list, n_clusters_list, label="all")
            plt.plot(coef_threshold_list, n_nonsingleton_clusters, "green", label="2+")
            plt.plot(coef_threshold_list, n_multiclusters, "red", label="3+")
            plt.plot(coef_threshold_list, n_five_plus_clusters, "black", label="5+")
            plt.title("%s %s" % (scaffold_prefix, region))
            plt.legend(loc='upper right')
            plt.ylabel("Number of clusters")
            plt.xlabel("Threshold")
            #plt.axvline(x=0.8, color="purple")
            #plt.axvline(x=1.1, color="purple")

            plt.ylim(ymin=0)
            index += 1
        for ext in extensions:
            plt.savefig("%s/clusters_%s.%s" % (testing_dir, extracting_method, ext))
        plt.close()

    def add_info(self, metadata_line, expression, info_name, info_value=None):
        self.metadata.add_metadata(metadata_line)
        for record in self:
            if expression(record):
                value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
                if info_name in record.info_dict:
                    record.info_dict[info_name] += value
                else:
                    record.info_dict[info_name] = value

    def parse_snpeff_info_record(self, string, snpeff_entry="ANN"):
        if snpeff_entry == "EFF":
            effect, parameters = string.split("(")
            # remove closing bracket and split
            parameters = parameters[:-1].split("|")
            return [effect] + parameters

        elif snpeff_entry == "ANN":
            return string.split("|")

    def extract_snpeff_info(self, output_file, snpeff_entry="ANN"):

        snpeff_info_dict_keys = "EFF", "LOS", "NMD"
        record_header_list = ["Chrom", "Pos", "Ref", "Alt", "Filter"]
        if snpeff_entry == "EFF":
            snpeff_header_list = ["Effect", "Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_Change",
                                  "Amino_Acid_Length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                                  "Transcript_ID", "Exon_Rank", "Genotype_Number", "ERRORS", "WARNINGS"]
        elif snpeff_entry == "ANN":
            snpeff_header_list = ["Allele", "Annotation",
                                  "Putative_impact", "Gene_Name",
                                  "Gene_ID", "Feature type",
                                  "Feature ID", "Transcript biotype",
                                  "Rank", "HGVS.c",
                                  "HGVS.p", "cDNA_position",
                                  "CDS_position", "Protein_position",
                                  "Distance_to_feature", "Errors_Warnings"
                                  ]
        else:
            raise ValueError("ERROR!!! Unknow SNPeff entry: %s. Only ANN or EFF are allowed..." % snpeff_entry)

        #print(output_file)
        with open(output_file, "w") as out_fd:
            header_string = "#" + "\t".join(record_header_list + snpeff_header_list) + "\n"
            out_fd.write(header_string)
            for scaffold in self.records:
                for record in self.records[scaffold]:
                    common_part = "%s\t%i\t%s\t%s\t%s" % (scaffold, record.pos, record.ref, ",".join(record.alt_list),
                                                          ",".join(record.filter_list))

                    if snpeff_entry not in record.info_dict:
                        continue

                    for effect in record.info_dict[snpeff_entry]:
                        effect_parameters = self.parse_snpeff_info_record(effect, snpeff_entry)
                        num_parameters = len(effect_parameters)
                        for i in range(0, num_parameters):
                            if effect_parameters[i] == "":
                                effect_parameters[i] = "."
                        if num_parameters < 14:
                            effect_parameters += ["." for i in range(num_parameters, 14)]
                        out_fd.write(common_part + "\t" + "\t".join(effect_parameters) + "\n")

    def count_strandness(self, prefix):
        count_dict = OrderedDict({})

        for scaffold in self.records:
            count_dict[scaffold] = np.zeros((2, 4), dtype=int)
        hor_coord_dict = {"C": 0, "G": 1}
        ver_coord_dict = {"N": 0, "P": 1, "M": 2, "B": 3}

        for scaffold in self.records:
            for record in self.records[scaffold]:
                count_dict[scaffold][hor_coord_dict[record.ref]][ver_coord_dict[record.info_dict["Fstrand"][0]]] += 1

        count_dict["all"] = sum(count_dict.values())

        for chromosome in count_dict:
            with open("%s_%s.t" % (prefix, chromosome), "w") as out_fd:
                out_list = count_dict[chromosome].tolist()
                for index, name in zip(range(0, len(out_list)), ["C", "G"]):
                    out_list[index].insert(0, name)
                out_list.insert(0, [".", "N", "P", "M", "B"])
                for string_list in out_list:
                    out_fd.write("\t".join([str(x) for x in string_list]) + "\n")

        return count_dict

    def variants_start_end(self, left, right, record_dict, skip_genes_without_five_utr=False,
                           min_five_utr_len=10):

        gene_variants_positions = []
        all_variant_start_positions = []
        all_variant_end_positions = []

        for record_id in record_dict:
            for feature in record_dict[record_id].features:
                if feature.type != "gene":
                    continue
                if skip_genes_without_five_utr:
                    for sub_feature in feature.sub_features:
                        if sub_feature.type == "five_prime_UTR" and len(sub_feature) >= min_five_utr_len:
                            break
                    else:
                        continue
                #print(feature.sub_features)
                for sub_feature in feature.sub_features:
                    if sub_feature.type != "CDS":
                        continue
                    chrom = record_id
                    strand = sub_feature.strand
                    CDS_start = sub_feature.location.start + 1 if strand == +1 else sub_feature.location.end
                    CDS_end = sub_feature.location.end if strand == +1 else sub_feature.location.start + 1
                    #region_start = CDS_start - (args.left * strand)
                    #region_end = CDS_start + (args.right * strand)

                    region_start_start = CDS_start - left if strand == +1 else CDS_start - right
                    region_start_end = CDS_start + right if strand == +1 else CDS_start + left

                    region_end_start = CDS_end - left if strand == +1 else CDS_end - right
                    region_end_end = CDS_end + right if strand == +1 else CDS_end + left
                    #print("aaa")
                    start_coordinates = []
                    end_coordinates = []
                    for variant in self.records[record_id]:
                        if region_start_start <= variant.pos <= region_start_end:
                            start_coordinates.append((variant.pos - CDS_start) * strand)
                        if region_end_start <= variant.pos <= region_end_end:
                            end_coordinates.append((variant.pos - CDS_end) * strand)
                    all_variant_start_positions += start_coordinates
                    all_variant_end_positions += end_coordinates
                    #print(feature.qualifiers)
                    gene_variants_positions.append([feature.qualifiers["Name"], strand, chrom, region_start_start,
                                                    region_start_end, start_coordinates,
                                                    region_end_start, region_end_end,
                                                    end_coordinates])
        return all_variant_start_positions, all_variant_end_positions, gene_variants_positions

    def find_location(self, record_dict, key="Ftype", strand_key="Fstrand", genes_key="Genes", genes_strand_key="Gstrand",
                      feature_type_black_list=[],
                      use_synonym=False, synonym_dict=None, add_intergenic_label=True):

        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Types of features\">" % key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=1,Type=String,Description=\"Strand of features\">" % strand_key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Names of genes\">" % genes_key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Strands of genes\">" % genes_strand_key)
        for record in self:
            record.find_location(record_dict, key=key, strand_key=strand_key,
                                 genes_key=genes_key, genes_strand_key=genes_strand_key,
                                 feature_type_black_list=feature_type_black_list,
                                 use_synonym=use_synonym, synonym_dict=synonym_dict,
                                 add_intergenic_label=add_intergenic_label)

    def draw_info_distribution(self, info_dict_key, expression, outfile_prefix,
                               extension_list=(".svg", ".png"), bins=None,):
        scaffold_distribution = OrderedDict()
        for scaffold in self.scaffold_list:
            scaffold_distribution[scaffold] = [[], []]
            for record in self.records[scaffold]:
                #print(scaffold)
                if expression(record):
                    scaffold_distribution[scaffold][0] += record.info_dict[info_dict_key]
                else:
                    scaffold_distribution[scaffold][1] += record.info_dict[info_dict_key]

        #print(scaffold_distribution[scaffold][0])
        side = int(sqrt(self.number_of_scaffolds))
        if side*side != self.number_of_scaffolds:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=150, figsize=(15, 15))
        fig.suptitle("Distribution of %s" % info_dict_key, fontsize=20, fontweight='bold')

        index = 1
        for scaffold in self.scaffold_list:
            #print(scaffold)
            #print(scaffold_distribution[scaffold][0])
            #print(scaffold_distribution[scaffold][1])
            sub_plot_dict[scaffold] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            #ax = plt.gca()
            #ax.set_xticks(np.arange(0.5, 2.2, 0.1))

            plt.grid()
            num_of_bins = bins if bins is not None else 20
            maximum = max(max(scaffold_distribution[scaffold][0]) if scaffold_distribution[scaffold][0] else 0,
                          max(scaffold_distribution[scaffold][1]) if scaffold_distribution[scaffold][1] else 0)
            if isinstance(bins, Iterable):
                if maximum > num_of_bins[-1]:
                    num_of_bins[-1] = maximum + 1
            plt.hist([scaffold_distribution[scaffold][0], scaffold_distribution[scaffold][1]], bins=num_of_bins)
            plt.xlim(xmin=0)
            plt.title("%s" % scaffold, fontweight='bold')
            #plt.legend(loc='upper right')
            plt.ylabel("Number of variants")
            plt.xlabel("%s" % info_dict_key)
            #plt.axvline(x=0.8, color="purple")
            #plt.axvline(x=1.1, color="purple")

            plt.ylim(ymin=0)
            index += 1
        plt.subplots_adjust(hspace=0.27, wspace=0.27, top=0.92, left=0.05, right=0.99, bottom=0.04)
        for extension in extension_list:
            plt.savefig("%s%s" % (outfile_prefix, extension if extension[0] == "." else ".%s" % extension))
        plt.close()

    def set_filter_for_indels_in_homopolymers(self, reference_dict, min_homopolymer_len=4,
                                              filter_name="indel_in_homopolymer"):

        def expression(scaffold, record):
            if not record.check_indel():
                False
            if scaffold not in reference_dict:
                raise ValueError("Scaffold %s is absent in reference" % scaffold)

            scaffold_length = len(reference_dict[scaffold])

            reference_pos = record.pos - 1
            left_flank = None if record.pos == 1 else reference_dict[scaffold].seq[max(0, record.pos-20):record.pos]
            right_flank = None if (record.pos + 1) == scaffold_length else reference_dict[scaffold].seq[record.pos+1: max(scaffold_length - 1,
                                                                                                                          record.pos++20)]
            ref_var_len = len(record.ref)
            indel_type_list = [None for variant in record.alt_list]

            for i in range(0, len(record.alt_list)):
                alt_variant = record.alt_list[i]
                alt_var_len = len(alt_variant)
                if len(alt_variant) < ref_var_len:
                    if record.ref[0:alt_var_len] == alt_variant:
                        indel_type_list[i] = "right"
                    elif record.ref[-alt_var_len:] == alt_variant:
                        indel_type_list[i] = "left"
                elif len(alt_variant) > ref_var_len:
                    if alt_variant[0:ref_var_len] == record.ref:
                        indel_type_list[i] = "right"
                    elif alt_variant[-ref_var_len:] == record.ref:
                        indel_type_list[i] = "left"
                else:
                    continue
                homo_len = 1
                if indel_type_list[i] == "right":
                    for letter in range(1, len(right_flank)):
                        if letter != right_flank[0]:
                            break
                        homo_len += 1
                elif indel_type_list[i] == "left":
                    for letter in range(-2, -len(right_flank)-1, -1):
                        if letter != right_flank[0]:
                            break
                        homo_len += 1
                if homo_len >= min_homopolymer_len:
                    return True
            return False
        self.set_filter(expression, filter_name)


if __name__ == "__main__":
    pass
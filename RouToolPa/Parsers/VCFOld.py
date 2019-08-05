#!/usr/bin/env python
"""
VCF Parser Module
"""
__author__ = 'Sergei F. Kliver'

import os
import re
from math import sqrt
from copy import deepcopy
from collections import OrderedDict, Iterable

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from RouToolPa.Parsers.Abstract import Record, Collection, Metadata, Header
from RouToolPa.General.GeneralCollections import IdList, IdSet, SynDict, TwoLvlDict
from MACE.Routines import DrawingRoutines

ref_alt_variants = {"deaminases": [("C", ["T"]), ("G", ["A"])]
                    }


class RecordVCF(Record):
    """
    RecordVCF class
    """
    __slots__ = ["pos", "id", "ref", "alt_list", "qual", "filter_list", "info_dict", "samples_list", "flags"]

    def __init__(self, pos, id, ref, alt_list, qual, filter_list, info_dict, samples_list,
                 flags=None):
        """
        Initializes record
        :param pos: coordinate of mutation in chromosome
        :param id: id of mutation
        :param ref: reference variant
        :param alt_list: list of alternative variants
        :param qual: quality of mutation
        :param filter_list: list of filters
        :param info_dict: dictionary containing non flag data from vcf file
        :param samples_list: list of samples
        :param flags: flags from INFO field of vcf file
        :return: None
        """

        self.pos = pos                                  #int
        self.id = id                                    #str
        self.ref = ref                                  #str
        self.alt_list = alt_list                        #list, entries are strings
        self.qual = qual                                #real or "."
        self.filter_list = sorted(filter_list)          #list, entries are strings
        self.info_dict = info_dict                      #dict
        self.samples_list = samples_list                #list entries are dicts with keys from format_list and
                                                        #values are lists
        #self.description = description if description else {}
        self.flags = set(flags) if flags is not None else set([])

    def __str__(self):
        """

        :return: string representation of record (vcf string without chromosome name)
        """
        alt_string = ",".join(self.alt_list)
        filter_string = ";".join(self.filter_list)
        if self.info_dict:
            key_string_list = []
            for key in sorted(list(self.info_dict.keys())):
                if self.info_dict[key]:
                    key_string_list.append(key + "=" + ",". join(map(lambda x: str(x), self.info_dict[key])))
                else:
                    key_string_list.append(key)
            info_string = ";".join(key_string_list)
            info_string = "%s;%s" % (info_string, ";".join(self.flags)) if self.flags else info_string
        else:
            info_string = "."
        for sample in self.samples_list:
            if len(sample.keys()) > 1:
                format_string = ":".join(sample.keys())
                break
        else:
            format_string = list(self.samples_list[0].keys())[0]

        samples_string = "\t".join([":".join([",".join(map(lambda x: str(x), sample[key])) for key in sample.keys()]) for sample in self.samples_list])
        return '\t'.join(map(lambda x: str(x), [self.pos, self.id, self.ref, alt_string,
                                                self.qual, filter_string, info_string, format_string, samples_string]))

    def check_ref_alt_list(self, ref_alt_list, flag):
        """
        Sets flag in record if mutations is in list
        :param ref_alt_list: list of references and corresponding to them alternatives
        :param flag: flag to set if mutation is in ref_alt_list
        :return: None
        """
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        self.set_flag(lambda record: (record.ref, record.alt_list) in ref_alt_list, flag)

    def count_samples(self):
        """
        Counts samples with mutations
        :return: number of samples with mutation
        """
        #
        number = 0
        for sample in self.samples_list:
            if sample["GT"][0] != "./." and sample["GT"][0] != "0/0":
                number += 1
        return number

    def set_filter(self, expression, filter_name):
        """
        Adds filter in RecordVCF.filter_list if expression is True
        :param expression: expression to check
        :param filter_name: filter to set
        :return: None
        """
        if expression(self):
            self.filter_list.append(filter_name)
            self.filter_list.sort()

    def add_info(self, info_name, info_value=None):
        """
        Adds parameter to RecordVCF.info_dict
        :param info_name: name of parameter
        :param info_value: value of parameter
        :return: None
        """
        value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
        if info_name in self.info_dict:
            self.info_dict[info_name] += value
        else:
            self.info_dict[info_name] = value

    def check_indel(self):
        """
        Checks if record is indel
        :return: True if at least one variant of alternatives is indel
        """
        if len(self.ref) > 1 or len("".join(self.alt_list)) > len(self.alt_list):
            return True
        return False

    def find_location(self, scaffold, annotation_dict, key="Ftype", strand_key="Fstrand", genes_key="Genes",
                      genes_strand_key="Gstrand", feature_type_black_list=[],
                      use_synonym=False, synonym_dict=None, add_intergenic_label=True):
        """
        Finds location of mutations in annotations. Adds four parameters to RecordVCF.info_dict. By default their names are "Ftype", "Fstrand", "Genes", "Gstrand"
        "Ftype" contains list types of annotation within mutation is located,
        "Fstrand" - summary of strands(N for no annotation, P or M if annotation is located in plus or minus strand respectively and B if annotations from both strands are overlapped)
        "Genes" - list of gene names within mutation is located,
        "Gstrand" - list of gene strands
        :param scaffold: scaffold of variant
        :param annotation_dict: dictionary of Biopython SeqRecord objects (keys are record ids, i.e. names of chromosomes)
        :param key: key to use for annotation type
        :param strand_key: key to use for summary strand of annotations
        :param genes_key: key to use for genes list
        :param genes_strand_key: key to use for list of gene strands
        :param feature_type_black_list: list of annotation types to skip
        :param use_synonym: use or not synonyms for annotations
        :param synonym_dict: dictionary of synonyms
        :param add_intergenic_label: label to use if mutation is located not in any gene
        :return:
        """
        """
        This method is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        id of one SeqRecord in record_dict must be equal to record.pos
        locations will be written to description dictionary of record using "key" as key
        """
        if key not in self.info_dict:
            self.info_dict[key] = set([])

        # strandness values:
        # N - not defined
        # B - both
        # P - plus
        # M - minus
        strands = ["B", "P", "M"]
        if strand_key not in self.info_dict:
            self.info_dict[strand_key] = ["N"]
        for flag_key in (genes_key, genes_strand_key):
            if flag_key not in self.info_dict:
                self.info_dict[flag_key] = []

        for feature in annotation_dict[scaffold].features:
            if feature.type in feature_type_black_list:
                continue

            if (self.pos - 1) in feature:
                if feature.type == "gene" or feature.type == "ncRNA":
                    self.info_dict[genes_key].append(feature.qualifiers["Name"][0])
                    self.info_dict[genes_strand_key].append(strands[feature.strand])

                self.info_dict[key].add(self.get_synonym(feature.type, use_synonym=use_synonym,
                                                         synonym_dict=synonym_dict))
                if self.info_dict[strand_key][0] == "N":
                    self.info_dict[strand_key][0] = strands[feature.strand]
                elif strands[feature.strand] != self.info_dict[strand_key][0]:
                    self.info_dict[strand_key][0] = "B"
            else:
                continue

            for sub_feature in feature.sub_features:
                if sub_feature.type in feature_type_black_list:
                    continue
                if (self.pos - 1) in sub_feature:
                    self.info_dict[key].add(self.get_synonym(sub_feature.type, use_synonym=use_synonym,
                                                             synonym_dict=synonym_dict))
                    if self.info_dict[strand_key][0] == "N":
                        self.info_dict[strand_key][0] = strands[sub_feature.strand]
                    elif strands[sub_feature.strand] != self.info_dict[strand_key][0]:
                        self.info_dict[strand_key][0] = "B"

        if not self.info_dict[genes_key]:
            self.info_dict.pop(genes_key)
            self.info_dict.pop(genes_strand_key)

        if add_intergenic_label and (not self.info_dict[key]): # or ("gene" not in self.info_dict[key])):
            # igc == intergenic
            self.info_dict[key].add("igc")

    def is_homozygous(self):
        """
        Checks if variant in all samples is homozygous
        :return: True if variant in all samples is homozygous, otherwise False
        """
        for sample_dict in self.samples_list:
            zyg = sample_dict["GT"][0].split("/")
            if zyg[0] != zyg[1]:
                return False
        return True

    def is_homozygous_sample(self, sample_index):
        zyg = self.samples_list[sample_index]["GT"][0].split("/")
        if zyg[0] != zyg[1]:
            return False
        return True

    def is_homozygous_list(self):
        zyg_list = [sample_dict["GT"][0].split("/") for sample_dict in self.samples_list]
        return [None if (zyg[0] == ".") or (zyg[1] == ".") else True if zyg[0] == zyg[1] else False for zyg in zyg_list]

    def no_reference_allel_and_multiallel(self, sample_index=None, max_allels=None):
        if max_allels:
            if len(self.alt_list) > max_allels:
                return False

        if sample_index:
            if "0" in self.samples_list[sample_index]["GT"][0].split("/"):
                return False
        else:
            for sample in self.samples_list:
                if "0" in sample["GT"][0].split("/"):
                    return False
        return True


class MetadataVCF(OrderedDict, Metadata):
    """
    MetadataVCF class
    """
    def read(self, in_file):
        with open(in_file, "r") as fd:
            for line in fd:
                if line[:2] != "##":
                    # self.header = HeaderVCF(line[1:].strip().split("\t"))   # line[1:].strip().split("\t")
                    # self.samples = self.header[9:]
                    break
                self.metadata.add_metadata(line)


    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            #if "=" is not present in string (case of flag type in INFO field)
            return string, None
        return string[:index], string[index+1:]

    @staticmethod
    def _split_by_comma_sign(string):
        index_list = [-1]
        i = 1
        while (i < len(string)):
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
            if not isinstance(self[key], dict):
                metadata_string += "##%s=%s\n" % (key, self[key])
            else:
                prefix = "##%s=<" % key
                suffix = ">\n"
                for att_id in self[key]:
                    middle = "ID=%s," % att_id + ",".join(["%s=%s" % (param, self[key][att_id][param])
                                                           for param in self[key][att_id]])
                    metadata_string += prefix + middle + suffix
        return metadata_string[:-1]


class HeaderVCF(list, Header):
    """
    HeaderVCF class
    """

    def __str__(self):
        """
        :return: vcf-like string representation of header
        """
        return "#" + "\t".join(self)


class CollectionVCF(Collection):
    """
    CollectionVCF class

    """

    def __init__(self, metadata=None, records_dict=None, header=None, in_file=None, samples=None,
                 from_file=True, external_metadata=None, threads=1, dont_parse_info_and_data=False, parse_only_coordinates=False):
        """
        Initializes collection. If from_file is True collection will be read from file (arguments other then in_file, external_metadata and threads are ignored)
        Otherwise collection will be initialize from meta, records_dict, header, samples
        :param metadata:
        :param records_dict:
        :param header:
        :param in_file:
        :param samples:
        :param from_file:
        :param external_metadata:
        :param threads:
        :return:
        """
        self.linkage_dict = None
        if from_file:
            self.read(in_file, external_metadata=external_metadata,
                      dont_parse_info_and_data=dont_parse_info_and_data,
                      parse_only_coordinates=parse_only_coordinates)
        else:
            self.metadata = metadata
            self.records = {} if records_dict is None else records_dict
            self.header = header
            self.samples = samples
        self.scaffold_list = self.scaffolds()
        self.scaffold_length = self.scaffold_len()
        self.number_of_scaffolds = len(self.scaffold_list)
        self.record_index = self.rec_index()
        self.threads = threads

    def read(self, in_file, external_metadata=None, dont_parse_info_and_data=False, parse_only_coordinates=False):
        """
        Reads collection from vcf file
        :param in_file: path to file
        :param external_metadata: external(not from input file) metadata that could be used to parse records
        :return: None
        """
        self.metadata = MetadataVCF()
        self.records = OrderedDict({})
        with open(in_file, "r") as fd:
            for line in fd:
                if line[:2] != "##":
                    self.header = HeaderVCF(line[1:].strip().split("\t"))   # line[1:].strip().split("\t")
                    self.samples = self.header[9:]
                    break
                self.metadata.add_metadata(line)
            for line in fd:
                scaffold, record = self.add_record(line, external_metadata=external_metadata,
                                                   dont_parse_info_and_data=dont_parse_info_and_data,
                                                   parse_only_coordinates=parse_only_coordinates)
                if scaffold not in self.records:
                    self.records[scaffold] = []
                self.records[scaffold].append(record)

    def add_record(self, line, external_metadata=None, dont_parse_info_and_data=False, parse_only_coordinates=False):
        """
        Adds record to collection from line
        :param line: record line from vcf file
        :param external_metadata: external(not from input file) metadata that could be used to parse records
        :return: None
        """
        line_list = line.strip().split("\t")
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
        position = int(line_list[1])
        if parse_only_coordinates:
            return line_list[0], RecordVCF(position, ".", ".",[], ".", [], {}, [], flags=set())
        quality = "."

        if quality != line_list[5]:
            quality = float(line_list[5])
        alt_list = line_list[4].split(",")
        filter_list = line_list[6].split(",")          # list, entries are strings

        if dont_parse_info_and_data:
            return line_list[0], RecordVCF(position, line_list[2], line_list[3],
                                           alt_list, quality, filter_list,
                                           {}, [], flags=set())
        info_dict = OrderedDict()
        metadata = self.metadata if self.metadata else external_metadata
        flag_set = set([])

        if line_list[7] != ".":
            info_tuple_list = [self._split_by_equal_sign(entry) for entry in line_list[7].split(";")]
            #print line_list
            for entry in info_tuple_list:
                if entry[0] not in metadata["INFO"]:
                    # do not parse data from INFO field that are not described in metadata
                    continue
                if metadata["INFO"][entry[0]]["Type"] == "Flag":
                    flag_set.add(entry[0]) #info_dict[entry[0]] = []
                elif metadata["INFO"][entry[0]]["Type"] == "Integer":
                    info_dict[entry[0]] = list(map(lambda x: int(x), entry[1].split(",")))
                elif metadata["INFO"][entry[0]]["Type"] == "Float":
                    #print entry
                    info_dict[entry[0]] = list(map(lambda x: float(x), entry[1].split(",")))
                else:
                    info_dict[entry[0]] = entry[1].split(",")
        samples_list = []
        for sample_string in line_list[9:]:
            sample_dict = OrderedDict()
            if sample_string == "./.":
                sample_dict["GT"] = ["./."]
            else:
                for key, value_list in zip(line_list[8].split(":"), sample_string.split(":")):
                    if metadata["FORMAT"][key]["Type"] == "Integer":
                        #print key, value_list
                        sample_dict[key] = list(map(lambda x: x if x == "." else int(x), value_list.split(",")))
                    elif metadata["FORMAT"][key]["Type"] == "Float":
                        sample_dict[key] = list(map(lambda x: x if x == "." else float(x), value_list.split(",")))
                    else:
                        sample_dict[key] = value_list.split(",")

            samples_list.append(sample_dict)
        return line_list[0], RecordVCF(position, line_list[2], line_list[3],
                                      alt_list, quality, filter_list,
                                      info_dict, samples_list, flags=flag_set)

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
        # ignores sign in "
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

    def filter(self, expression):
        """
        Splits collection based on expression. Expression should be a function with one argument - record entry
        :param expression: filtering expression
        :return: tuple of two CollectionVCF. First contains records for which expression is True, second - False.
        """
        #
        filtered_records, filtered_out_records = self.filter_records(expression)
        return CollectionVCF(metadata=self.metadata, records_dict=filtered_records,
                             header=self.header, samples=self.samples, from_file=False),\
               CollectionVCF(metadata=self.metadata, records_dict=filtered_out_records,
                             header=self.header, samples=self.samples, from_file=False)

    def no_reference_allel_and_multiallel(self, record, sample_index=None, max_allels=None):
        return record.no_reference_allel_and_multiallel(sample_index=sample_index, max_allels=max_allels)

    def filter_variants_with_reference_allel_and_multiallelic(self, sample_index=None, max_allels=None):

        def expression(record):
            return self.no_reference_allel_and_multiallel(record, sample_index=sample_index, max_allels=max_allels)

        return self.filter(expression)

    def count_records(self, expression):
        """
        Counts records in collection based on expression. Expression should be a function with one argument - record entry
        :param expression: filtering expression
        :return: tuple of two numbers. First is number of records for which expression is True, second - False.
        """

        true_records = 0
        false_records = 0
        for scaffold in self.scaffold_list:
            for record in self.records[scaffold]:
                if expression(record):
                    true_records += 1
                else:
                    false_records += 1
        return true_records, false_records

    def count_zygoty(self):
        """
        :return: tuple - (N of homozygotes, N of heterozygotes)
        """
        return self.count_records(self.filter_zygoty_expression)

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

    def rainfall_plot(self, plot_name, base_colors=[], single_fig=True, dpi=300, figsize=(40, 40), facecolor="#D6D6D6",
                      ref_genome=None, masked_regions=None, min_gap_length=10, draw_gaps=False, suptitle=None,
                      gaps_color="#777777", masked_regions_color="#aaaaaa", logbase=2,
                      extension_list=("svg", "eps", "pdf", "png", "jpg"),
                      scaffold_black_list=None, scaffold_white_list=None,
                      scaffold_ordered_list=None, sort_scaffolds=False):
        """

        :param plot_name:
        :param base_colors:
        :param single_fig:
        :param dpi:
        :param figsize:
        :param facecolor:
        :param ref_genome:
        :param masked_regions:
        :param min_gap_length:
        :param draw_gaps:
        :param suptitle:
        :param gaps_color:
        :param masked_regions_color:
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
        reference_colors = {"A": "#FBFD2B",    # yellow
                            "C": "#FF000F",     # red
                            "G": "#000FFF",     # blue
                            "T": "#4ED53F",     # green
                            "INDEL": "#000000"  # black
                            }
        if base_colors:
            reference_colors = base_colors

        num_of_regions = self.number_of_scaffolds
        positions_dict = self.get_positions()
        distances_dict = {}
        region_reference_dict = {}
        os.system("mkdir -p %s" % plot_dir)
        if single_fig:
            fig = plt.figure(1, dpi=dpi, figsize=figsize, facecolor=facecolor)
            fig.suptitle(suptitle if suptitle else "Rainfall plot", fontsize=40, fontweight='bold', y=0.94)
            sub_plot_dict = OrderedDict({})
        index = 1

        final_scaffold_list = DrawingRoutines.get_filtered_scaffold_list(self.records,
                                                                         scaffold_black_list=scaffold_black_list,
                                                                         sort_scaffolds=sort_scaffolds,
                                                                         scaffold_ordered_list=scaffold_ordered_list,
                                                                         scaffold_white_list=scaffold_white_list,
                                                                         sample_level=False)

        for region in final_scaffold_list: # self.records

            # np.ediff1d return differences between consecutive elements in array, then 0 is added to the beginning
            distances_dict[region] = np.insert(np.ediff1d(positions_dict[region]), 0, 0)
            region_reference_dict[region] = OrderedDict({"A": [[], []],
                                                         "C": [[], []],
                                                         "G": [[], []],
                                                         "T": [[], []],
                                                         "INDEL": [[], []]})
            for i in range(0, len(self.records[region])):
                region_reference_dict[region][self._reference(self.records[region][i])][0].append(positions_dict[region][i])
                region_reference_dict[region][self._reference(self.records[region][i])][1].append(distances_dict[region][i])
            if single_fig:
                if not sub_plot_dict:
                    sub_plot_dict[region] = plt.subplot(num_of_regions, 1, index, axisbg=facecolor)
                else:
                    keys = list(sub_plot_dict.keys())
                    sub_plot_dict[region] = plt.subplot(num_of_regions, 1, index,
                                                        sharex=sub_plot_dict[keys[0]],
                                                        sharey=sub_plot_dict[keys[0]],
                                                        facecolor=facecolor)
                                                        #axisbg=facecolor)

                index += 1
                if draw_gaps:
                    if ref_genome:
                        for gap in ref_genome.gaps_dict[region]:
                            plt.gca().add_patch(plt.Rectangle((gap.location.start, 1),
                                                              gap.location.end - gap.location.start,
                                                              1024*32, facecolor=gaps_color, edgecolor='none'))
                # masked regions should be SeqRecord dict
                if masked_regions:
                    for feature in masked_regions[region].features:
                        plt.gca().add_patch(plt.Rectangle((int(feature.location.start)+1, 1),
                                                           feature.location.end - feature.location.start,
                                                           1024*32, facecolor=masked_regions_color, edgecolor='none'))

                for reference in region_reference_dict[region]:
                    plt.plot(region_reference_dict[region][reference][0],
                             region_reference_dict[region][reference][1],
                             color=reference_colors[reference],
                             marker='.', linestyle='None', label=reference)

                plt.text(-0.08, 0.5, region, rotation=0, fontweight="bold", transform=sub_plot_dict[region].transAxes,
                         fontsize=30,
                         horizontalalignment='center',
                         verticalalignment='center')
                plt.ylabel("Distanse")
                plt.axhline(y=100, color="#000000")
                plt.axhline(y=1000, color="#000000")
                plt.axhline(y=500, color="purple")
                plt.axhline(y=10, color="#000000")

        if single_fig:
            for region in sub_plot_dict:
                sub_plot_dict[region].set_yscale('log', basey=logbase)
            for extension in extension_list:
                plt.savefig("%s/%s_log_scale.%s" % (plot_dir, plot_name, extension))
            plt.close()

    def count_variants_in_windows(self, window_size, window_step, reference_scaffold_length_dict,
                                  ignore_scaffolds_shorter_than_window=True, output_prefix=None,
                                  skip_empty_windows=False, expression=None, per_sample_output=False):

        window_stepppp = window_size if window_step is None else window_step

        if window_stepppp > window_size:
            raise ValueError("ERROR!!! Window step can't be larger then window size")
        elif (window_size % window_stepppp) != 0:
            raise ValueError("ERROR!!! Window size is not a multiple of window step...")

        steps_in_window = window_size / window_stepppp

        short_scaffolds_ids = IdList()

        vcf_scaffolds = set(self.scaffold_list)
        reference_scaffolds = set(reference_scaffold_length_dict.keys())

        scaffolds_absent_in_reference = IdSet(vcf_scaffolds - reference_scaffolds)
        scaffolds_absent_in_vcf = IdSet(reference_scaffolds - vcf_scaffolds)
        if per_sample_output:
            count_dict = TwoLvlDict()
            for sample in self.samples:
                count_dict[sample] = SynDict()
        else:
            count_dict = SynDict()

        uncounted_tail_variants_number_dict = SynDict()

        if scaffolds_absent_in_reference:
            raise ValueError("ERROR!!! Some scaffolds from vcf file are absent in reference...")

        for scaffold_id in self.scaffold_list + list(scaffolds_absent_in_vcf):
            number_of_windows = self.count_number_of_windows(reference_scaffold_length_dict[scaffold_id],
                                                             window_size,
                                                             window_stepppp)
            if scaffold_id not in self.records:
                continue
            if number_of_windows == 0:
                short_scaffolds_ids.append(scaffold_id)
                if ignore_scaffolds_shorter_than_window:
                    continue

            if scaffold_id in scaffolds_absent_in_vcf:
                if skip_empty_windows:
                    continue
            if per_sample_output:
                for sample in self.samples:
                    #print scaffold_id
                    count_dict[sample][scaffold_id] = np.zeros(number_of_windows, dtype=np.int64)
                    #print count_dict[sample][scaffold_id]
            else:
                count_dict[scaffold_id] = np.zeros(number_of_windows, dtype=np.int64)

            uncounted_tail_variants_number_dict[scaffold_id] = 0

            variant_index = 0
            #print list(count_dict.keys())
            for variant in self.records[scaffold_id]:
                step_size_number = ((variant.pos - 1)/window_stepppp)

                if step_size_number - steps_in_window + 1 >= number_of_windows:
                    #print scaffold_id
                    #print self.scaffold_length[scaffold_id]
                    #print variant_index
                    #print("\n")
                    uncounted_tail_variants_number_dict[scaffold_id] = self.scaffold_length[scaffold_id] - variant_index
                    break
                if per_sample_output:

                    for i in range(max(step_size_number - steps_in_window + 1, 0),
                                    step_size_number + 1 if step_size_number < number_of_windows else number_of_windows):
                        #print step_size_number, steps_in_window
                        for sample_index in range(0, len(self.samples)):
                            sample_id = self.samples[sample_index]
                            #print sample_id, scaffold_id
                            if "GT" not in variant.samples_list[sample_index]:
                                print("WARNING: no genotype for sample %s for variant %s!!! Skipping..." % (sample_id, str(variant)))
                                continue
                            else:
                                if (variant.samples_list[sample_index]["GT"][0] == "0/0") or (variant.samples_list[sample_index]["GT"][0] == "./."):
                                    continue
                            if expression:
                                #print("AAAAAAAAAAAAAAAAAAAAAAAAA")
                                count_dict[sample_id][scaffold_id][i] += (1 if expression(variant, sample_index) else 0)
                                #if expression(variant, sample_index):
                                #print sample_id, sample_index, scaffold_id, i, count_dict[sample_id][scaffold_id][i]
                                #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBb"
                            else:
                                 count_dict[sample_id][scaffold_id][i] += 1

                else:
                    for i in range(max(step_size_number - steps_in_window + 1, 0),
                                   step_size_number + 1 if step_size_number < number_of_windows else number_of_windows):
                        if expression:
                            count_dict[scaffold_id][i] += 1 if expression(variant) else 0
                        else:
                            count_dict[scaffold_id][i] += 1

                variant_index += 1
        #print count_dict[self.samples[0]][list(count_dict[self.samples[0]].keys())[5]]
        #print "BBBBBBBBBBBBBB"
        #print count_dict[self.samples[0]]
        if output_prefix:
            scaffolds_absent_in_reference.write("%s.scaffolds_absent_in_reference.ids" % output_prefix)
            scaffolds_absent_in_vcf.write("%s.scaffolds_absent_in_vcf.ids" % output_prefix)
            uncounted_tail_variants_number_dict.write("%s.uncounted_tail_variant_number.tsv" % output_prefix)
            if per_sample_output:
                for sample in count_dict:
                    count_dict[sample].write("%s.%s.variant_number.tsv" % (output_prefix, sample), splited_values=True)
            else:
                count_dict.write("%s.variant_number.tsv" % output_prefix, splited_values=True)

        return count_dict

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

        print "Parsing reference and..."
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
                print normalized_data
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


class ReferenceGenome(object):
    """
    ReferenceGenome class
    """
    def __init__(self, ref_gen_file, masked_regions=None, index_file="refgen.idx", filetype="fasta", mode="index_db",
                 black_list=(), masking_gff_list=None, feature_mask_list=None, sort_scaffolds_by_length=True):
        """

        :param ref_gen_file: file with sequences of genome.
        :param masked_regions: dictionary with SeqRecords of masked regions.
        :param index_file: index file (created if absent) of genome used for parsing in index_db mode. Ignored  in other modes.
        :param filetype: format of file with sequences. Allowed formats - fasta, genbank.
        :param mode: mode of parsing reference genome. Allowed variants - to_dict, index, index_db.
        :param black_list: list of records to be not parsed.
        :return: None
        """
        self.ref_gen_file = ref_gen_file
        if mode == "index_db":
            self.reference_genome = SeqIO.index_db(index_file, [ref_gen_file], filetype)
        elif mode == "index":
            self.reference_genome = SeqIO.index(ref_gen_file, filetype)
        else:
            self.reference_genome = SeqIO.to_dict(SeqIO.parse(ref_gen_file, filetype))

        for record_id in self.reference_genome.keys():
            if record_id in black_list:
                self.reference_genome.pop(record_id, None)

        self.region_list = list(self.reference_genome.keys())
        self.length = 0
        self.feature_dict = OrderedDict()
        self.region_length = OrderedDict()

        for region in self.reference_genome:
            length = len(self.reference_genome[region])
            self.length += length
            self.region_length[region] = length

        self.length_to_region_dict = OrderedDict()
        for region in self.region_length:
            if self.region_length[region] in self.length_to_region_dict:
                self.length_to_region_dict[self.region_length[region]].append(region)
            else:
                self.length_to_region_dict[self.region_length[region]] = [region]
        #print self.length_to_region_dict
        lengths_list = sorted(list(self.length_to_region_dict.keys()), reverse=True)

        self.region_sorted_by_length_list = []
        #print self.length_to_region_dict
        for length in lengths_list:
            #print self.length_to_region_dict
            #type(self.region_sorted_by_length_list)
            #type(self.length_to_region_dict[length])
            self.region_sorted_by_length_list += self.length_to_region_dict[length]

        self.region_index = self.rec_index()
        self.gaps_dict = OrderedDict()
        self.masked_regions = masked_regions if masked_regions else OrderedDict()

        if masking_gff_list:
            masking_gffs = [masking_gff_list] if isinstance(masking_gff_list, str) else masking_gff_list
            for gff in masking_gffs:
                with open(gff, "r") as gff_fd:
                    for line in gff_fd:
                        if line[0] == "#":
                            continue
                        line_list = line.split("\t")
                        if feature_mask_list:
                            if line_list[2] not in feature_mask_list:
                                continue
                        if line_list[0] not in self.masked_regions:
                            self.masked_regions[line_list[0]] = []

                        self.masked_regions[line_list[0]].append(SeqFeature(FeatureLocation(int(line_list[3]) - 1,
                                                                                            int(line_list[4])),
                                                                            type="masked_region", strand=None))

    def merge_features_by_coordinate(self, feature_dict_list, return_seq_fetures_dict=True, feature_type='masked_region'):

        unified_dict = OrderedDict()
        merged_dict = OrderedDict()

        region_set = set()
        #print feature_dict_list

        for feature_dict in feature_dict_list:
            region_set |= set(feature_dict.keys())

        for region in region_set:
            unified_dict[region] = []
            merged_dict[region] = []

        for feature_dict in feature_dict_list:
            for region in feature_dict:
                for feature in feature_dict[region]:
                    unified_dict[region].append([feature.location.start, feature.location.end])

        for region in unified_dict:
            if unified_dict[region]:
                unified_dict[region].sort()
            if unified_dict[region] is None:
                print region

        for region in unified_dict:
            number_of_records = len(unified_dict[region])
            if number_of_records == 0:
                continue

            # [a, b) [c, d), a < b, c < d
            # after sorting c >= a
            i = 1

            prev_coordinates = deepcopy(unified_dict[region][0])

            while i < number_of_records:
                if unified_dict[region][i][0] > prev_coordinates[1]: # c > b
                    merged_dict[region].append(deepcopy(prev_coordinates))
                    prev_coordinates = deepcopy(unified_dict[region][i])
                elif unified_dict[region][i][1] > prev_coordinates[1]: # d > b; c<=b
                    prev_coordinates[1] = deepcopy(unified_dict[region][i][1])
                else: # d <= b
                    pass
                i += 1

            if unified_dict[region]:
                if prev_coordinates != unified_dict[region][-1]:
                    merged_dict[region].append(prev_coordinates)
            else:
                merged_dict[region].append(prev_coordinates)

        if return_seq_fetures_dict:
            feature_dict = OrderedDict()
            for region in merged_dict:
                feature_dict[region] = []
                for (start, stop) in merged_dict[region]:
                    feature_dict[region].append(SeqFeature(FeatureLocation(start, stop),
                                                           type=feature_type,
                                                           strand=None))

            return feature_dict
        else:
            return merged_dict

    def write_feature_dict_as_gff(self, feature_dict, output_gff):
        feature_id = 1
        with open(output_gff, "w") as out_gff:
            for region in feature_dict:
                for feature in feature_dict[region]:
                    out_gff.write("%s\tsource\t%s\t%i\t%i\t.\t%s\t.\t%s\n" % (region,
                                                                              feature.type,
                                                                              feature.location.start + 1,
                                                                              feature.location.end,
                                                                              "." if feature.location.strand == None else feature.location.strand,
                                                                              "ID=%i" % feature_id))
                    feature_id += 1

    def write_coords_dict_as_gff(self, coords_dict, output_gff, feature_type="masked_region"):
        feature_id = 1
        with open(output_gff, "w") as out_gff:
            for region in coords_dict:
                for coords in coords_dict[region]:
                    out_gff.write("%s\tsource\t%s\t%i\t%i\t.\t%s\t.\t%s\n" % (region,
                                                                              feature_type,
                                                                              coords[0] + 1,
                                                                              coords[1],
                                                                              ".",
                                                                              "ID=%i" % feature_id))
                    feature_id += 1

    def __len__(self):
        """

        :return: length of genome.
        """
        return self.length

    def rec_index(self):
        """
        Region order is based on descending region length
        :return:
        """
        index_dict = OrderedDict({})
        index_dict[self.region_sorted_by_length_list[0]] = [0, self.region_length[self.region_sorted_by_length_list[0]] - 1]
        for index in range(1, self.number_of_regions()):
            index_dict[self.region_sorted_by_length_list[index]] = [index_dict[self.region_sorted_by_length_list[index-1]][1] + 1,
                                            index_dict[self.region_sorted_by_length_list[index-1]][1] + self.region_length[self.region_sorted_by_length_list[index]]]
        return index_dict

    def get_position(self, coordinate):
        # coordinate should be 0-based
        tmp_coordinate = self.region_index[self.region_list[-1]][1] + coordinate + 1 if coordinate < 0 else coordinate
        if tmp_coordinate < 0:
            raise ValueError("Coordinate %i is too small for this genome" % tmp_coordinate)
        for region in self.region_list:
            start, end = self.region_index[region]
            if start <= tmp_coordinate <= end:
                coordinate_region = region
                break
        else:
            raise ValueError("Coordinate %i is too large for this genome" % tmp_coordinate)
        shift = tmp_coordinate - start
        return coordinate_region, shift

    def number_of_regions(self):
        """

        :return: number of regions/scaffolds/chromosomes in genome
        """
        return len(self.region_length)

    def find_gaps(self, min_gap_length=1):
        """
        Finds gaps (N) in reference genome and writes them as SeqFeatures to self.gaps_dict.
        Keys of dict are region names.
        :return: None
        """
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        for region in self.reference_genome:
            self.gaps_dict[region] = []
            gaps = gap_reg_exp.finditer(str(self.reference_genome[region].seq))  # iterator with
            for match in gaps:
                if (match.end() - match.start()) >= min_gap_length:
                    self.gaps_dict[region].append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                                             type="gap", strand=None))

    def generate_snp_set(self, size, substitution_dict=None, zygoty="homo", out_vcf="synthetic.vcf"):
        """
        Generates set of mutations
        :param size: size of snp set to be generated.
        :param substitution_dict: dictionary of substitutions.
        :param zygoty: zygoty of mutations in set.
        :param out_vcf: output .,vcf file with mutations
        """

        multiplier = (5 - len(substitution_dict)) if substitution_dict else 1
        unique_set = np.array([], dtype=np.int64)
        print("Generating...")

        index = 1
        while len(unique_set) < size:
            print("Iteration %i..." % index)
            print("    Generating raw set %i..." % index)
            raw_set = np.random.random_integers(0, high=self.length-1, size=size*multiplier)
            print("    Generated %i..." % len(raw_set))
            print("    Removing duplicates from raw set %i..." % index)
            raw_set = np.hstack((unique_set, raw_set))
            unique_set = np.unique(raw_set)
            np.random.shuffle(unique_set)
            np.random.shuffle(unique_set)
            print("    After removing  duplicates left %i.." % len(unique_set))
            print("    Checking filtered set %i..." % index)
            mutation_count = 0
            report_count = 1
            if substitution_dict is not None:
                tmp_list = []
                references = list(substitution_dict.keys())
                for coordinate in unique_set:
                    region, interregion_pos = self.get_position(coordinate)
                    if self.reference_genome[region][interregion_pos] in references:
                        if self.masked_regions and (region in self.masked_regions):
                            for masked_region in self.masked_regions[region].features:
                                #print(self.masked_regions[region].features)
                                #print(masked_region)
                                if interregion_pos in masked_region:
                                    break
                            else:
                                tmp_list.append(coordinate)
                                mutation_count += 1
                        else:
                            tmp_list.append(coordinate)
                            mutation_count += 1
                    if mutation_count/report_count == 500:
                        print("    Generated %i" % mutation_count)
                        report_count += 1
                unique_set = np.array(tmp_list)
            index += 1
        print("    Generated %i" % mutation_count)
        nucleotides = {"A", "C", "G", "T"}
        alt_dict = {}
        if not substitution_dict:
            for nuc in nucleotides:
                alt_dict[nuc] = list(nucleotides - {nuc})
        else:
            for ref in substitution_dict:
                if substitution_dict[ref]:
                    alt_dict[ref] = list(set(substitution_dict[ref]) - {ref})
                else:
                    alt_dict[ref] = list(nucleotides - {ref})
        print("Writing vcf...")

        coordinates_dict = dict((record_id, []) for record_id in self.region_list)
        check = 0
        for coordinate in unique_set[:size]:
            region, interregion_pos = self.get_position(coordinate)
            coordinates_dict[region].append(interregion_pos)

        with open(out_vcf, "w") as out_fd:
            out_fd.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            out_fd.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSynthetic_uniform_sample\n")
            #mutation_count = 0
            for region in coordinates_dict:
                if not coordinates_dict[region]:
                    continue
                num_of_snps = len(coordinates_dict[region])
                check += num_of_snps

                coordinates_dict[region].sort()
                #region, interregion_pos = self.get_position(coordinate)
                for interregion_pos in coordinates_dict[region]:
                    ref = self.reference_genome[region][interregion_pos]

                    if len(substitution_dict[ref]) != 1:
                        alt = alt_dict[ref][np.random.randint(0, len(alt_dict[ref]))]
                    else:
                        alt = alt_dict[ref][0]
                    if zygoty == "homo":
                        zyg = "1/1"
                    elif zygoty == "hetero":
                        zyg = "0/1"
                    out_fd.write("%s\t%i\t.\t%s\t%s\t.\t.\t.\tGT\t%s\n" %
                                 (region, interregion_pos + 1, ref, alt, zyg))
                print("    %s : %i SNPs" % (region, num_of_snps))
                    #mutation_count += 1
                    #if mutation_count == size:
                    #    break
            print ("Totaly %i mutations were written" % check)

    @staticmethod
    def count_number_of_windows(scaffold_length, window_size, window_step):
        if scaffold_length < window_size:
            return 0
        return int((scaffold_length - window_size)/window_step) + 1

    def count_gaped_and_masked_positions_in_windows(self, window_size, window_step,
                                                    ignore_scaffolds_shorter_than_window=True,
                                                    output_prefix=None,
                                                    min_gap_len=1):

        window_stepppp = window_size if window_step is None else window_step
        if window_stepppp > window_size:
            raise ValueError("ERROR!!! Window step can't be larger then window size")
        elif (window_size % window_stepppp) != 0:
            raise ValueError("ERROR!!! Window size is not a multiple of window step...")

        if not self.gaps_dict:
            self.find_gaps(min_gap_len)

        steps_in_window = window_size / window_stepppp

        short_scaffolds_ids = IdList()

        count_dict = SynDict()

        uncounted_tail_variants_number_dict = SynDict()

        masked_region_dict = self.merge_features_by_coordinate([self.masked_regions,
                                                                self.gaps_dict],
                                                               return_seq_fetures_dict=False)

        self.write_coords_dict_as_gff(masked_region_dict, "%s.merged_masked_regions.gff" % output_prefix,
                                      feature_type="masked_region")
        self.write_feature_dict_as_gff(self.gaps_dict, "%s.gaps.gff" % output_prefix)

        for scaffold_id in self.region_length:
            number_of_windows = self.count_number_of_windows(self.region_length[scaffold_id],
                                                             window_size,
                                                             window_stepppp)
            if number_of_windows == 0:
                short_scaffolds_ids.append(scaffold_id)
                if ignore_scaffolds_shorter_than_window:
                    continue

            scaffold_windows_list = np.zeros(number_of_windows, dtype=np.int64) #[0 for i in range(0, number_of_windows)]

            uncounted_tail_variants_number_dict[scaffold_id] = 0

            for masked_region_location in masked_region_dict[scaffold_id]:
                max_start_step = masked_region_location[0] / window_stepppp
                min_start_step = max(max_start_step - steps_in_window + 1, 0)
                max_end_step = (masked_region_location[1] - 1) / window_stepppp
                min_end_step = max(max_end_step - steps_in_window + 1, 0)

                if min_start_step >= number_of_windows:
                    break

                for i in range(min_start_step, min(max_start_step + 1, min_end_step, number_of_windows)):
                    scaffold_windows_list[i] += (i * window_step) + window_size - masked_region_location[0]

                for i in range(min_end_step, min(max_start_step + 1, number_of_windows)):
                    scaffold_windows_list[i] += masked_region_location[1] - masked_region_location[0]

                for i in range(max_start_step + 1, min(min_end_step, number_of_windows)):
                    scaffold_windows_list[i] += window_size

                for i in range(max(max_start_step + 1, min_end_step), min(max_end_step + 1, number_of_windows)):
                    scaffold_windows_list[i] += masked_region_location[1] - (i * window_stepppp)

            count_dict[scaffold_id] = scaffold_windows_list

        if output_prefix:
            count_dict.write("%s.gapped_and_masked_site_counts.tsv" % output_prefix, splited_values=True)

        return count_dict


if __name__ == "__main__":
    pass
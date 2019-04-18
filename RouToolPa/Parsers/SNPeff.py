#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import re
from math import sqrt
from copy import deepcopy
from collections import OrderedDict, Iterable

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from MACE.Parsers.Abstract import Record, Collection, Metadata, Header


class RecordSNPeff(Record):
    def __init__(self, pos, ref, alt_list, filter_list, effect, effect_impact, functional_class, codon_change,
                 amino_acid_change, amino_acid_len, gene_name, transcript_biotype, gene_coding, transcript_id, exon_rank,
                 genotype_number, errors, warnings, gene_name_alias_list=[], gene_function_list=[], gene_description=None,
                 info_dict=None, flags=None, biochemical_pathway=[]):

        Record.__init__(self, pos=pos, info_dict=info_dict, flags=flags)

        self.ref = ref
        self.alt_list = alt_list
        self.filter_list = filter_list
        self.effect = effect
        self.effect_impact = effect_impact
        self.functional_class = functional_class
        self.codon_change = codon_change
        self.amino_acid_change = amino_acid_change
        self.amino_acid_len = amino_acid_len
        self.gene_name = gene_name
        self.gene_name_alias_list = gene_name_alias_list
        self.gene_function_list = gene_function_list
        self.gene_description = gene_description
        self.transcript_biotype = transcript_biotype
        self.gene_coding = gene_coding
        self.transcript_id = transcript_id
        self.exon_rank = exon_rank
        self.genotype_number = genotype_number
        self.errors = errors
        self.warnings = warnings
        self.biochemical_pathway = biochemical_pathway

    def __str__(self):
        vcf_part = "%i\t%s\t%s\t%s" % (self.pos, self.ref, ",".join(self.alt_list) if self.alt_list else ".",
                                       ",".join(self.filter_list) if self.filter_list else ".")
        snpeff_part = ("\t" + "%s\t" * 17 + "%s") % (self.effect,
                                                     self.effect_impact,
                                                     self.functional_class,
                                                     self.codon_change,
                                                     self.amino_acid_change,
                                                     str(self.amino_acid_len) if self.amino_acid_len else ".",
                                                     self.gene_name,
                                                     self.transcript_biotype,
                                                     self.gene_coding,
                                                     self.transcript_id,
                                                     self.exon_rank,
                                                     self.genotype_number,
                                                     self.errors,
                                                     self.warnings,
                                                     ",".join(self.gene_name_alias_list) if self.gene_name_alias_list else ".",
                                                     ",".join(self.gene_function_list) if self.gene_function_list else ".",
                                                     ",".join(self.gene_description) if self.gene_description else ".",
                                                     ";".join(self.biochemical_pathway) if self.biochemical_pathway else ".")

        return vcf_part + snpeff_part


class HeaderSNPeff(list, Header):
    """
    HeaderVCF class
    """

    def __str__(self):
        """
        :return: tab-separatd representation of header
        """
        return "#" + "\t".join(self)


class CollectionSNPeff(Collection):

    def __init__(self, metadata=None, records_dict=None, header=None, input_file=None, from_file=False, filetype='tab'):
        # metadata should be Metadata class
        # header should be Header class
        # record_list should be list of Record class
        if from_file:
            self.read(input_file, filetype=filetype)
        else:
            self.records = OrderedDict() if records_dict is None else records_dict
            self.header = header
            self.metadata = metadata
        self.scaffold_list = self.scaffolds()
        self.scaffold_length = self.scaffold_len()
        self.number_of_scaffolds = len(self.scaffold_list)
        self.record_index = self.rec_index()

    def read(self, in_file, filetype='tab'):
        """
        Reads collection from vcf file
        :param in_file: path to file
        :param external_metadata: external(not from input file) metadata that could be used to parse records
        :return: None
        """
        self.metadata = None
        self.records = OrderedDict({})
        with open(in_file, "r") as fd:
            if filetype == 'tab':
                self.header = HeaderSNPeff(fd.readline()[1:].strip().split("\t"))
                gene_alias_presence = "Gene_Name_Aliases" in self.header
                gene_function_presence = "Gene_Functions" in self.header
                gene_description_presence = "Gene_Description" in self.header
                biochemical_pathway_presence = "Biochemical_Pathway" in self.header
                for line in fd:
                    scaffold, record = self.add_record_from_tab_line(line,
                                                                     gene_alias_presence=gene_alias_presence,
                                                                     gene_function_presence=gene_function_presence,
                                                                     gene_description_presence=gene_description_presence,
                                                                     biochemical_pathway_presence=biochemical_pathway_presence)
                    if scaffold not in self.records:
                        self.records[scaffold] = []
                    self.records[scaffold].append(record)

            elif filetype == 'vcf':
                raise ValueError("Parsing from vcf file annotated by SNPeff is not supported yet")

            else:
                raise ValueError("Wrong filetype for parsing")

        if not gene_alias_presence:
            self.header.append("Gene_Name_Aliases")
        if not gene_function_presence:
            self.header.append("Gene_Functions")
        if not gene_description_presence:
            self.header.append("Gene_Description")
        if not biochemical_pathway_presence:
            self.header.append("Biochemical_Pathway")

    def add_record_from_tab_line(self, line,
                                 gene_alias_presence,
                                 gene_function_presence,
                                 gene_description_presence,
                                 biochemical_pathway_presence):
        """
        Adds record to collection from line
        :param line: record line from tab file
        :return: None
        """
        line_list = line.strip().split("\t")
        pos = int(line_list[1])
        ref = line_list[2]
        alt_list = line_list[3].split(",")
        filter_list = line_list[4].split(",")
        effect = line_list[5]
        effect_impact = line_list[6]
        functional_class = line_list[7]
        codon_change = line_list[8]
        amino_acid_change = line_list[9]
        amino_acid_len = int(line_list[10]) if line_list[10] != "." else None
        gene_name = line_list[11]
        transcript_biotype = line_list[12]
        gene_coding = line_list[13]
        transcript_id = line_list[14]
        exon_rank = line_list[15]
        genotype_number = line_list[16]
        errors = line_list[17]
        warnings = line_list[18]
        gene_name_alias_list = line_list[self.header.index("Gene_Name_Aliases")].split(",") if gene_alias_presence else []
        gene_function_list = line_list[self.header.index("Gene_Functions")].split(",") if gene_function_presence else []
        gene_description = line_list[self.header.index("Gene_Description")].split(",") if gene_description_presence else []
        biochemical_pathway = line_list[self.header.index("Biochemical_Pathway")].split(",") if biochemical_pathway_presence else []

        return line_list[0], RecordSNPeff(pos, ref, alt_list, filter_list, effect, effect_impact,
                                          functional_class, codon_change,
                                          amino_acid_change, amino_acid_len, gene_name, transcript_biotype,
                                          gene_coding, transcript_id, exon_rank,
                                          genotype_number, errors, warnings,
                                          gene_name_alias_list=gene_name_alias_list,
                                          gene_function_list=gene_function_list,
                                          gene_description=gene_description,
                                          info_dict=None, flags=None,
                                          biochemical_pathway=biochemical_pathway)

    def add_gene_name_aliases(self, alias_dict):
        for scaffold in self.records:
            for record in self.records[scaffold]:
                if record.gene_name in alias_dict:
                    tmp_aliases = [alias_dict[record.gene_name]] if isinstance(alias_dict[record.gene_name], str) else alias_dict[record.gene_name]
                    new_aliases = []
                    for alias in tmp_aliases:
                        if alias != ".":
                            new_aliases.append(alias)
                    record.gene_name_alias_list = (record.gene_name_alias_list + new_aliases) if record.gene_name_alias_list else new_aliases

    def add_gene_functions(self, function_dict):
        for scaffold in self.records:
            for record in self.records[scaffold]:
                for gene_name in [record.gene_name] + record.gene_name_alias_list:
                    if gene_name in function_dict:
                        new_functions = [function_dict[gene_name]] if isinstance(function_dict[gene_name], str) else function_dict[gene_name]
                        record.gene_function_list = (record.gene_function_list + new_functions) if record.gene_function_list else new_functions
                        break

    def add_gene_description(self, description_dict):
        for scaffold in self.records:
            for record in self.records[scaffold]:
                for gene_name in [record.gene_name] + record.gene_name_alias_list:
                    if gene_name in description_dict:
                        new_descriptions = [description_dict[gene_name]] if isinstance(description_dict[gene_name], str) else description_dict[gene_name]
                        record.gene_description = (record.gene_description + new_descriptions) if record.gene_description else new_descriptions
                        break

    def add_biochemical_pathway(self, biochemical_pathway_dict):
        for scaffold in self.records:
            for record in self.records[scaffold]:
                for gene_name in [record.gene_name] + record.gene_name_alias_list:
                    if gene_name in biochemical_pathway_dict:
                        new_pathway = [biochemical_pathway_dict[gene_name]] if isinstance(biochemical_pathway_dict[gene_name], str) else biochemical_pathway_dict[gene_name]
                        record.biochemical_pathway = (record.biochemical_pathway + new_pathway) if record.biochemical_pathway else new_pathway
                        break

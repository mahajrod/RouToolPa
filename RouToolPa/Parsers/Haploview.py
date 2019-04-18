#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from collections import Iterable, OrderedDict
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from RouToolPa.Parsers.Abstract import Header, Collection


class HaploviewHeader(list, Header):

    def __str__(self):
        """
        :return: vcf-like string representation of header
        """
        return "#" + "\t".join(self)


class HaploviewCollection(Collection):

    def __init__(self, metadata=None, records_dict=None, header=None, input_file=None, from_file=False):
        # metadata should be Metadata class
        # header should be Header class
        # record_list should be list of Record class
        if from_file:
            self.read(input_file)
        else:
            self.records = OrderedDict() if records_dict is None else records_dict
            self.metadata = metadata
            self.header = header
        self.scaffold_list = self.scaffolds()
        self.scaffold_length = self.scaffold_len()
        self.number_of_scaffolds = len(self.scaffold_list)
        self.record_index = self.rec_index()

    def __str__(self):
        pass

    def __len__(self):
        return self.record_index[self.scaffold_list[-1]][1] + 1

    def read(self, input_file, parsing_mode="full"):
        """
        :param input_file:
        :param parsing_mode: full - parse all fields, short - parse coordinates, D' and r2
        :return:
        """
        with open(input_file, "r") as in_fd:
            if parsing_mode == "full":
                self.header = HaploviewHeader(in_fd.readline().strip().split("\t"))
            else:
                self.header = HaploviewHeader(["L1", "L2", "D'", "r^2"])
            current_chr = None
            for line in in_fd:
                tmp = line.strip().split()
                chr, first_snp = tmp[0].split(":")
                if (current_chr is None) or (chr != current_chr):
                    self.records[chr] = []
                    if current_chr is not None:
                        self.records[current_chr] = np.array(self.records[current_chr])
                    current_chr = chr
                second_snp = tmp[1].split(":")[1]
                if parsing_mode == "full":
                    tmp = [int(first_snp), int(second_snp)] + map(float, tmp[2:])
                    self.parsing_mode = "full"
                elif parsing_mode == "short":
                    tmp = [int(first_snp), int(second_snp), float(tmp[2]), float(tmp[4])]
                    self.parsing_mode = "short"
                self.records[chr].append(tmp)

    def filter(self, expression):
        filtered_dict = OrderedDict()
        filtered_out_dict = OrderedDict()

        for scaffold in self.records:
            filtered_dict[scaffold] = []
            for record in self.records[scaffold]:
                if expression(record):
                    filtered_dict[scaffold].append(record)
                else:
                    filtered_out_dict[scaffold].append(record)

            if filtered_dict[scaffold]:
                filtered_dict[scaffold] = np.array(filtered_dict[scaffold])
            else:
                filtered_dict.pop(scaffold, None)

            if filtered_out_dict[scaffold]:
                filtered_out_dict[scaffold] = np.array(filtered_out_dict[scaffold])
            else:
                filtered_out_dict.pop(scaffold, None)

        return HaploviewCollection(records_dict=filtered_dict), HaploviewCollection(records_dict=filtered_out_dict)

    @staticmethod
    def filter_record_by_snp_distance(record, max_dist):
        if (record[1] - record[0]) <= max_dist:
            return True
        return False

    def filter_by_snp_distance(self, max_dist):

        def expression(record):
            self.filter_record_by_snp_distance(record, max_dist)

        return self.filter(expression)

    def write_as_tsv(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#Scaffold\t" + "\t".join(self.header) + "\n")
            for scaffold in self.records:
                for record in self.records[scaffold]:
                    out_fd.write("%s\t%s\n" % (scaffold, "\t".join(map(str, record))))

"""
    def filter_records(self, expression):
        # expression should be a function with one argument - record entry
        filtered_records = OrderedDict({})
        filtered_out_records = OrderedDict({})
        for scaffold in self.scaffold_list:
            filtered_records[scaffold] = []
            filtered_out_records[scaffold] = []
            for record in self.records[scaffold]:
                if expression(record):
                    filtered_records[scaffold].append(record)
                else:
                    filtered_out_records[scaffold].append(record)
        return filtered_records, filtered_out_records

    def write(self, output_file, desired_scaffolds=None):
        with open(output_file, "w") as out_fd:
            if self.metadata:
                out_fd.write(str(self.metadata) + "\n")
            if self.header:
                out_fd.write(str(self.header) + "\n")
            scaffolds = desired_scaffolds if desired_scaffolds else self.scaffold_list
            for scaffold in scaffolds:
                for record in self.records[scaffold]:
                    out_fd.write(scaffold + "\t" + str(record) + "\n")

    def get_location(self, annotation_dict, key="Loc", use_synonym=False, strand_key="strand",
                     synonym_dict=None, feature_type_black_list=[]):
        for scaffold in self.records:
            for record in self.records[scaffold]:
                record.get_location(scaffold, annotation_dict, key=key, use_synonym=use_synonym,
                                    synonym_dict=synonym_dict, feature_type_black_list=feature_type_black_list,
                                    strand_key=strand_key)

    def set_location_flag(self, annotation_dict, expression, flag):
        for scaffold in self.records:
            for record in self.records[scaffold]:
                record.set_location_flag(scaffold, annotation_dict, expression, flag)

    def split_records_by_flags(self, flag_set, mode="all"):
        # possible modes:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list
        flags = set(flag_set)
        return self.filter_records(lambda record: flags & record.flags == flags if mode == "all" else flags & record.flags)

    def split_by_flags(self, flag_set, mode="all"):
        # possible modes:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list
        flags = set(flag_set)
        return self.filter(lambda record: flags & record.flags == flags if mode == "all" else flags & record.flags)

    def count_locations(self, annotation_black_list=[],
                        allow_several_counts_of_record=False,
                        out_filename="location_counts.t",
                        write=True,
                        count_dir="location_counts"):
        os.system("mkdir -p %s" % count_dir)
        region_counts_dict = TwoLvlDict({})
        for region in self.records:
            count_locations_dict = {"igc": 0, "unknown": 0}

            for record in self.records[region]:
                if "Loc" not in record.info_dict:
                    count_locations_dict["unknown"] += 1
                    continue
                if not record.info_dict["Loc"]:
                    count_locations_dict["unknown"] += 1
                    continue
                if allow_several_counts_of_record:
                    for location in record.info_dict["Loc"]:
                        if location in annotation_black_list:
                            continue
                        if location not in count_locations_dict:
                            count_locations_dict[location] = 1
                        else:
                            count_locations_dict[location] += 1
                else:
                    full_location = []
                    for location in record.info_dict["Loc"]:
                        if location in annotation_black_list:
                            continue
                        full_location.append(location)
                    if not full_location:
                        continue
                    full_location.sort()
                    full_location = "/".join(full_location)
                    if full_location not in count_locations_dict:
                        count_locations_dict[full_location] = 1
                    else:
                        count_locations_dict[full_location] += 1

            labels = []
            counts = []
            for location in count_locations_dict:
                if count_locations_dict[location] == 0 or location in annotation_black_list:
                    continue
                labels.append(location)
                counts.append(count_locations_dict[location])
            region_counts_dict[region] = OrderedDict([(label, count) for label, count in zip(labels, counts)])

        if write:
            region_counts_dict.write("%s/%s" % (count_dir, out_filename))
        return region_counts_dict

    def location_pie(self, pie_name="Location of variants", annotation_colors=[],
                     dpi=150, figsize=(30, 30),
                     explode=True, annotation_black_list=[],
                     allow_several_counts_of_record=False,
                     pie_prefix="variant_location_pie1",
                     full_genome_pie_prefix="variant_location_pie_full_genome",
                     counts_filename="location_counts.t",
                     plot_dir="variant_location_pie",
                     counts_dir="location_counts",
                     radius = 1,
                     draw_percents_on_single_pie=False,
                     combine_mixed=False,
                     extension_list=["svg", "eps", "pdf", "png", "jpg"]):
        print("Drawing location pie...")

        os.system("mkdir -p %s" % plot_dir)
        reference_colors = {"CDS": "#FBFD2B",    # yellow
                            "5'_UTR": "#FF000F",
                            "five_prime_UTR": "#FF000F",     # red
                            "3'_UTR": "#000FFF",
                            "three_prime_UTR": "#000FFF",     # blue
                            "igc": "#4ED53F",     # green
                            "ncRNA": 'cyan',
                            "other": "#ADB696",
                            "tRNA": "magenta"
                            }

        if annotation_colors:
            reference_colors = annotation_colors
            reference_colors["other"] = "#ADB696"

        count_locations_dict = self.count_locations(annotation_black_list=annotation_black_list,
                                                    allow_several_counts_of_record=allow_several_counts_of_record,
                                                    out_filename=counts_filename,
                                                    write=True,
                                                    count_dir=counts_dir)
        num_of_regions = len(count_locations_dict)
        side = int(sqrt(num_of_regions))
        if side*side != num_of_regions:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=dpi, figsize=figsize)
        fig.suptitle(pie_name, fontsize=20)

        index = 1
        all_labels = []
        all_counts = []
        all_colors = []

        for region in count_locations_dict:
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            labels = []
            counts = []
            colors = []
            if combine_mixed:
                labels.append("other")
                counts.append(0)
                colors.append(reference_colors["other"])
            for label in count_locations_dict[region]:
                if count_locations_dict[region][label] == 0 or label in annotation_black_list:
                    continue
                if combine_mixed and "/" in label:
                    counts[0] += count_locations_dict[region][label]
                    continue

                labels.append(label)
                counts.append(count_locations_dict[region][label])
                if label not in reference_colors:
                    colors.append(reference_colors["other"])
                else:
                    colors.append(reference_colors[label])

            explodes = np.zeros(len(counts))
            if explode and counts:
                max_count_index = counts.index(max(counts))
                explodes[max_count_index] = 0.1
            plt.pie(counts, explode=explodes, labels=labels, colors=colors,
                    autopct='%1.1f%%', shadow=True, startangle=90)
            plt.title("Region %s" % region, y=1.08)
            # Set aspect ratio to be equal so that pie is drawn as a circle.
            plt.axis('equal')
            index += 1
            for i in range(0, len(labels)):
                if labels[i] not in all_labels:
                    all_labels.append(labels[i])
                    all_counts.append(counts[i])
                    all_colors.append(colors[i])
                else:
                    label_index = all_labels.index(labels[i])
                    all_counts[label_index] += counts[i]

        for extension in extension_list:
            plt.savefig("%s/%s.%s" % (plot_dir, pie_prefix, extension))
        plt.close()

        fig = plt.figure(3, dpi=dpi, figsize=(4, 4))
        fig.suptitle(pie_name, fontsize=20)
        plt.subplot(1, 1, 1, axisbg="#D6D6D6")
        all_explodes = np.zeros(len(all_counts))
        if explode and all_counts:
            max_count_index = all_counts.index(max(all_counts))
            all_explodes[max_count_index] = 0.1
        if len(all_labels) > 0:
            max_label_length = max([len(x) for x in all_labels])
            max_count = max(all_counts)
            max_letters = 1
            while int(max_count / 10**max_letters) != 0:
                max_letters += 1
        if draw_percents_on_single_pie:
            patches, texts, autotexts = plt.pie(all_counts, explode=all_explodes, colors=all_colors,
                                                autopct='%1.1f%%', shadow=True, startangle=90, radius=2)
        else:
            patches, texts = plt.pie(all_counts, explode=all_explodes, colors=all_colors,
                                     shadow=True, startangle=90, radius=2)
        percent = 100 * np.array(all_counts).astype(np.float32, copy=False)/sum(all_counts)

        labels = ['{0}  -  {1} ({2:1.2f}%)'.format(i.ljust(max_label_length), str(j).ljust(max_letters), k) for i, j, k in zip(all_labels, all_counts, percent)]

        if len(all_labels) > 0:
            sort_legend = True
            if sort_legend:
                patches, labels, dummy = zip(*sorted(zip(patches, labels, all_counts),
                                                     key=lambda x: x[2],
                                                     reverse=True))

        plt.legend(patches, labels, loc='center left',  fontsize=13, bbox_to_anchor=(-1.0, 0.5))
        # plt.title("Full genome")
        plt.axis('equal')
        for extension in extension_list:
            plt.savefig("%s/%s.%s" % (plot_dir, full_genome_pie_prefix, extension), bbox_inches='tight')
        plt.close()
"""
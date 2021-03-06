#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import datetime
from collections import OrderedDict
import numpy as np
import pandas as pd
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines import MathRoutines




class GenomeCov(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bedtools genomecov", path=path, max_threads=max_threads)

    @staticmethod
    def parse_options(input_bam=None, report_zero_coverage=None,
                      one_based_coordinates=None, scale=None, bedgraph_output=False,
                      genome_bed=None):

        options = " -ibam %s" % input_bam
        if bedgraph_output and report_zero_coverage:
            options += " -bga"
        elif bedgraph_output:
            options += " -bg"
        else:
            options += " -d" if one_based_coordinates else " -dz"
        options += " -scale %f" % scale if scale else ""
        options += " -g %s" % genome_bed if genome_bed else ""

        return options

    def get_coverage(self, output, input_bam=None, report_zero_coverage=None,
                     one_based_coordinates=None, scale=None, bedgraph_output=False,
                     genome_bed=None):

        options = self.parse_options(input_bam=input_bam, report_zero_coverage=report_zero_coverage,
                                     one_based_coordinates=one_based_coordinates, scale=scale,
                                     bedgraph_output=bedgraph_output, genome_bed=genome_bed)
        options += " > %s" % output

        self.execute(options)

    def get_coverage_for_gff(self, input_file, genome_bed, output=None):

        options = " -i %s" % input_file
        options += " -g %s" % genome_bed
        options += " > %s" % output if output else ""

        self.execute(options=options)

    def get_bam_coverage_stats(self, input_bam, output_prefix, genome_bed=None,
                               max_coverage=None, min_coverage=None,
                               verbose=True, calc_stats=True):
        #options_list = []

        each_position_coverage_file = "%s.tab" % output_prefix
        coverage_stat_file = "%s.stat" % output_prefix

        each_position_options = self.parse_options(input_bam=input_bam, report_zero_coverage=True,
                                                   one_based_coordinates=True, genome_bed=genome_bed)
        each_position_options += " | gzip > %s.gz " % each_position_coverage_file

        region_options = self.parse_options(input_bam=input_bam, report_zero_coverage=True,
                                            bedgraph_output=True, genome_bed=genome_bed)

        region_options += "  | gzip > %s.bedgraph.gz " % output_prefix

        options_list = [each_position_options, region_options]

        self.parallel_execute(options_list, threads=2)

        #~/Soft/MAVR/scripts/math/get_stats_from_numer_data.py -i A_ventralis_pe.nodup.q20.tab -o A_ventralis_pe.nodup.q20.stat -l 2
        if calc_stats:
            MathRoutines.get_stats_from_file(each_position_coverage_file, minimum=min_coverage, maximum=max_coverage,
                                             dtype=float, comments="#", delimiter="\t", converters=None, skiprows=0,
                                             usecols=2, unpack=False, ndmin=0, output_file=coverage_stat_file,
                                             verbose=verbose)

    @staticmethod
    def mean_from_dict(coverage_dict):
        sum = 0
        total_sites = 0
        for coverage in coverage_dict:
            sum += coverage * coverage_dict[coverage]
            total_sites += coverage_dict[coverage]

        return float(sum) / float(total_sites)

    @staticmethod
    def median_from_dict(coverage_dict):
        total_sites = 0
        for coverage in coverage_dict:
            total_sites += coverage_dict[coverage]
        if total_sites % 2 == 0:
            left_half_sites = total_sites / 2
            right_half_sites = left_half_sites + 1
        else:
            left_half_sites = right_half_sites = int(total_sites / 2) + 1

        sorted_coverage = list(coverage_dict.keys())
        sorted_coverage.sort()

        count = 0
        for i in range(0, len(sorted_coverage)):
            count += coverage_dict[sorted_coverage[i]]

            if count >= right_half_sites:
                #print ("1111")
                return float(sorted_coverage[i])
            elif count == left_half_sites:
                #print("11111")
                return float(sorted_coverage[i] + sorted_coverage[i+1]) / 2

    def get_stats_from_coverage_file_stream_version(self, coverage_file, output_prefix, verbose=True,
                                                    scaffold_column=0, coverage_column=1,
                                                    separator="\t", buffering=None):

        stats = OrderedDict()
        summary_stats = OrderedDict()
        with self.metaopen(coverage_file, "r", buffering=buffering) as in_fd:
            line_list = in_fd.readline().strip().split(separator)
            scaffold, coverage = line_list[scaffold_column], int(line_list[coverage_column])
            coverage_dict = OrderedDict([(coverage, 1)])
            summary_coverage_dict = OrderedDict([(coverage, 1)])
            current_scaffold = scaffold
            line_counter = 1
            for line in in_fd:
                line_list = line.strip().split(separator)
                scaffold, coverage = line_list[scaffold_column], int(line_list[coverage_column])
                if coverage in summary_coverage_dict:
                    summary_coverage_dict[coverage] += 1
                else:
                    summary_coverage_dict[coverage] = 1
                line_counter += 1
                if line_counter % 1000000 == 0:
                    print("%s\tProcessed %i lines" % (str(datetime.datetime.now()), line_counter))
                if scaffold != current_scaffold:
                    #print(scaffold)
                    print("%s\tCalculating stats for %s" % (str(datetime.datetime.now()), current_scaffold))
                    stats[current_scaffold] = [sum(list(coverage_dict.values())),
                                               min(list(coverage_dict.keys())),
                                               max(list(coverage_dict.keys())),
                                               self.mean_from_dict(coverage_dict),
                                               self.median_from_dict(coverage_dict)]
                    coverage_dict = OrderedDict([(coverage, 1)])
                    current_scaffold = scaffold

                else:
                    if coverage in coverage_dict:
                        coverage_dict[coverage] += 1
                    else:
                        coverage_dict[coverage] = 1
            else:
                #print("END")
                #print(scaffold)
                stats[current_scaffold] = [sum(list(coverage_dict.values())),
                                           min(list(coverage_dict.keys())),
                                           max(list(coverage_dict.keys())),
                                           self.mean_from_dict(coverage_dict),
                                           self.median_from_dict(coverage_dict)]

        summary_stats["all"] = [sum(list(summary_coverage_dict.values())),
                                min(list(summary_coverage_dict.keys())),
                                max(list(summary_coverage_dict.keys())),
                                self.mean_from_dict(summary_coverage_dict),
                                self.median_from_dict(summary_coverage_dict)]

        #print(stats)
        stats = pd.DataFrame.from_dict(stats, orient="index", columns=["length", "min", "max", "mean", "median"])
        summary_stats = pd.DataFrame.from_dict(summary_stats, orient="index", columns=["length", "min", "max", "mean", "median"])
        stats.to_csv("%s.per_scaffold.stat" % output_prefix, sep="\t", index_label="#scaffold")
        summary_stats.to_csv("%s.all.stat" % output_prefix, sep="\t", index_label="#scaffold")
        if verbose:
            print(stats)

    def get_coverage_stats_in_windows(self, coverage_file, window_size, output, window_step=None,
                                      buffering=None):
        win_step = window_size if window_step is None else window_step
        stats = []

        with self.metaopen(coverage_file, "r", buffering=buffering) as in_fd:
            prev_scaffold, position, coverage = in_fd.readline().strip().split("\t")
            coverage_list = [int(coverage)]
            for line in in_fd:
                current_scaffold, position, coverage = line.strip().split("\t")
                coverage = int(coverage)
                if current_scaffold == prev_scaffold:
                    coverage_list.append(coverage)
                else:
                    scaffold_length = len(coverage_list)
                    if scaffold_length >= window_size:
                        number_of_windows = int((scaffold_length - window_size) / win_step) + 1
                        for i in range(0, number_of_windows):
                            win_start = i * win_step
                            window_coverage_list = coverage_list[win_start:win_start + window_size]
                            uncovered = window_coverage_list.count(0)
                            stats.append([prev_scaffold,
                                          scaffold_length,
                                          i,
                                          np.mean(window_coverage_list),
                                          np.median(window_coverage_list),
                                          np.min(window_coverage_list),
                                          np.max(window_coverage_list),
                                          uncovered,
                                          float(uncovered) / float(window_size)], )

                    prev_scaffold = current_scaffold
                    coverage_list = [coverage]

            scaffold_length = len(coverage_list)
            if scaffold_length >= window_size:
                number_of_windows = int((scaffold_length - window_size) / win_step) + 1
                for i in range(0, number_of_windows):
                    win_start = i * win_step
                    window_coverage_list = coverage_list[win_start:win_start + window_size]
                    uncovered = window_coverage_list.count(0)
                    stats.append([prev_scaffold,
                                  scaffold_length,
                                  i,
                                  np.mean(window_coverage_list),
                                  np.median(window_coverage_list),
                                  np.min(window_coverage_list),
                                  np.max(window_coverage_list),
                                  uncovered,
                                  float(uncovered)/float(window_size)],)

        stats = pd.DataFrame.from_records(stats, index=("scaffold", "window"),
                                          columns=("scaffold", "scaffold_length", "window", "mean",
                                                   "median", "min", "max", "uncovered", "uncovered,fraction"))
        stats.to_csv(output, sep="\t", header=True, index=True)

    @staticmethod
    def get_coverage_stats(coverage_file, output, verbose=True):
        print("Reading...")
        coverage_df = pd.read_csv(coverage_file, sep='\t', header=None, index_col=(0, 1), 
                                  names=("scaffold", "position", "coverage"))

        print("Calculating stats...")
        stat_df = coverage_df.groupby(level=0).agg(["min", "max", "mean", "median"])

        stat_df.to_csv(output, sep="\t")
        if verbose:
            print(stat_df)
        
        return stat_df
    
    def collapse_coverage_file(self, coverage_file, output_file):

        awk_string = "awk -F'\\t' 'BEGIN {SCAF=\"\"; LEN=\"\"; COV=\"\"} {if (($1 != SCAF)) {if (NR > 1) {printf \"%%s\\t%%s\\t%%s\\n\",SCAF,LEN, COV}; SCAF=$1; LEN=$2; COV=$3} else {LEN=$2; COV=COV\",\"$3}} ; END {printf \"%%s\t%%s\t%%s\n\",SCAF,LEN, COV} %s > %s" % (coverage_file, output_file)

        self.execute(options="", cmd=awk_string)

    def convert_per_base_coverage_file_to_vcf(self, coverage_file, output_vcf, sample_name, reference_file,
                                              parsing_mode="parse", format="fasta"):
        vcf_header = "##fileformat=VCFv4.2\n"
        vcf_header += '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">\n'
        vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name

        reference_dict = self.parse_seq_file(reference_file, parsing_mode, format=format)

        with self.metaopen(output_vcf, "w") as out_fd:
            out_fd.write(vcf_header)

            for line_list in self.file_line_as_list_generator(coverage_file):
                out_fd.write("%s\t%s\t.\t%s\t.\t.\t.\t.\tDP\t%s\n" % (line_list[0],
                                                                      line_list[1],
                                                                      str(reference_dict[line_list[0]].seq[int(line_list[1])-1]),
                                                                      line_list[2]))

    @staticmethod
    def extract_data_for_cds_from_collapsed_coverage_file(collapsed_coverage_file, cds_bed_file, output_file,
                                                          skip_transcript_with_no_cds=False):
        cds_dict = OrderedDict()

        with open(cds_bed_file, "r") as cds_fd:
            for line in cds_fd:
                line_list = line.strip().split("\t")
                # convert coordinates to python notation
                cds_dict[line_list[0]] = (int(line_list[1]) - 1, int(line_list[2]))
        with open(collapsed_coverage_file, "r") as col_cov_fd:
            with open(output_file, "w") as out_fd:
                for line in col_cov_fd:
                    transcript_id, length, coverage_array = line.strip().split("\t")

                    if transcript_id not in cds_dict:
                        print("No CDS for transcript %s. %s" % (transcript_id, "Skipping..." if skip_transcript_with_no_cds else ""))
                        if skip_transcript_with_no_cds:
                            continue
                        out_fd.write(line)
                        continue
                    coverage_array = coverage_array.split(",")
                    cds_len = cds_dict[transcript_id][1] - cds_dict[transcript_id][0]
                    out_fd.write("%s\t%s\t%s\n" % (transcript_id, str(cds_len),
                                                   ",".join(coverage_array[cds_dict[transcript_id][0]:cds_dict[transcript_id][1]])))

    @staticmethod
    def analyze_collapsed_coverage_file(collapsed_file, output_file):
        line_number = 0

        def int_through_float(string):
            return int(float(string))

        with open(collapsed_file, "r") as in_fd:
            with open(output_file, "w") as out_fd:

                output_header = "#Scaffold\tLength\tMean_coverage\tMedian_coverage\tMin_coverage\tMax_coverage\tCoveraged_positions\tZero_coverage_position_number\tZero_coveraged_region_number\tLeading_zero_covarage_len\tTrailing_zero_covarage_len\tZero_coverage_region_coordinates\tMean_coverage(without_zerocoveraged_ends)\tMedian_coverage(without_zerocoveraged_ends)\tLength_without_zero_coverage_ends\n"

                out_fd.write(output_header)
                for line in in_fd:
                    line_number += 1
                    #print (line_number)
                    #print [line]
                    if line == "\n" or line == "":    # skip blank lines
                        continue
                    tmp = line.strip().split("\t")
                    #print tmp
                    record_id = tmp[0]
                    record_len = int(tmp[1])

                    try:
                        coverage_array = map(int, tmp[2].split(","))
                    except ValueError:
                        coverage_array = map(int_through_float, tmp[2].split(","))

                    if record_len != len(coverage_array):
                        raise ValueError("Malformed line %i" % line_number)

                    if sum(coverage_array) == 0:
                        continue

                    mean_coverage = float(np.mean(coverage_array))
                    median_coverage = float(np.median(coverage_array))
                    min_coverage = float(np.min(coverage_array))
                    max_coverage = float(np.max(coverage_array))
                    coveraged_position_number = np.count_nonzero(coverage_array)

                    zero_coverage_position_number, zero_coverage_regions_list = MathRoutines.find_flat_regions_in_array(coverage_array, value=0)

                    if zero_coverage_position_number > 0:

                        if zero_coverage_regions_list[0][0] == 0:
                            leading_zero_coverage_len = zero_coverage_regions_list[0][1]
                            start_coverage_coordinate = zero_coverage_regions_list[0][1]
                        else:
                            leading_zero_coverage_len = 0
                            start_coverage_coordinate = 0

                        if zero_coverage_regions_list[-1][0] + zero_coverage_regions_list[-1][1] == record_len:
                            trailing_zero_coverage_len = zero_coverage_regions_list[-1][1]
                            end_coverage_coordinate = zero_coverage_regions_list[-1][0]
                        else:
                            trailing_zero_coverage_len = 0
                            end_coverage_coordinate = record_len

                        zero_coveraged_region_number = len(zero_coverage_regions_list)
                        zero_coverage_coordinates_list = []
                        for (start, length) in zero_coverage_regions_list:
                            zero_coverage_coordinates_list.append("%i-%i" % (start + 1, start + length))

                    else:
                        leading_zero_coverage_len = 0
                        trailing_zero_coverage_len = 0
                        zero_coveraged_region_number = 0
                        start_coverage_coordinate = 0
                        end_coverage_coordinate = record_len
                        zero_coverage_coordinates_list = ["."]

                    length_without_zerocoveraged_ends = end_coverage_coordinate - start_coverage_coordinate
                    mean_coverage_without_zero_coverage_ends = float(np.mean(coverage_array[start_coverage_coordinate:end_coverage_coordinate]))
                    median_coverage_without_zero_coverage_ends = float(np.median(coverage_array[start_coverage_coordinate:end_coverage_coordinate]))

                    out_fd.write("%s\t%i\t%.2f\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\t%.2f\t%.2f\t%i\n" % (record_id,
                                                                                                           record_len,
                                                                                                           mean_coverage,
                                                                                                           median_coverage,
                                                                                                           min_coverage,
                                                                                                           max_coverage,
                                                                                                           coveraged_position_number,
                                                                                                           zero_coverage_position_number,
                                                                                                           zero_coveraged_region_number,
                                                                                                           leading_zero_coverage_len,
                                                                                                           trailing_zero_coverage_len,
                                                                                                           ",".join(zero_coverage_coordinates_list),
                                                                                                           mean_coverage_without_zero_coverage_ends,
                                                                                                           median_coverage_without_zero_coverage_ends,
                                                                                                           length_without_zerocoveraged_ends
                                                                                                           ))

    def collapse_bedgraph(self, input_file, output_file):
        with self.metaopen(output_file, "w") as out_fd:
            line_list_generator = self.file_line_as_list_generator(input_file, separator="\t")

            prev_scaffold, prev_start, prev_end = line_list_generator.readline()[:3]
            prev_start = int(prev_start)
            prev_end = int(prev_end)

            for line in line_list_generator:

                scaffold = line[0]
                start = int(line[1])
                end = int(line[2])

                if (scaffold != prev_scaffold) or (start != prev_end):
                    out_fd.write("%s\t%i\t%i\n" % (prev_scaffold, prev_start, prev_end))
                    prev_scaffold = scaffold
                    prev_start = start
                    prev_end = end
                else:
                    prev_end = end

            out_fd.write("%s\t%i\t%i\n" % (prev_scaffold, prev_start, prev_end))







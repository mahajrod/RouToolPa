#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import numpy as np
from RouToolPa.Routines import MathRoutines
from RouToolPa.Tools.Abstract import Tool


class PurgeDups(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "augustus", path=path, max_threads=max_threads)

    def convert_coverage_file_to_bed(self, input_file, output_prefix):
        length_dict = {}
        coverage_dict = {}
        mean_coverage_dict = {}
        median_coverage_dict = {}

        with self.metaopen(input_file, "r", buffering=100000000) as in_fd, \
                self.metaopen(output_prefix + ".bed", "w", buffering=100000000) as out_fd:
            scaffold, length = in_fd.readline()[1:].split()
            length_dict[scaffold] = int(length)
            coverage_dict[scaffold] = {}

            for line in in_fd:
                if line[0] == ">":
                    scaffold, length = line[1:].split()
                    length_dict[scaffold] = int(length)
                    coverage_dict[scaffold] = {}

                    continue

                #print(line)
                value_list = list(map(int, line.strip().split()))
                value_list[0] -= 1  # convert to zero-based and  half open coordinates
                out_fd.write("{0}\t{1}\n".format(scaffold, "\t".join(map(str, value_list))))
                #print(value_list)
                if value_list[-1] not in coverage_dict[scaffold]:
                    coverage_dict[scaffold][value_list[-1]] = value_list[1] - value_list[0]
                else:
                    coverage_dict[scaffold][value_list[-1]] += value_list[1] - value_list[0]

        for scaffold in coverage_dict:
            median_coverage_dict[scaffold] = MathRoutines.median_from_dict(coverage_dict[scaffold])
            mean_coverage_dict[scaffold] = MathRoutines.mean_from_dict(coverage_dict[scaffold])
        stat_df = pd.DataFrame.from_dict(length_dict, columns=["length", ], orient='index').sort_values(by=["length"], ascending=False)
        stat_df.index.name = "scaffold"
        stat_df["mean_cov"] = pd.Series(mean_coverage_dict)
        stat_df["median_cov"] = pd.Series(median_coverage_dict)
        stat_df.to_csv(output_prefix + ".stat", sep="\t", header=False, index=True)
        stat_df[["length"]].to_csv(output_prefix + ".len", sep="\t", header=False, index=True)

        return stat_df

    def add_lengths_to_dups_bed(self, input_file, length_file, output_file):
        if isinstance(length_file, (str, Path)):
            length_df = pd.read_csv(length_file, sep="\t", header=None, index_col=0, names=["scaffold", "length"])
        else:
            length_df = length_file
        dups_bed_df = pd.read_csv(input_file, sep="\t", header=None, index_col=0, names=["scaffold", "start", "end", "type", "overlapping_scaffold"])

        dups_bed_df["overlap_len"] = dups_bed_df["end"] - dups_bed_df["start"]

        dups_bed_df["scaffold_len"] = length_df["length"]
        dups_bed_df["overlapping_scaffold_len"] = dups_bed_df["overlapping_scaffold"].apply(lambda s: length_df.loc[s, "length"])

        dups_bed_df["overlap_faction,%"] = dups_bed_df["overlap_len"] / dups_bed_df["scaffold_len"]
        dups_bed_df["overlap_faction_overlapping_scaffold,%"] = dups_bed_df["overlap_len"] / dups_bed_df["overlapping_scaffold_len"]
        with open(output_file, "w") as out_fd:
            out_fd.write("#{0}\n".format("\t".join(["scaffold", "start", "end", "type", "overlapping_scaffold",
                                                    "overlap_len", "scaffold_len", "overlapping_scaffold_len",
                                                    "overlap_faction,%", "overlap_faction_overlapping_scaffold,%"])))
            dups_bed_df.to_csv(out_fd, sep="\t", header=False, index=True, na_rep=".")

        def count_fraction(df):
            scaffold_len = df["scaffold_len"].iloc[0]

            sorted_df = df[["start", "end"]].sort_values(by=["start", "end"])

            fraction_df = list(sorted_df.iloc[0])
            print(sorted_df)
            for row in sorted_df.itertuples(index=False):
                if row[0] <= fraction_df[-1][1]:
                    if row[1] > fraction_df[-1][1]:
                        fraction_df[-1][1] = row[1]
                else:
                    fraction_df.append(list(row))

            fraction_df = pd.DataFrame(fraction_df, columns=["start", "end"])
            fraction_df["fraction"] = (fraction_df["end"] - fraction_df["start"]) / scaffold_len

            return sum(fraction_df["fraction"])

        haplo_fraction_df = dups_bed_df[["start", "end", "scaffold_len"]].groupby(by='scaffold').apply(count_fraction)
        print(haplo_fraction_df)

        return dups_bed_df
    """
    def count_contig_fraction_in_haplotype(self, input_file_with_len, output_file):
        if isinstance(input_file_with_len, (str, Path)):
            haplo_df = pd.read_csv(input_file_with_len, sep="\t", header=None, index_col=0,
                                   names=["scaffold", "start", "end", "type", "overlapping_scaffold",
                                          "overlap_len", "scaffold_len", "overlapping_scaffold_len",
                                          "overlap_faction,%", "overlap_faction_overlapping_scaffold,%"])
        else:
            haplo_df = input_file_with_len

        print(haplo_df)
    """


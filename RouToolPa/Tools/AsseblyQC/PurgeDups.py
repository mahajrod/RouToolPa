#!/usr/bin/env python

from pathlib import Path
import pandas as pd
from RouToolPa.Tools.Abstract import Tool


class PurgeDups(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "augustus", path=path, max_threads=max_threads)

    def convert_coverage_file_to_bed(self, input_file, output_prefix):
        length_dict = {}
        with self.metaopen(input_file, "r") as in_fd, self.metaopen(output_prefix + ".bed", "w") as out_fd:
            scaffold, length = in_fd.readline()[1:].split()
            length_dict[scaffold] = length
            for line in in_fd:
                if line[0] == ">":
                    scaffold = in_fd.readline()[1:].split()[0]
                    continue
                out_fd.write(scaffold + "\t" + line)
        length_df = pd.DataFrame(length_dict, columns=["scaffold", "length"])
        #length_df.index.name = "scaffold"
        length_df.to_csv(output_prefix + ".len", sep="\t", header=False, index=True)
        return length_df

    def add_lengths_to_dups_bed(self, input_file, length_file, output_file):
        if isinstance(length_file, [str, Path]):
            length_df = pd.read_csv(length_file, sep="\t", header=None, index_col=0, names=["scaffold", "length"])
        else:
            length_df = length_file
        dups_bed_df = pd.read_csv(input_file, sep="\t", header=None, index_col=0, names=["scaffold", "start", "end", "type", "overlapping_scaffold"])

        dups_bed_df["overlap_len"] = dups_bed_df["end"] - dups_bed_df["start"]

        dups_bed_df["scaffold_len"] = length_df["length"]
        dups_bed_df["overlapping_scaffold_len"] = length_df["length"]

        dups_bed_df.to_csv(output_file, sep="\t", header=False, index=True, na_rep=".")




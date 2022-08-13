__author__ = 'mahajrod'
from pathlib import Path
from functools import partial
from collections import OrderedDict

import pandas as pd

from RouToolPa.Routines.MultipleAlignment import MultipleAlignmentRoutines

class MGIRoutines:
    def __init__(self):

        self.sample_file_headers = OrderedDict({"sample_id": "#SampleID",
                                                "i5": "i5",
                                                "i7": "i7",
                                               })

        pass
    @staticmethod
    def parse_samples(sample_file):
        return pd.read_csv(sample_file, sep="\t", header=0, index_col=0)

    @staticmethod
    def extend_indexes(index_row, max_i5, max_i7):
        return index_row[0] + ("N" * (max_i5 - len(index_row[0]))), \
               index_row[1] + ("N" * (max_i7 - len(index_row[1])))

    @staticmethod
    def get_hamming_distances_for_indexes(index_array, index_label):
        index_number = len(index_array)
        hamming_df = pd.DataFrame(pd.NA, index=index_array, columns=index_array, dtype="Int8")
        hamming_df.index.name = index_label
        for index_1 in range(0, index_number):
            for index_2 in range(index_1 + 1, index_number):
                hamming_df.loc[index_array[index_1], index_array[index_2]] = MultipleAlignmentRoutines.hamming(index_array[index_1],
                                                                                                                       index_array[index_2],
                                                                                                                        any_symbol="N")
                hamming_df.loc[index_array[index_2], index_array[index_1]] = hamming_df.loc[index_array[index_1],
                                                                                                index_array[index_2]]

        return hamming_df

    def analyze_sample_df(self, sample_df, output_prefix, max_index_len=10, max_index_errors=2):
        sample_index_len_df = sample_df.applymap(len)
        sample_index_len_stats_df = sample_index_len_df.agg(['min','max'], axis=0) # calculate min and max values for i5 and i7 indexes

        sample_index_len_df.to_csv(output_prefix + ".sample_index_len.tsv", sep="\t", header=True, index=True)
        sample_index_len_stats_df.to_csv(output_prefix + ".sample_index_len_stats.tsv", sep="\t", header=True, index=True)

        if (sample_index_len_stats_df > 10).values.sum() > 0: # check if indexes are too long
            raise ValueError ("Error!!! One or more indexes are longer than allowed (%i)!" % max_index_len)

        max_i5_len = max(sample_index_len_stats_df[self.sample_file_headers["i5"]])
        max_i7_len = max(sample_index_len_stats_df[self.sample_file_headers["i7"]])
        sample_extended_df = sample_df.apply(partial(self.extend_indexes, max_i5=max_i5_len, max_i7=max_i7_len),
                                             axis=1,
                                             result_type="broadcast")
        sample_extended_df.to_csv(output_prefix + ".extended_index.tsv", sep="\t", header=True,index=True)

        # create df with samples indexes used as indexes for dataframe
        index_df = sample_extended_df.reset_index().set_index([self.sample_file_headers["i5"], self.sample_file_headers["i7"]])
        index_df.to_csv(output_prefix + ".index.tsv", sep="\t", header=True,index=True)

        # calculate hamming distances for
        #print(sample_extended_df[self.sample_file_headers["i5"]].unique())
        unique_i5 = sorted(sample_extended_df[self.sample_file_headers["i5"]].unique())
        unique_i7 = sorted(sample_extended_df[self.sample_file_headers["i7"]].unique())

        hamming_i5_df = self.get_hamming_distances_for_indexes(unique_i5, index_label="i5")
        hamming_i7_df = self.get_hamming_distances_for_indexes(unique_i7, index_label="i7")
        hamming_i5_df.to_csv(output_prefix + ".i5.unique.hamming.tsv", sep="\t", header=True, index=True, na_rep=".")
        hamming_i7_df.to_csv(output_prefix + ".i7.unique.hamming.tsv", sep="\t", header=True, index=True, na_rep=".")

        stats_hamming_i5_df = hamming_i5_df.agg(['min'], axis=1)
        stats_hamming_i5_df.index.name = "i5"
        stats_hamming_i5_df.columns = ["min_pw_dist"]
        stats_hamming_i5_df["max_errors_allowed"] = ((stats_hamming_i5_df["min_pw_dist"] - 1) // 2).apply(lambda x: x if x <= max_index_errors else max_index_errors)
        stats_hamming_i7_df = hamming_i7_df.agg(["min"], axis=1)
        stats_hamming_i7_df.index.name = "i7"
        stats_hamming_i7_df.columns = ["min_pw_dist"]
        stats_hamming_i7_df["max_errors_allowed"] = ((stats_hamming_i7_df["min_pw_dist"] - 1) // 2).apply(lambda x: x if x <= max_index_errors else max_index_errors)
        stats_hamming_i5_df.to_csv(output_prefix + ".i5.unique.stats_hamming.tsv", sep="\t", header=True, index=True)
        stats_hamming_i7_df.to_csv(output_prefix + ".i7.unique.stats_hamming.tsv", sep="\t", header=True, index=True)

        return index_df, max_i5_len, max_i7_len, unique_i5, unique_i7, stats_hamming_i5_df, stats_hamming_i7_df

    def demultiplex_fastq_with_barcodes(self, forward_reads, reverse_reads, sample_file, output_prefix, max_index_len=10, max_index_errors=2,
                                        i5_index_start=150, i7_index_start=160):

        sample_df = self.parse_samples(sample_file)
        index_df, max_i5_len, max_i7_len, unique_i5, unique_i7, stats_hamming_i5_df, stats_hamming_i7_df = self.analyze_sample_df(sample_df,
                                                                                                                                  output_prefix,
                                                                                                                                  max_index_len=max_index_len,
                                                                                                                                  max_index_errors=max_index_errors)
        output_prefix_path = Path(output_prefix)
        output_dir = output_prefix_path.parent
        sample_fd_dict = OrderedDict()
        for sample_id in list(sample_df.index) + ["Undetermined"]:
            sample_fd_dict[sample_id]["forward"] = open(output_dir / (sample_id + "_1.fastq"), "w")
            sample_fd_dict[sample_id]["reverse"] = open(output_dir / (sample_id + "_2.fastq"), "w")
            sample_fd_dict[sample_id]["index"] = open(output_dir / (sample_id + "_I.fastq"), "w")
        with open(forward_reads, "r") as in_for_fd, open(reverse_reads, "r") as in_rev_fd:
            for forward_line in in_for_fd:
                reverse_line = in_rev_fd.readline()
                reverse_seq = in_rev_fd.readline()
                i5_index = reverse_seq[i5_index_start:i5_index_start + max_i5_len]
                i7_index = reverse_seq[i7_index_start:i7_index_start + max_i7_len]

            pass

        for sample_id in sample_fd_dict:
            sample_fd_dict[sample_id]["forward"].close()
            sample_fd_dict[sample_id]["reverse"].close()
            sample_fd_dict[sample_id]["index"].close()

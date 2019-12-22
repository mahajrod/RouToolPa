#!/usr/bin/env python

import numpy as np
import pandas as pd
from collections import OrderedDict

from RouToolPa.Parsers.Abstract import Parser


class CollectionCDHIT(Parser):

    def __init__(self, input_file=None):
        self.records = None
        if input_file:
            self.read(input_file)

    def read(self, input_file):
        # column names
        column_names = ("cluster_id", "sequence_id", "order", "length", "identity")
        line_list = []
        with open(input_file, "r") as in_fd:
            for line in in_fd:

                if line[0] == ">":
                    current_cluster_id = line.strip().split()[-1]
                else:
                    tmp_list = line.strip().split()
                    line_list.append([current_cluster_id,
                                      tmp_list[2][1:-3],
                                      int(tmp_list[0]),
                                      int(tmp_list[1][:-3]),
                                      float(tmp_list[-1][:-1]) if tmp_list[-1] != "*" else np.nan])

        self.records = pd.DataFrame.from_records(line_list,
                                                 index="cluster_id",
                                                 columns=("cluster_id", "sequence_id", "order", "length", "identity"))

    def write(self, output, format="tab"):
        if format == "tab":
            self.records.to_csv(output, sep="\t", na_rep="*")
        elif format == "fam":
            self.write_fam(output)
        else:
            raise ValueError("ERROR!!! Unrecognized output format!")

    def write_fam(self, output):
        with open(output, "w") as out_fd:
            for cluster_id in self.records.index.get_level_values(level=0).unique():
                seq_ids = self.records["sequence_id"].loc[cluster_id]
                if isinstance(seq_ids, str):
                    seq_ids = [seq_ids]
                out_fd.write("%s\t%s\n" % (cluster_id, ",".join(seq_ids)))

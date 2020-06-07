#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
from _collections import OrderedDict

import pandas as pd

from RouToolPa.Tools.Abstract import Tool


class Salmon(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "salmon", path=path, max_threads=max_threads)

    @staticmethod
    def combine_per_sample_counts(quant_filelist, sample_list, output_prefix, tpm_threshold_list=(1.0, 0.5, 0.1)):

        data_dict = OrderedDict()
        columns_dict = OrderedDict()
        tpm_columns_list = []

        for filename, sample_name in zip(quant_filelist, sample_list):
            data_dict[sample_name] = pd.read_csv(filename, sep="\t", index_col="Name")
            columns_dict[sample_name] = ["%s.EffectiveLength" % sample_name,
                                         "%s.TPM" % sample_name,
                                         "%s.NumReads" % sample_name]
            tpm_columns_list.append("%s.TPM" % sample_name)
            data_dict[sample_name].columns = ["Length",] + columns_dict[sample_name]

        combined_df = data_dict[sample_list[0]][data_dict[sample_list[0]].columns]
        if len(sample_list) > 1:
            for sample_name in sample_list[1:]:
                combined_df[columns_dict[sample_name]] = data_dict[sample_name][columns_dict[sample_name]]

        combined_df.to_csv("%s.combined.quant.sf" % output_prefix, sep="\t")

        for threshold in tpm_threshold_list:
            output_file = "%s.quant.TPM%fplus.ids" % threshold
            combined_df[(combined_df[tpm_columns_list] >= threshold).sum(axis=1) >= 1].to_csv(output_file,
                                                                                              sep="\t",
                                                                                              columns=[],
                                                                                              header=False)

        return combined_df



#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
from RouToolPa.Tools.Abstract import Tool
t


class BCL2FASTQ(Tool):
    """
    http://hmmer.janelia.org/
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bcl2fastq", path=path, max_threads=max_threads)

    def hmmbuild(self, input, output, sample_name, summary_output="hmmbuild_summary.t",
                 data_type="protein", max_insert_len=None, window_length=None,
                 alt_model_constact_strategy=None, alt_rel_seq_weigthing_strategy=None,
                 alt_eff_weighting_strategy=None, alt_prior_strategy=None, ):

        options = " --cpu %i" % self.threads
        options += " --dna" if data_type == "dna" else " --rna" if data_type == "rna" else " --amino"
        options += " -n %s" % sample_name
        options += " --maxinsertlen %i" % max_insert_len if max_insert_len else ""
        options += " --w_length %i" % window_length if window_length else ""
        options += " --%s" % alt_model_constact_strategy if alt_model_constact_strategy else ""
        options += " --%s" % alt_rel_seq_weigthing_strategy if alt_rel_seq_weigthing_strategy else ""
        options += " --%s" % alt_eff_weighting_strategy if alt_eff_weighting_strategy else ""
        options += " --%s" % alt_prior_strategy if alt_prior_strategy else ""
        options += " -o %s" % summary_output
        options += " %s" % output
        options += " %s" % input

        self.execute(options, cmd="hmmbuild")



if __name__ == "__main__":
    pass

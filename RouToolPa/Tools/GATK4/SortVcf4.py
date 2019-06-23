#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from RouToolPa.Tools.Abstract import Tool


class SortVcf4(Tool):

    def __init__(self, max_threads=4,  max_memory="4g"):
        Tool.__init__(self, "gatk SortVcf", max_threads=max_threads,
                      max_memory=max_memory)

    def sort_vcf(self, input_file, output_file, seq_dict=None,
                 handling_mode="local",
                 max_memory_per_node=None,
                 job_name=None,
                 log_prefix=None,
                 error_log_prefix=None,
                 modules_list=None,
                 environment_variables_dict=None,
                 max_running_time=None):

        input_file_list = self.make_list_of_path_to_files_from_string(input_file) if isinstance(input_file, str) else self.make_list_of_path_to_files(input_file)

        options = ""

        for filename in input_file_list:
            options += " -I %s" % filename

        options += " -O %s" % output_file
        options += " -SD %s" % seq_dict if seq_dict else ""

        if handling_mode == 'local':
            self.execute(options=options,
                         cmd=("gatk --java-options -Xmx%s SortVcf" % self.max_memory) if self.max_memory else None)
        elif handling_mode == 'slurm':
            self.timelog = None
            slurm_cmd = self.execute(options=options, generate_cmd_string_only=True,
                                     cmd=("gatk --java-options -Xmx%s SortVcf" % self.max_memory) if self.max_memory else None)

            job_id = self.slurm_run_job(job_name, log_prefix, slurm_cmd, error_log_prefix,
                                        "".join(self.split_filename(output_file)[:2]) + ".slurm",
                                        modules_list=modules_list,
                                        environment_variables_dict=environment_variables_dict,
                                        max_running_time=max_running_time,
                                        max_memory_per_node=max_memory_per_node)

            return job_id

if __name__ == "__main__":
    pass

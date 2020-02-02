#!/usr/bin/env python
import os
from RouToolPa.Tools.Abstract import Tool


class FastqDump(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "fastq-dump", path=path, max_threads=max_threads)

    def parse_options(self, split_pe=True, retain_original_ids=True):

        options = ""
        options += " --split-3" if split_pe else ""
        options += " --origfmt" if retain_original_ids else ""

        return options

    def parallel_download(self, sra_id_list, output_dir, split_pe=True, retain_original_ids=True):

        common_options = self.parse_options(split_pe=split_pe, retain_original_ids=retain_original_ids)

        sra_id_list = [sra_id_list] if isinstance(sra_id_list, str) else sra_id_list

        option_list = []
        for sra_id in sra_id_list:
            sra_id_dir = "%s/%s/" % (output_dir, sra_id)
            sra_id_prefix = "%s/%s" % (sra_id_dir, sra_id)
            try:
                os.mkdir(sra_id_dir)
            except OSError:
                pass

            option = common_options
            option += " -O %s" % sra_id_dir
            option += " %s" % sra_id
            option += " > %s.stats 2>%s.error" % (sra_id_prefix, sra_id_prefix) 

            option_list.append(option)

        self.parallel_execute(option_list)

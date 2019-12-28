#!/usr/bin/env python
import os
from RouToolPa.Tools.Abstract import Tool


class FasterqDump(Tool):
    def __init__(self, path="", max_threads=4, tmp_dir=None):
        """
        Spliting sra into forward and reverse fastqs is default now.

        Here are some important differences to fastq-dump:

            The -Z|--stdout option does not work for split-3 and split-files. The tool will fall back to producing files in these cases.

            There is no --gzip|--bizp2 option, you have to compress your files explicitly after they have been written.

            There is no -A option for the accession; just specify the accession or the absolute path directly.

            fasterq-dump does not take multiple accessions, just one.

            There is no -N|--minSpotId and no -X|--maxSpotId option. fasterq-dump version 2.9.1 processes always the whole accession, although it may support partial access in future versions.

            fasterq-dump.2.9.1 was released only on Linux platforms, and is intended primarily for use under larger installations. Our next release should be available for Windows and MacOS.


        :param path:
        :param max_threads:
        """
        Tool.__init__(self, "fasterq-dump", path=path, max_threads=max_threads, tmp_dir=tmp_dir)

    def parse_options(self, tmp_dir=None):

        options = " --threads %i" % self.threads if self.threads else ""
        options += " -t %s" % self.tmp_dir if tmp_dir else ""

        return options

    def parallel_download(self, sra_id_list, output_dir, tmp_dir=None):

        common_options = self.parse_options(tmp_dir=tmp_dir)

        sra_id_list = [sra_id_list] if isinstance(sra_id_list, str) else sra_id_list

        option_list = []
        for sra_id in sra_id_list:
            sra_id_dir = "%s/%s/" % (output_dir, sra_id)
            try:
                os.mkdir(sra_id_dir)
            except OSError:
                pass

            option = common_options
            option += " -O %s" % sra_id_dir
            option += " %s" % sra_id

            option_list.append(option)

        self.parallel_execute(option_list)

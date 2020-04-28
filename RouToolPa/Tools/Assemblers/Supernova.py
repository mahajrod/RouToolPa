#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from RouToolPa.Tools.Abstract import Tool


class Supernova(Tool):
    def __init__(self, path="", max_threads=4, max_memory=512):
        Tool.__init__(self, "supernova", path=path, max_threads=max_threads, max_memory=max_memory)

    def generate_fasta(self, assembly_dir, output_prefix, min_length=1000, header_style="full"):

        options = " mkoutput"
        options += " --asmdir %s" % assembly_dir
        options += " --minsize %i" % min_length
        options += " --headers %s" % header_style

        options_list = []

        for out_type in "raw", "megabubbles", "pseudohap", "pseudohap2":
            tmp_options = options

            if out_type == "pseudohap2":
                tmp_options += " --index"

            tmp_options += " --style %s" % out_type
            tmp_options += " --outprefix %s.min_%i.%s" % (output_prefix, min_length, out_type)

            options_list.append(tmp_options)

        self.parallel_execute(options_list=options_list)

    def assembly(self, input_fastq_dir, output_dir, output_prefix, max_reads="all", disable_ui=True,
                 min_seq_length=150):

        options = " --localcores %i" % self.threads
        options += " --localmem %i" % self.max_memory
        options += " --maxreads %s" % str(max_reads)

        options += " --disable-ui" if disable_ui else ""

        options += " --id %s" % output_dir
        options += " --fastqs %s" % input_fastq_dir

        self.execute(options=options, cmd="supernova run")

        fasta_dir = "%s/fasta/" % output_dir
        fasta_prefix = "%s/%s" % (fasta_dir, output_prefix)
        assembly_dir = "%s/outs/assembly/" % output_dir
        self.safe_mkdir(fasta_dir)

        self.generate_fasta(assembly_dir, fasta_prefix, min_length=1000, header_style="full")
        self.generate_fasta(assembly_dir, fasta_prefix, min_length=min_seq_length, header_style="full")

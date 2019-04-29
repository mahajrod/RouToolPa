#!/usr/bin/env python

from RouToolPa.Tools.Abstract import Tool
import RouToolPa.Formats.AlignmentFormats as AlignmentFormats

class LAST(Tool):
    """
    Class for samtools 1.0+
    Several subcommands are not implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "last", path=path, max_threads=max_threads)

    def parse_lastdb_options(self, db_prefix, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                             verbose=True, keep_preliminary_masking=True,
                             mask_simple_repeats=True):
        options = " -P %i" % self.threads
        options += " -c" if softmasking else ""
        options += " -u %s" % seeding_scheme if seeding_scheme else ""
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)

        options += " -v" if verbose else ""

        options += " %s" % db_prefix
        options += " %s" % (input_fasta_list if isinstance(input_fasta_list, str) else " ".join(input_fasta_list))

        return options

    def create_last_db(self, db_prefix, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                       verbose=True, keep_preliminary_masking=True, mask_simple_repeats=True):

        options = self.parse_lastdb_options(db_prefix, input_fasta_list=input_fasta_list,
                                            softmasking=softmasking,
                                            seeding_scheme=seeding_scheme, verbose=verbose,
                                            keep_preliminary_masking=keep_preliminary_masking,
                                            mask_simple_repeats=mask_simple_repeats)

        self.execute(options=options, cmd="lastdb")

    def parse_lastal_options(self, lastdb, query, output, verbose=True,
                             keep_preliminary_masking=True, mask_simple_repeats=True,
                             output_format="MAF", per_thread_memory="4G"):

        options = " -P %i" % self.threads

        options += " -v" if verbose else ""
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)
        options += " -f %s" % output_format if output_format else ""
        options += " -i %s" % per_thread_memory if per_thread_memory else ""

        options += " %s" % lastdb
        options += " %s" % query
        if output_format == "MAF":
            output_filename_list = self.split_filename(output)
            tab_filename = output + ".tab" if output_filename_list[-1] != "maf" else output_filename_list[0] + output_filename_list[1] + ".tab"
            maf_filename = output + ".maf" if output_filename_list[-1] != "maf" else output
            options += " | tee %s | maf-convert tab > %s" % (maf_filename, tab_filename)
        else:
            options += " > %s" % output
        return options

    def lastal(self, lastdb, query, output, verbose=True,
               keep_preliminary_masking=True, mask_simple_repeats=True,
               output_format="MAF", per_thread_memory="4G"):

        options = self.parse_lastal_options(lastdb, query, output, verbose=verbose,
                                            keep_preliminary_masking=keep_preliminary_masking,
                                            mask_simple_repeats=mask_simple_repeats,
                                            output_format=output_format,
                                            per_thread_memory=per_thread_memory)

        self.execute(options=options, cmd="lastal")

    def convert_maf(self, maf_file, output_file, format="TAB"):

        options = " %s" % format
        options += " %s" % maf_file
        options += " > %s" % output_file

        self.execute(options=options, cmd="maf-convert")

    def filter_tab_output_by_length(self, input_tab, min_hit_len, output_tab, mode="target", verbose=False):

        input_generator = self.file_line_as_list_generator(input_tab)
        read_line_counter = 0
        written_line_counter = 0
        with self.metaopen(output_tab, "w") as out_fd:
            for line_list in input_generator:
                read_line_counter += 1
                if int(line_list[AlignmentFormats.LAST_TAB_COLS["%s_hit_len" % mode]]) >= min_hit_len:
                    out_fd.write("\t".join(line_list))
                    out_fd.write("\n")
                    written_line_counter += 1
        if verbose:
            print("Read %i lines" % read_line_counter)
            print("Written %i lines" % written_line_counter)

    def plot(self, alignment, output, first_genome_seq_id_list=None,
             second_genome_seq_id_list=None,
             first_genome_seq_order=None,
             second_genome_seq_order=None,
             xsize=None, ysize=None, ):

        options = ""

        if first_genome_seq_id_list:
            for seq_id in first_genome_seq_id_list:
                options += " -1 %s" % seq_id

        if second_genome_seq_id_list:
            for seq_id in second_genome_seq_id_list:
                options += " -2 %s" % seq_id

        options += " -x %i" % xsize
        options += " -y %i" % ysize
        options += " %s" % alignment
        options += " %s" % output

        sequence_order_dict = {
                               "input":     "0",
                               "name":      "1",
                               "length":    "2",
                               "alignment": "3"
                               }

        if first_genome_seq_order:

            options += " --sort1 %s" % sequence_order_dict[first_genome_seq_order]

        if second_genome_seq_order:
            options += " --sort2 %s" % sequence_order_dict[second_genome_seq_order]

        self.execute(options=options, cmd="last-dotplot")



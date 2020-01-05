#!/usr/bin/env python

import datetime
from RouToolPa.Tools.Abstract import Tool
import RouToolPa.Formats.AlignmentFormats as AlignmentFormats


class LAST(Tool):
    """
    Class for samtools 1.0+
    Several subcommands are not implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "last", path=path, max_threads=max_threads)

        self.sequence_order_dict = {
                                    "input":     "0",
                                    "name":      "1",
                                    "length":    "2",
                                    "alignment": "3"
                                    }

    def parse_lastdb_options(self, db_prefix, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                             verbose=True, keep_preliminary_masking=True,
                             mask_simple_repeats=True):
        options = " -P %i" % self.threads
        options += " -c" if softmasking else ""
        options += " -u %s" % seeding_scheme if seeding_scheme else ""
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)

        options += " -v" if verbose else ""

        options += " %s.%s.R%i%i.%s" % (db_prefix, seeding_scheme,
                                        1 if keep_preliminary_masking else 0,
                                        1 if mask_simple_repeats else 0,
                                        "soft" if softmasking else "no_soft")
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

    def parse_lastal_options(self, lastdb, query, output_prefix, verbose=True,
                             keep_preliminary_masking=True, mask_simple_repeats=True,
                             output_format="MAF", per_thread_memory="4G", match_score_matrix=None,
                             eg2_threshold=None, discard_limit=None):

        prefix = "%s.R%i%i.%s%s" % (output_prefix,
                                    1 if keep_preliminary_masking else 0,
                                    1 if mask_simple_repeats else 0,
                                    "EG2_%f." % eg2_threshold if eg2_threshold else "",
                                    "C_%i." % discard_limit if discard_limit else "")

        tab_filename = "%stab" % prefix
        maf_filename = "%smaf" % prefix

        options = " -P %i" % self.threads

        options += " -v" if verbose else ""
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)
        options += " -f %s" % output_format if output_format else ""
        options += " -i %s" % per_thread_memory if per_thread_memory else ""
        options += " -p %s" % match_score_matrix if match_score_matrix else ""

        options += " -E %f" % eg2_threshold if eg2_threshold else ""
        options += " -C %i" % discard_limit if discard_limit else ""

        options += " %s" % lastdb
        options += " %s" % query
        if output_format == "MAF":
            options += " | tee %s | maf-convert tab > %s" % (maf_filename, tab_filename)
        else:
            options += " > %s" % tab_filename
        return options

    def parse_last_train_options(self, lastdb, query, output_prefix, verbose=True,
                                 reverse_complement_symmetry=True, directional_symmetry=True,
                                 insertion_deletion_symmetry=True, eg2_threshold=0.05,
                                 discard_limit=2):

        output = "%s.%s%strain.out" % (output_prefix,
                                       "EG2_%f." % eg2_threshold if eg2_threshold else "",
                                       "C_%i." % discard_limit if discard_limit else "")

        options = " -P %i" % self.threads

        options += " -v" if verbose else ""
        options += " --revsym" if reverse_complement_symmetry else ""
        options += " --matsym" if directional_symmetry else ""
        options += " --gapsym" if insertion_deletion_symmetry else ""
        options += " -E %f" % eg2_threshold if eg2_threshold else ""
        options += " -C %i" % discard_limit if discard_limit else ""

        options += " %s" % lastdb
        options += " %s" % query
        options += " > %s" % output

        return options

    def lastal(self, lastdb, query, output_prefix, verbose=True,
               keep_preliminary_masking=True, mask_simple_repeats=True,
               output_format="MAF", per_thread_memory="4G", match_score_matrix=None,
               eg2_threshold=None, discard_limit=None):

        options = self.parse_lastal_options(lastdb, query, output_prefix, verbose=verbose,
                                            keep_preliminary_masking=keep_preliminary_masking,
                                            mask_simple_repeats=mask_simple_repeats,
                                            output_format=output_format,
                                            per_thread_memory=per_thread_memory, match_score_matrix=match_score_matrix,
                                            eg2_threshold=eg2_threshold, discard_limit=discard_limit)

        self.execute(options=options, cmd="lastal")

    def convert_maf(self, maf_file, output_file, format="TAB"):

        options = " %s" % format
        options += " %s" % maf_file
        options += " > %s" % output_file

        self.execute(options=options, cmd="maf-convert")

    def train(self, lastdb, query, output, verbose=True,
              reverse_complement_symmetry=True, directional_symmetry=True,
              insertion_deletion_symmetry=True, eg2_threshold=0.05,
              discard_limit=2):

        options = self.parse_last_train_options(lastdb, query, output, verbose=verbose,
                                                reverse_complement_symmetry=reverse_complement_symmetry,
                                                directional_symmetry=directional_symmetry,
                                                insertion_deletion_symmetry=insertion_deletion_symmetry,
                                                eg2_threshold=eg2_threshold,
                                                discard_limit=discard_limit)

        self.execute(options=options, cmd="last-train")

    def extract_one_to_one_alignments(self, input_maf, output_prefix, max_mismap_probability=1):
        options = ""
        options += " -m %f" % max_mismap_probability if max_mismap_probability else ""
        options += " %s" % input_maf
        options += " | tee %s.first_split.maf" % output_prefix

        options += " | maf-swap | last-split"
        options += " -m %f" % max_mismap_probability if max_mismap_probability else ""
        options += " | maf-swap"
        options += " > %s.one_to_one.maf" % output_prefix

        self.execute(options, cmd="last-split")

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

        if first_genome_seq_order:

            options += " --sort1 %s" % self.sequence_order_dict[first_genome_seq_order]

        if second_genome_seq_order:
            options += " --sort2 %s" % self.sequence_order_dict[second_genome_seq_order]

        options += " %s" % alignment
        options += " %s" % output

        self.execute(options=options, cmd="last-dotplot")

    def filter_lastall_maf(self, input_maf, output_maf, max_eg2=None, max_e=None, min_score=None,
                           buffering=1000000000):

        aln_line_starts = set(["s", "q", "i"])
        kept_count = 0
        removed_count = 0
        alignment_count = 0
        with self.metaopen(input_maf, "r", buffering=buffering) as in_fd, self.metaopen(output_maf, "w", buffering=buffering) as out_fd:
            for line in in_fd:
                if line[0] == "a":
                    alignment_count += 1
                    if alignment_count % 10000:
                        print("%s\t %i alignments handled" % (str(datetime.datetime.now()), alignment_count))

                    #print line
                    for entry in list(map(lambda s: s.split("="), line.strip().split())):
                        #print entry
                        if entry[0] == "EG2":
                            eg2 = float(entry[1])
                        elif entry[0] == "E":
                            e = float(entry[1])
                        elif entry[0] == "score":
                            score = float(entry[1])

                    if max_eg2 and eg2 > max_eg2:
                        store = False
                        removed_count += 1
                        continue
                    if max_e and e > max_e:
                        store = False
                        removed_count += 1
                        continue
                    if min_score and score < min_score:
                        store = False
                        removed_count += 1
                        continue
                    kept_count += 1
                    store = True
                    out_fd.write(line)

                elif line[0] in aln_line_starts:
                    if store:
                        out_fd.write(line)
                else:
                    out_fd.write(line)

        print("\n%i alignments handled\n %i kept\n %i removed" % (alignment_count, kept_count, removed_count))

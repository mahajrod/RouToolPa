__author__ = 'mahajrod'
import os
import sys
from collections import OrderedDict

if sys.version_info[0] == 3:
    izip = zip
else:
    from itertools import izip

import pandas as pd
from Bio import SearchIO, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from RouToolPa.Routines.Sequence import SequenceRoutines
#from Data.Nucleotides import back_degenerate_nucleotides
from RouToolPa.Collections.General import SynDict, IdList

import numpy as np

class AlignmentRoutines(SequenceRoutines):
    def __init__(self):
        self.psl_query_id_column = 9           # zero-based
        self.psl_target_id_column = 13         # zero-based
        pass

    """
    @staticmethod
    def get_db_ids(search_dict):
        id_set = set()
        for query_id in search_dict:
            for hit in search_dict[query_id]:
                id_set.add(hit.id)
        return id_set

    @staticmethod
    def get_codon_alignment(protein_alignment, nucleotide_seq_dict, codon_alignment_file,
                            protein2cds_accordance_dict=None):
        codon_alignment = {}
        if protein2cds_accordance_dict:
            for record in protein_alignment:
                if record.id not in protein2cds_accordance_dict:
                    print("Corresponding CDS was not found for %s protein" % record.id)
                    return -1
        for record in protein_alignment:
            nucleotide_seq = ""
            i = 0
            for aminoacid in record.seq:
                if aminoacid == "-":
                    nucleotide_seq += "---"
                    continue
                else:
                    nucleotide_seq += str(nucleotide_seq_dict[protein2cds_accordance_dict[record.id] if protein2cds_accordance_dict else record.id].seq[3*i:3*(i+1)])
                    i += 1
            codon_alignment[record.id] = SeqRecord(Seq(nucleotide_seq),
                                                   id=record.id,
                                                   description=record.description,
                                                   name=record.name)
            #print(record.id, record.seq)
        SeqIO.write(list(codon_alignment.values()), codon_alignment_file, "fasta")
        return codon_alignment

    def get_codon_alignment_from_files(self, protein_aln_file, nucleotide_seq_file, codon_alignment_file,
                                       cds2protein_accordance_file=None,
                                       alignment_format="fasta", nucleotide_sequence_format="fasta",
                                       cds_index_file=None, retain_cds_index=False):
        protein_aln_dict = AlignIO.read(protein_aln_file, format=alignment_format)
        nucleotide_seq_dict = SeqIO.index_db(cds_index_file if cds_index_file else "nuc_tmp.idx", nucleotide_seq_file,
                                             format=nucleotide_sequence_format)

        protein2cds_accordance_dict = None
        if cds2protein_accordance_file:
            protein2cds_accordance_dict = SynDict()
            protein2cds_accordance_dict.read(cds2protein_accordance_file, key_index=1, value_index=0)

        self.get_codon_alignment(protein_aln_dict, nucleotide_seq_dict, codon_alignment_file,
                                 protein2cds_accordance_dict=protein2cds_accordance_dict)
        if (not cds_index_file) and (not retain_cds_index):
            os.remove("nuc_tmp.idx")

    @staticmethod
    def merge_alignment(alignment_file_list, merged_alignment_file, coordinates_file, format="fasta"):
        #print("Merging alignments...")
        alignment_list = []
        sequence_lengthes = []
        #print(alignment_file_list)

        alignment_file_list_sorted = sorted(alignment_file_list)
        #print(alignment_file_list_sorted)
        for alignment_file in alignment_file_list_sorted:
            #alignment_file.sort()
            parsed = AlignIO.read(alignment_file, format=format)
            parsed.sort()
            alignment_list.append(parsed)
        merged_alignment = None
        for alignment in alignment_list:
            if not merged_alignment:
                sequence_lengthes.append(alignment.get_alignment_length())
                merged_alignment = alignment
                continue
            #print(alignment)
            sequence_lengthes.append(alignment.get_alignment_length())
            merged_alignment += alignment
        SeqIO.write(merged_alignment, merged_alignment_file, "fasta")
        sequence_coordinates = []
        #
        for seq_length in sequence_lengthes:
            if not sequence_coordinates:
                sequence_coordinates.append((1, seq_length))
                continue
            sequence_coordinates.append((sequence_coordinates[-1][1]+1, sequence_coordinates[-1][1]+seq_length))
        #print(sequence_coordinates)
        with open(coordinates_file, "w") as coord_fd:
            coord_fd.write("#length\tstart\tend\n")
            for coord_tuple in sequence_coordinates:
                coord_fd.write("%i\t%i\t%i\n" % (coord_tuple[1] - coord_tuple[0] + 1, coord_tuple[0], coord_tuple[1]))
        return merged_alignment, sequence_lengthes, sequence_coordinates

    @staticmethod
    def extract_degenerate_sites_from_codon_alignment(alignment, genetic_code_table=1):
        degenerate_codon_set = SequenceRoutines.get_degenerate_codon_set(genetic_code_table)
        number_of_alignments = len(alignment)
        alignment_length = len(alignment[0])
        if alignment_length % 3 > 0:
            raise(ValueError, "Length of alignment is not divisible by 3")
        else:
            number_of_codons = int(alignment_length / 3)
        degenerate_columns = []
        for i in range(0, number_of_codons):
            position_strings = []
            for j in range(0, 3):
                position_strings.append(list(set(alignment[:, 3*i + j])))
            if (len(position_strings[0]) > 1) or (len(position_strings[1]) > 1):
                continue
            ambigious_codon = position_strings[0][0] + position_strings[1][0] + "N"

            #if Seq(ambigious_codon).translate(table=genetic_code_table) == "X":
            #    continue
            #else:
            #    degenerate_columns.append(alignment[:, 3*i + 2])
            #    #print(i*3 +3)

            if ambigious_codon in degenerate_codon_set:
                degenerate_columns.append(alignment[:, 3*i + 2])

        number_of_degenerate_columns = len(degenerate_columns)
        record_list = []
        for i in range(0, number_of_alignments):
            string = ""
            for j in range(0, number_of_degenerate_columns):
                string += degenerate_columns[j][i]
            record = SeqRecord(seq=Seq(string), id=alignment[i].id)
            record_list.append(record)

        #print(number_of_alignments, alignment_length)
        degenerate_alignment = MultipleSeqAlignment(record_list)

        return degenerate_alignment

    def extract_degenerate_sites_from_codon_alignment_from_file(self, alignment_file, output_alignment_file,
                                                                genetic_code_table=1, format="fasta"):
        alignment = AlignIO.read(alignment_file, format=format)
        degenerate_alignment = self.extract_degenerate_sites_from_codon_alignment(alignment,
                                                                                  genetic_code_table=genetic_code_table)

        AlignIO.write([degenerate_alignment], output_alignment_file, format=format)

    @staticmethod
    def sequences_from_alignment_generator(alignments, gap_symbol="-"):
        for record in alignments:
            record.seq = Seq(str(record.seq).replace(gap_symbol, ""))
            yield record

    def extract_sequences_from_alignment(self, alignment_file, output_file, alignment_format="fasta",
                                         output_format="fasta", gap_symbol="-"):
        alignments = AlignIO.read(alignment_file, format=alignment_format)
        SeqIO.write(self.sequences_from_alignment_generator(alignments, gap_symbol=gap_symbol),
                    output_file, format=output_format)
    """

    def filter_psl_by_ids(self, psl_file, output_file,
                          white_query_id_list=(),
                          black_query_id_list=(),
                          white_target_id_list=(),
                          black_target_id_list=()):

        with self.metaopen(psl_file, "r") as in_fd, \
             self.metaopen(output_file, "w") as out_fd:

            for i in range(0, 20):
                line = in_fd.readline()
                if line[:3] == "---":
                    break
            else:
                in_fd.seek(0)

            for line in in_fd:
                tmp_line = line.split("\t")
                if self.check_id(tmp_line[self.psl_query_id_column],
                                 white_query_id_list,
                                 black_query_id_list) and \
                   self.check_id(tmp_line[self.psl_target_id_column],
                                 white_target_id_list,
                                 black_target_id_list):

                    out_fd.write(line)

    def filter_psl_by_ids_from_file(self, psl_file, output_file,
                                    white_query_id_file=None,
                                    black_query_id_file=None,
                                    white_target_id_file=None,
                                    black_target_id_file=None):

        self.filter_psl_by_ids(psl_file, output_file,
                               white_query_id_list=IdList(filename=white_query_id_file) if white_query_id_file else (),
                               black_query_id_list=IdList(filename=black_query_id_file) if black_query_id_file else (),
                               white_target_id_list=IdList(filename=white_target_id_file) if white_target_id_file else (),
                               black_target_id_list=IdList(filename=black_target_id_file) if black_target_id_file else ())

    def filter_psl_by_target_len(self, psl_file, output_file, min_target_len=None, max_target_len=None):
        pass

    @staticmethod
    def parse_search_file(input_file, mode, format="hmmer3-text", index_file=None):
        if mode == "index_db" or ((not isinstance(input_file, str)) and (len(input_file) > 1)):
            index = index_file if index_file else "tmp.idx"
            seq_dict = SearchIO.index_db(index, [input_file] if isinstance(input_file, str) else input_file, format=format)
        elif mode == "index":
            seq_dict = SearchIO.index(input_file if isinstance(input_file, str) else input_file[0], format=format)
        elif mode == "parse":
            seq_dict = OrderedDict()
            for record in SearchIO.parse(input_file if isinstance(input_file, str) else input_file[0], format=format):
                seq_dict[record.id] = record
            #seq_dict = SeqIO.to_dict(SeqIO.parse(input_file if isinstance(input_file, str) else input_file[0], format=format))

        return seq_dict

    def calculate_masking_from_coverage_files(self, coverage_file_list, mean_coverage_list, output_file,
                                              sample_labels=None,
                                              max_threshold=2.5, min_threshold=None,min_sample_number=1,
                                              scaffold_column=0, position_column=1, coverage_column=2):
        """

        :param coverage_file_list:
        :param mean_coverage_list:
        :param output_file:
        :param sample_labels:
        :param max_threshold:
        :param min_threshold:
        :param min_sample_number:
        :param scaffold_column: 0-based
        :param position_column: 0-based
        :param coverage_column: 0-based
        :return:
        """
        file_number = len(coverage_file_list)
        index_list = [i for i in range(0, file_number)]
        if min_sample_number > file_number:
            raise ValueError("ERROR!!! Minimun number of files necessary for treating position as masked is greater "
                             "than total number of files")

        if max_threshold and min_threshold:
            max_coverage_list = np.array(mean_coverage_list) * max_threshold
            min_coverage_list = np.array(mean_coverage_list) * min_threshold

            def check_pos(entry):
                if np.sum(np.logical_and(np.less(entry, min_coverage_list),
                                         np.greater(entry, max_coverage_list))) >= min_sample_number:
                    return True
                return False

        elif max_threshold:
            max_coverage_list = np.array(mean_coverage_list) * max_threshold

            def check_pos(entry):
                if np.sum(np.greater(entry, max_coverage_list)) >= min_sample_number:
                    return True
                return False

        elif min_threshold:
            min_coverage_list = np.array(mean_coverage_list) * min_threshold

            def check_pos(entry):
                if np.sum(np.less(entry, min_coverage_list)) >= min_sample_number:
                    return True
                return False
        else:
            raise ValueError("ERROR!!! Neither minimum nor maximum threshold was set!")

        line_index = 1
        with open(output_file, "w") as out_fd:
            out_fd.write("#scaffold\tposition\t%s\n" % (",".join(sample_labels if sample_labels else coverage_file_list)))

            for line_list in izip(*[self.file_line_as_list_generator(filename) for filename in coverage_file_list]):

                coverage_list = [int(line_list[i][coverage_column]) for i in index_list]
                if check_pos(coverage_list):
                    out_fd.write("%s\t%s\t%s\n" % (line_list[0][scaffold_column],
                                                   line_list[0][position_column],
                                                   ",".join(map(str, coverage_list))))

                if line_index % 100000000 == 0:
                    print("Handled %i million positions..." % (line_index/1000000))

                line_index += 1

    def collapse_per_base_coverage_mask(self, mask_file, out_file,
                                        scaffold_column=0,
                                        position_column=1,
                                        comments_prefix="#",
                                        output_format="0-based",
                                        in_memory=True):

        line_list_generator = self.file_line_as_list_generator(mask_file, comments_prefix=comments_prefix)
        with self.metaopen(out_file, "w") as out_fd:
            tmp = next(line_list_generator)

            prev_scaffold = tmp[scaffold_column]
            prev_start = int(tmp[position_column])
            prev_end = int(tmp[position_column])

            coordinates_df = []
            if in_memory:
                for line_list in line_list_generator:
                    pos = int(line_list[position_column])
                    if (prev_scaffold != line_list[scaffold_column]) or (pos != prev_end + 1):
                        coordinates_df.append((prev_scaffold, prev_start, prev_end))

                        prev_scaffold = line_list[scaffold_column]
                        prev_start = pos
                        prev_end = pos
                    else:
                        prev_end += 1

                coordinates_df.append((prev_scaffold, prev_start, prev_end))
                coordinates_df = pd.DataFrame(coordinates_df, columns=("scaffold", "start", "end"), index="scaffold")

                if output_format == "0-based":
                    coordinates_df["start"] -= 1

                coordinates_df.to_csv(out_fd, sep="\t", index=True)
            else:
                for line_list in line_list_generator:
                    pos = int(line_list[position_column])
                    if (prev_scaffold != line_list[scaffold_column]) or (pos != prev_end + 1):
                        out_fd.write("%s\t%i\t%i\n" % (prev_scaffold, (prev_start - 1) if output_format == "0-based" else prev_start, prev_end))

                        prev_scaffold = line_list[scaffold_column]
                        prev_start = pos
                        prev_end = pos
                    else:
                        prev_end += 1

                out_fd.write("%s\t%i\t%i\n" % (prev_scaffold, (prev_start - 1) if output_format == "0-based" else prev_start, prev_end))

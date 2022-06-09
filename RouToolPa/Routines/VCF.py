import os
import sys

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

from copy import deepcopy
from collections import OrderedDict

import pandas as pd

from RouToolPa.Collections.General import IdList, SynDict
from RouToolPa.Routines.Sequence import SequenceRoutines
from RouToolPa.Formats.VariantFormats import VCF_COLS


class VCFRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)

    def combine_same_samples_vcfs(self, output, vcf_list=None,
                                  order_vcf_files=False, sort=False, chunk_folder=None, chunk_prefix=None,
                                  chunk_suffix=None, starting_chunk=None, chunk_number_list=None,
                                  close_fd_after=False, extension_list=[".vcf", ]):

        output_fd = output if isinstance(output, file) else open(output, "w")

        if vcf_list:
            vcf_files = self.make_list_of_path_to_files_by_extension(vcf_list, extension_list=extension_list,
                                                                     recursive=False, return_absolute_paths=True)

            #print(vcf_files)
            #print( type(vcf_files))
        else:
            pass

        if order_vcf_files:

            vcf_files.sort(key=self.natural_keys_int)

        if sort:
            unsorted_file = "%s.unsorted.tmp" % output

            string = "cat %s > %s" % (vcf_files[0], unsorted_file)
            print(string)
            os.system(string)
            #with open(unsorted_file, "w") as out_fd:
            #    pass

            for filename in vcf_files[1:]:
                string = "sed -n '/^[^#]/p' %s >> %s" % (filename, unsorted_file)
                print(string)
                os.system(string)

            sorting_string = "awk -F'\\t' '{if (substr($1,1,1) == \"#\") {print $0} else {print $0 | \"sort -k1,1 -k2,2n\"} }' %s > %s" % (unsorted_file, output)
            """
            sorting_string = "(sed '/^[^#]/Q' %s; sort -k1,1 -k2,2n %s) > %s" % (vcf_files[0],
                                                                                 unsorted_file,
                                                                                 output)
            """
            print(sorting_string)

            os.system(sorting_string)
            os.remove(unsorted_file)

        else:
            print("Handling file %s ..." % vcf_files[0])
            with open(vcf_files[0], "r") as in_fd:
                for line in in_fd:
                    output_fd.write(line)

            if len(vcf_files) > 1:
                for filename in vcf_files[1:]:
                    print("Handling file %s ..." % filename)
                    with open(filename, "r") as in_fd:
                        for line in in_fd:
                            if line[0] == "#":
                                continue
                            else:
                                output_fd.write(line)
                                break
                        for line in in_fd:
                            output_fd.write(line)

        if not isinstance(output, file) or close_fd_after:
            output_fd.close()

    def check_gvcf_integrity(self, gvcf_file, output_prefix, reference=None, length_dict=None, parsing_mode="parse"):
        len_dict = length_dict if length_dict else self.get_lengths(record_dict=self.parse_seq_file(reference,
                                                                                                    mode=parsing_mode),
                                                                    out_file=None,
                                                                    close_after_if_file_object=False)

        scaffold_dict = OrderedDict()

        with self.metaopen(gvcf_file, "r") as gvcf_fd:
            prev_scaffold = ""

            for line in gvcf_fd:
                #print line
                if line[0] == "#":
                    continue

                line_list = line.split("\t")
                scaffold = line_list[0]
                start = int(line_list[1])
                format = line_list[7].split(";")

                if (len(format) == 1) and (format[0][0:3] == "END"):
                    end = int(format[0].split("=")[1])
                else:
                    end = start + len(line_list[3]) - 1
                #print line_list
                #print scaffold, start, end, format

                if scaffold not in scaffold_dict:
                    scaffold_dict[scaffold] = []

                if scaffold != prev_scaffold:
                    scaffold_dict[scaffold].append([deepcopy(start), deepcopy(end)])
                else:
                    #print scaffold_dict[scaffold][-1][1]
                    if scaffold_dict[scaffold][-1][1] + 1 >= start:
                        scaffold_dict[scaffold][-1][1] = deepcopy(max(end, scaffold_dict[scaffold][-1][1]))
                    else:
                        print(scaffold_dict[scaffold])
                        print(line)
                        scaffold_dict[scaffold].append([deepcopy(start), deepcopy(end)])
                prev_scaffold = scaffold

        complete_scaffolds = IdList()
        fragmented_scaffolds = IdList()
        scaffolds_with_absent_fragments = IdList()

        with open("%s.scaffold_regions" % output_prefix, "w") as scaf_reg_fd:

            for scaffold in scaffold_dict:
                if len(scaffold_dict[scaffold]) > 1:
                    fragmented_scaffolds.append(scaffold)

                scaffold_length = sum(map(lambda s: s[1] - s[0] + 1, scaffold_dict[scaffold]))
                if scaffold_length != len_dict[scaffold]:
                    scaffolds_with_absent_fragments.append(scaffold)
                else:
                    complete_scaffolds.append(scaffold)
                scaf_reg_fd.write("%s\t%s\n" % (scaffold, ",".join(map(lambda s: "-".join(map(str,s)), scaffold_dict[scaffold]))))

        complete_scaffolds.write("%s.complete_scaffolds" % output_prefix)
        fragmented_scaffolds.write("%s.fragmented_scaffolds" % output_prefix)
        scaffolds_with_absent_fragments.write("%s.scaffolds_with_absent_fragments" % output_prefix)

    def extract_heterozygous_variants(self, input_vcf, output_prefix, mode="one", verbose=True):
        with self.metaopen(input_vcf, "r") as in_fd, \
             self.metaopen("%s.hetero.vcf" % output_prefix, "w") as het_fd, \
             self.metaopen("%s.homo.vcf" % output_prefix, "w") as homo_fd:

            het_counter = 0
            homo_counter = 0

            def test_heterozygosity(s):
                if len(s) == 1:
                    return False
                if (s[0] != s[1]) and (s[0] != ".") and (s[1] != "."):
                    return True
                else:
                    return False

            for line in in_fd:
                if line[0] == "#":
                    het_fd.write(line)
                    homo_fd.write(line)
                else:
                    line_list = line.strip().split("\t")
                    info_list = line_list[VCF_COLS["FORMAT"]].split(":")
                    samples = map(lambda s: s.split(":"), line_list[9:])
                    for index in range(0, len(info_list)):
                        if info_list[index] == "GT":
                            genotype_index = index
                            break
                    heterozygous_sample_number = sum(map(test_heterozygosity,
                                                     map(lambda s: s[genotype_index].split("/") if "/" in s[genotype_index] else s[genotype_index].split("|") if "|" in s[genotype_index] else s[genotype_index], samples)))

                    if mode == "one":
                        if heterozygous_sample_number > 0:
                            het_fd.write(line)
                            het_counter += 1
                        else:
                            homo_fd.write(line)
                            homo_counter += 1
                    elif mode == "all":
                        if heterozygous_sample_number == len(samples):
                            het_fd.write(line)
                            het_counter += 1
                        else:
                            homo_fd.write(line)
                            homo_counter += 1

            if verbose:
                print("Heterozygous variants: %i" % het_counter)
                print("Homozygous variants: %i" % homo_counter)
                print("Total variants: %i" % (het_counter + homo_counter))

    @staticmethod
    def extract_per_sample_freq_df(collection_vcf):

        freq_df_list = [collection_vcf.records[["POS", "REF", "ALT"]].copy()]
        freq_df_list[0]["POS"] += 1
        freq_df_list.append(collection_vcf.records[[("INFO", "DP", 0)]].copy())
        for sample in sorted(collection_vcf.samples):
            freq_df_list.append(collection_vcf.records.loc[:, (sample, "AD", slice(None))])
            freq_df_list.append(collection_vcf.records.loc[:, (sample, "AD", slice(None))].divide(
                collection_vcf.records[(sample, "DP", 0)], axis='index').rename(columns={"AD": "ADF"}))
        return pd.concat(freq_df_list, axis="columns")

    @staticmethod
    def filter_freq_df(freq_df, max_alt_allel_freq_minimum, min_total_coverage):
        samples = set(freq_df.columns.levels[0]) - set(["POS", "REF", "ALT", "DP"])
        max_alt_allel_freq = freq_df.loc[:,
                                         pd.IndexSlice[samples, "ADF", range(1, len(freq_df[("ALT", "ALT")].columns) + 1)]].max(
            axis=1)

        return freq_df[
            (max_alt_allel_freq >= max_alt_allel_freq_minimum) & (freq_df[("INFO", "DP", 0)] >= min_total_coverage)]

    @staticmethod
    def save_freq_df_to_xlsx(freq_df, sheet_name, output_file):

        sample_number = len(freq_df.columns.levels[0]) - len(["POS", "REF", "ALT", "DP"])
        alt_allel_column_number = len(freq_df["ALT"].columns)
        allel_column_number = alt_allel_column_number + 1

        hor_shift = len(["POS", "REF", "ALT", "DP"]) + allel_column_number
        ver_shift = 4
        data_row_number = len(freq_df)
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')

        freq_df.to_excel(writer, sheet_name=sheet_name, header=True, index=True, freeze_panes=(ver_shift, hor_shift))
        worksheet = writer.sheets[sheet_name]

        for sample_index in range(0, sample_number):
            sample_ref_freq_column = hor_shift + allel_column_number + sample_index * allel_column_number * 2
            worksheet.write(ver_shift - 2, sample_ref_freq_column, "REF")
            worksheet.conditional_format(ver_shift, sample_ref_freq_column,
                                         ver_shift + data_row_number, sample_ref_freq_column,
                                         {'type': '3_color_scale',
                                          'min_value': 0,
                                          'mid_value': 0.50,
                                          'max_value': 1.00,
                                          'min_type': 'value',
                                          'mid_type': 'value',
                                          'max_type': 'value',
                                          'min_color': '#F8696B',  # red
                                          'mid_color': '#FFEB84',  # yellow
                                          'max_color': '#63BE7B'  # green
                                          })
            worksheet.conditional_format(ver_shift,
                                         sample_ref_freq_column + 1,
                                         ver_shift + data_row_number,
                                         sample_ref_freq_column + alt_allel_column_number,
                                         {'type': '3_color_scale',
                                          'min_value': 0,
                                          'mid_value': 0.50,
                                          'max_value': 1.00,
                                          'min_type': 'value',
                                          'mid_type': 'value',
                                          'max_type': 'value',
                                          'min_color': '#63BE7B',  # green
                                          'mid_color': '#FFEB84',  # yellow
                                          'max_color': '#F8696B'  # red
                                          })

        writer.save()
        return writer

    def add_variant_ids_to_vcf(self, input_vcf, output_vcf, id_prefix=None, retain_old_id=True):
        with self.metaopen(input_vcf, "r") as in_fd, self.metaopen(output_vcf, "w") as out_fd:
            for line in in_fd:
                out_fd.write(line)
                if line[:6] == "#CHROM":
                    break
            variant_counter = 1
            for line in in_fd:
                line_list = line.split("\t")

                if retain_old_id and (line_list[2] != "."):
                    line_list[2] = "%s%s,%s" % (id_prefix, variant_counter, line_list[2])
                else:
                    line_list[2] = "%s%s" % (id_prefix, variant_counter)

                out_fd.write("\t".join(line_list))

                variant_counter += 1



    """ 
    @staticmethod
    def convert_gvcf_to_coverage_file(self, gvcf_file, coverage_file):
        with self.metaopen(coverage_file, "w") as out_fd:
            for line_list in self.file_line_as_list_generator(gvcf_file, comments_prefix="#", separator="\t"):
                scaffold = line
    """
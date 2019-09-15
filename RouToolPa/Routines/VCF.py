import os
import sys

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

from copy import deepcopy
from collections import OrderedDict
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
             self.metaopen("%s.hetero.vcf" % output_prefix) as het_fd, \
             self.metaopen("%s.homo.vcf") as homo_fd:

            het_counter = 0
            homo_counter =0

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
                    heterozygous_sample_number = sum(map(lambda s: True if s[0] != s[1] else False,
                                                     map(lambda s: s[genotype_index].split("/"), samples)))

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
    """ 
    @staticmethod
    def convert_gvcf_to_coverage_file(self, gvcf_file, coverage_file):
        with self.metaopen(coverage_file, "w") as out_fd:
            for line_list in self.file_line_as_list_generator(gvcf_file, comments_prefix="#", separator="\t"):
                scaffold = line
    """
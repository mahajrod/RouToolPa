#!/usr/bin/env python
__author__ = 'mahajrod'

import re
import os
import sys
import bz2
import gzip
import shutil
from collections import OrderedDict
from collections.abc import Iterable

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

from RouToolPa.Collections.General import IdSet,  IdList, SynDict


class FileRoutines:
    def __init__(self):
        self.filetypes_dict = {"fasta": [".fa", ".fasta", ".fa", ".pep", ".cds"],
                               "fastq": [".fastq", ".fq"],
                               "genbank": [".gb", ".genbank"],
                               "newick": [".nwk"],
                               "gz": [".gz"],
                               "bzip": [".bz2"],
                               "bam": [".bam"],
                               "sam": [".sam"],
                               "cram": [".cram"]}

    @staticmethod
    def metaopen(filename, flags, buffering=None, compresslevel=5):
        if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
            if isinstance(filename, file):
                return filename
            else:
                raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
        elif filename[-3:] == ".gz":
            return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
        elif filename[-4:] == ".bz2":
            return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
        else:
            if buffering is not None:
                return open(filename, flags, buffering=buffering)
            else:
                return open(filename, flags)

    @staticmethod
    def add_external_extraction_to_filename(filename):
        if filename[-3:] == ".gz":
            return "<(gunzip -c %s)" % filename
        elif filename[-4:] == ".bz2":
            return "<(bunzip2 -c %s)" % filename
        else:
            return filename

    @staticmethod
    def add_external_extraction_to_filelist(filelist):
        new_filelist = []

        for filename in filelist:
            if filename[-3:] == ".gz":
                new_filelist.append("<(gunzip -c %s)" % filename)
            elif filename[-4:] == ".bz2":
                new_filelist.append("<(bunzip2 -c %s)" % filename)
            else:
                new_filelist.append(filename)

        return new_filelist

    @staticmethod
    def safe_mkdir(dirname, 
                   description_filename=None, description_text=None,
                   readme_filename=None, readme_text=None):
        try:
            os.mkdir(dirname)
        except OSError:
            pass
        
        if not(description_filename is None):
            description_filename = "%s/%s" % (dirname, description_filename)
    
            if not os.path.isfile(description_filename):
                with open(description_filename, "w") as descr_fd:
                    if not (description_text is None):
                        descr_fd.write(description_text)
                        
        if not(readme_filename is None):
            readme_filename = "%s/%s" % (dirname, readme_filename)
    
            if not os.path.isfile(readme_filename):
                with open(readme_filename, "w") as descr_fd:
                    if not (readme_text is None):
                        descr_fd.write(readme_text)

    def recursive_mkdir(self, dir_dict=None, out_dir=None,
                        description_filename=None, description_text=None,
                        readme_filename=None, readme_text=None):
        if not(out_dir is None):
            self.safe_mkdir(out_dir)

        if dir_dict:
            for directory in dir_dict:
                dirname = directory if out_dir is None else "%s/%s" % (out_dir, directory)
                self.safe_mkdir(dirname,
                                description_filename=description_filename,
                                description_text=description_text,
                                readme_filename=readme_filename,
                                readme_text=readme_text)

                if isinstance(dir_dict[directory], dict):
                    self.recursive_mkdir(dir_dict[directory],
                                         out_dir=dirname,
                                         description_filename=description_filename,
                                         description_text=description_text,
                                         readme_filename=readme_filename,
                                         readme_text=readme_text)

    def detect_filetype_by_extension(self, filename, filetypes_dict=None):
        filetypes = filetypes_dict if filetypes_dict else self.filetypes_dict
        directory, prefix, extension = self.split_filename(filename)
        for filetype in filetypes:
            if extension in filetypes[filetype]:
                return filetype
        return None

    @staticmethod
    def check_path(path_to_check):
        #print (path_to_check)
        #returns path with / at end or blank path
        if path_to_check != "":
            if path_to_check[-1] != "/":
                return path_to_check + "/"
        return path_to_check

    @staticmethod
    def check_id(id, white_list=(), black_list=()):
        if white_list and black_list:
            if (id in white_list) and (id not in black_list):
                return True
            return False
        elif white_list:
            return True if (id in white_list) else False
        elif black_list:
            return False if (id in black_list) else True
        else:
            return True

    @staticmethod
    def check_dir_path(path_to_check):
        #print (path_to_check)
        #returns path with / at end or blank path
        if path_to_check != "":
            if path_to_check[-1] != "/":
                return path_to_check + "/"
        return path_to_check

    @staticmethod
    def split_filename(filepath):
        directory, basename = os.path.split(filepath)
        if directory:
            if directory[-1] != "/":
                directory += "/"
        prefix, extension = os.path.splitext(basename)
        return directory, prefix, extension if filepath else None

    @staticmethod
    def get_basename(filepath):
        return os.path.split(filepath)[-1] if filepath else None

    def make_list_of_path_to_files(self, list_of_dirs_and_files, expression=None, recursive=False,
                                   return_absolute_paths=True):
        file_list = []
        for entry in [list_of_dirs_and_files] if isinstance(list_of_dirs_and_files, str) else list_of_dirs_and_files:
            if os.path.isdir(entry):
                files_in_dir = ["%s%s" % (self.check_path(entry), filename)
                                for filename in sorted(filter(expression, os.listdir(entry))
                                                       if expression else os.listdir(entry))]
                if recursive:
                    for filename in files_in_dir:
                        if os.path.isdir(filename):
                            file_list += self.make_list_of_path_to_files([filename],
                                                                         expression=expression,
                                                                         recursive=recursive)
                        else:
                            file_list.append(filename)
                else:
                    file_list += files_in_dir
            elif os.path.exists(entry):
                if expression:
                    if expression(os.path.abspath(entry)):
                        file_list.append(os.path.abspath(entry))

                else:
                    file_list.append(os.path.abspath(entry))
            else:
                print("%s does not exist" % entry)
        # direct conversion to list was added for compatibility with python3
        # in which map function returns map object instead of list
        return list(map(os.path.abspath, file_list)) if return_absolute_paths else file_list

    def make_list_of_path_to_files_from_string(self, input_string, file_separator=",",
                                               expression=None, recursive=False):
        return self.make_list_of_path_to_files(input_string.split(file_separator),
                                               expression=expression,
                                               recursive=recursive)

    @staticmethod
    def split_string_by_separator(string, separator):
        return string.split(separator)

    def split_string_by_comma(self, string):
        return self.split_string_by_separator(string, separator=",")

    def split_string_by_colon(self, string):
        return self.split_string_by_separator(string, separator=":")

    @staticmethod
    def check_extension(filename, extension_list=[]):
        if extension_list:
            for extension in extension_list:
                #print extension, filename[-len(extension):]
                if extension == filename[-len(extension):]:
                    return True
        else:
            #print "CCCCCCCCCCCCCCCCC"
            #print filename
            return True

        return False

    def make_list_of_path_to_files_by_extension(self, list_of_dirs_and_files, extension_list=[], recursive=False,
                                                return_absolute_paths=True):

        def check_extension(filename):
            return self.check_extension(filename, extension_list)

        return self.make_list_of_path_to_files(list_of_dirs_and_files, expression=check_extension,
                                               recursive=recursive, return_absolute_paths=return_absolute_paths)

    @staticmethod
    def read_synonyms_dict(filename, header=False, separator="\t",
                           split_values=False, values_separator=",", key_index=0, value_index=1):
        # reads synonyms from file
        synonyms_dict = OrderedDict()
        with open(filename, "r") as in_fd:
            if header:
                header_str = in_fd.readline().strip()
            for line in in_fd:
                tmp = line.strip().split(separator)
                #print line
                key, value = tmp[key_index], tmp[value_index]
                if split_values:
                    value = value.split(values_separator)
                synonyms_dict[key] = value
        return synonyms_dict

    @staticmethod
    def read_ids(filename, header=False, close_after_if_file_object=False):
        #reads ids from file or file object with one id per line
        id_list = []

        in_fd = filename if isinstance(filename, file) else open(filename, "r")

        if header:
            header_str = in_fd.readline().strip()
        for line in in_fd:
            id_list.append(line.strip())
        if (not isinstance(filename, file)) or close_after_if_file_object:
            in_fd.close()

        return id_list

    @staticmethod
    def read_tsv_as_rows_list(filename, header=False, header_prefix="#", separator="\t"):
        tsv_list = []
        with open(filename, "r") as tsv_fd:
            if header:
                header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
            for line in tsv_fd:
                tsv_list.append(line.strip().split(separator))
        return header_list, tsv_list if header else tsv_list

    @staticmethod
    def read_tsv_as_columns_dict(filename, header_prefix="#", separator="\t"):
        tsv_dict = OrderedDict()
        with open(filename, "r") as tsv_fd:
            header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
            number_of_columns = len(header_list)
            for column_name in header_list:
                tsv_dict[column_name] = []
            for line in tsv_fd:
                tmp_line = line.strip().split(separator)
                for i in range(0, number_of_columns):
                    tsv_dict[header_list[i]].append(tmp_line[i])
        return tsv_dict

    @staticmethod
    def tsv_split_by_column(tsv_file, column_number, separator="\t", header=False, outfile_prefix=None):
        # column number should start from 0
        header_string = None
        splited_name = tsv_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd_dict = {}

        with open(tsv_file, "r") as in_fd:
            if header:
                header_string = in_fd.readline()
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] not in out_fd_dict:
                    print (line_str[column_number])
                    suffix = line_str[column_number].replace(" ", "_")
                    out_name = "%s_%s.%s" % (out_prefix, suffix, extension)
                    out_fd_dict[line_str[column_number]] = open(out_name, "w")
                    if header:
                        out_fd_dict[line_str[column_number]].write(header_string)
                out_fd_dict[line_str[column_number]].write(line)

        for entry in out_fd_dict:
            out_fd_dict[entry].close()

    def tsv_extract_by_column_value(self, tsv_file, column_number, column_value,
                                    output_file, comments_prefix="#", header=False, separator="\t"):
        # column number should start from 0
        # column_value should be string or list of strings

        comments_prefix_len = len(comments_prefix)

        column_value_list = column_value if isinstance(column_value, Iterable) else []
        with self.metaopen(tsv_file, "r") as in_fd, self.metaopen(output_file, "w") as out_fd:
            if header:
                out_fd.write(in_fd.readline())
            for line in in_fd:
                if line[:comments_prefix_len] == comments_prefix:
                    out_fd.write(line)
                line_str = line.strip().split(separator)
                if line_str[column_number] in column_value_list:
                    out_fd.write(line)
        out_fd.close()

    @staticmethod
    def tsv_remove_by_column_value(tsv_file, column_number, column_value, separator="\t",
                                   header=False, outfile_prefix=None):
        # column number should start from 0
        # column_value should be string or list of strings

        splited_name = tsv_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd = open("%s.%s" % (out_prefix, extension), "w")

        column_value_list = column_value if isinstance(column_value, Iterable) else []
        with open(tsv_file, "r") as in_fd:
            if header:
                out_fd.write(in_fd.readline())
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] not in column_value_list:
                    out_fd.write(line)
        out_fd.close()

    @staticmethod
    def p_distance(seq_a, seq_b, seq_len):
        dist = 0
        for i in range(0, seq_len):
            if seq_a[i] != seq_b[i]:
                dist += 1
        return dist

    @staticmethod
    def make_lists_forward_and_reverse_files(sample_dir, filename_fragment_to_mark_se_reads=".se.", input_is_se=False):
        """Legacy function. Use FastQRoutines.make_fastq_lists() instead"""

        file_list = sorted(os.listdir(sample_dir))
        filtered_filelist = []
        filetypes = set()
        for entry in file_list:
            if entry[-3:] == ".fq" or entry[-6:] == ".fastq":
                filetypes.add("fq")
            elif entry[-6:] == ".fq.gz" or entry[-9:] == ".fastq.gz":
                filetypes.add("fq.gz")
            elif entry[-7:] == ".fq.bz2" or entry[-10:] == ".fastq.bz2":
                filetypes.add("fq.bz2")
            else:
                continue
            filtered_filelist.append("%s/%s" % (sample_dir, entry))

        if input_is_se:
            return filetypes, [], [], filtered_filelist

        single_end_filelist = []
        paired_end_filelist = []

        for entry in filtered_filelist:
            if filename_fragment_to_mark_se_reads in entry:
                single_end_filelist.append(entry)
            else:
                paired_end_filelist.append(entry)

        forward_files = paired_end_filelist[::2]
        reverse_files = paired_end_filelist[1:][::2]
        if len(forward_files) != len(reverse_files):
            raise

        if len(filetypes) > 1:
            print("WARNING: mix of archives of different types and/or uncompressed files")

        return filetypes, forward_files, reverse_files, single_end_filelist

    @staticmethod
    def get_sample_list(samples_directory):
        samples = sorted(os.listdir(samples_directory))
        sample_list = []
        for sample in samples:
            if os.path.isdir("%s/%s" % (samples_directory, sample)):
                sample_list.append(sample)
        return sample_list

    @staticmethod
    def extract_ids_from_file(input_file, output_file=None, header=False, column_separator="\t",
                              comments_prefix="#", column_number=None):
        id_list = IdList()
        id_list.read(input_file, column_separator=column_separator, comments_prefix=comments_prefix,
                     column_number=column_number, header=header)
        if output_file:
            id_list.write(output_file, header=header)
        return id_list

    @staticmethod
    def intersect_ids_from_files(files_with_ids_from_group_a, files_with_ids_from_group_b,
                                 result_file=None, mode="common", case_insensitive=False):
        a = IdSet()
        b = IdSet()

        if mode == "common":
            expression = lambda a, b: a & b
        elif mode == "only_a":
            expression = lambda a, b: a - b
        elif mode == "only_b":
            expression = lambda a, b: b - a
        elif mode == "not_common":
            expression = lambda a, b: a ^ b
        elif mode == "combine":
            expression = lambda a, b: a | b

        #print(files_with_ids_from_group_a)
        for filename in [files_with_ids_from_group_a] if isinstance(files_with_ids_from_group_a, str) else files_with_ids_from_group_a:
            id_set = IdSet()
            id_set.read(filename, comments_prefix="#", expression=(lambda s: s.upper()) if case_insensitive else None)
            a = a | id_set

        for filename in [files_with_ids_from_group_b] if isinstance(files_with_ids_from_group_b, str) else files_with_ids_from_group_b:
            id_set = IdSet()
            id_set.read(filename, comments_prefix="#", expression=(lambda s: s.upper()) if case_insensitive else None)
            b = b | id_set
        result_fd = open(result_file, "w") if result_file else sys.stdout
        if mode != "count":
            final_set = IdSet(expression(a, b))
            final_set.write(result_fd)
        else:
            result_fd.write("Group_A\t%i\nGroup_B\t%i\nCommon\t%i\nOnly_group_A\t%i\nOnly_group_B\t%i\nNot_common\t%i\nAll\t%i\n" %
                            (len(a), len(b), len(a & b), len(a - b), len(b - a), len(a ^ b), len(a | b)))

    @staticmethod
    def intersect_ids(list_of_group_a, list_of_group_b, mode="common"):
        # possible modes: common, only_a, only_b, not_common,  combine, count
        a = IdSet()
        b = IdSet()

        if mode == "common":
            expression = lambda a, b: a & b
        elif mode == "only_a":
            expression = lambda a, b: a - b
        elif mode == "only_b":
            expression = lambda a, b: b - a
        elif mode == "not_common":
            expression = lambda a, b: a ^ b
        elif mode == "combine":
            expression = lambda a, b: a | b

        for id_list in list_of_group_a:
            a = a | IdSet(id_list)

        for id_list in list_of_group_b:
            b = b | IdSet(id_list)

        if mode != "count":
            return IdSet(expression(a, b))
        else:
            return len(a), len(b), len(a & b), len(a - b), len(b - a), len(a ^ b), len(a | b)

    @staticmethod
    def split_by_column(input_file, column_number, separator="\t", header=False, outfile_prefix=None,
                        use_column_value_as_prefix=False, sorted_input=False):
        # column number should start from 0
        # use sorted input to reduce number of simalteniously open files
        header_string = None
        splited_name = input_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd_dict = {}
        previous_value = None
        with open(input_file, "r") as in_fd:
            if header:
                header_string = in_fd.readline()
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] not in out_fd_dict:
                    print (line_str[column_number])
                    if previous_value and sorted_input:
                        out_fd_dict[previous_value].close()
                    suffix = line_str[column_number].replace(" ", "_")
                    out_name = "%s.%s" % (suffix, extension) if use_column_value_as_prefix else "%s.%s.%s" % (out_prefix, suffix, extension)
                    out_fd_dict[line_str[column_number]] = open(out_name, "w")
                    if header:
                        out_fd_dict[line_str[column_number]].write(header_string)

                out_fd_dict[line_str[column_number]].write(line)
                previous_value = line_str[column_number]
        if sorted_input:
            out_fd_dict[previous_value].close()
        else:
            for entry in out_fd_dict:
                out_fd_dict[entry].close()

    @staticmethod
    def get_sample_list(sample_dir, sample_list=None):
        samples = []
        if sample_list:
            return [sample_list] if isinstance(sample_list, str) else sample_list
        else:
            dir_list = os.listdir(sample_dir)
            for directory in dir_list:
                if os.path.isdir("%s/%s" % (sample_dir, directory)):
                    samples.append(directory)

            return samples

    @staticmethod
    def replace_column_value_by_syn(input_file, syn_file, out_file, column=0, comment_prefix=None, separator="\t",
                                    syn_header=False, syn_separator="\t",
                                    syn_key_index=0, syn_value_index=1, syn_comment_prefix=None):
        syn_dict = SynDict(filename=syn_file, header=syn_header, separator=syn_separator, key_index=syn_key_index,
                           value_index=syn_value_index, comments_prefix=syn_comment_prefix)
        if comment_prefix:
            comment_prefix_len = len(comment_prefix)
        line_number = 0
        replaced = 0
        not_replaced = 0
        with open(input_file, "r") as in_fd:
            with open(out_file, "w") as out_fd:
                for line in in_fd:
                    line_number += 1
                    if comment_prefix:
                        if line[0:comment_prefix_len] == comment_prefix:
                            out_fd.write(line)
                            continue
                    line_list = line.strip("\n").split(separator)
                    if len(line_list) < column + 1:
                        sys.stderr.write("WARNING!!! Line %i doesn't have column %i\n" % (line_number, column))
                    if line_list[column] in syn_dict:
                        replaced += 1
                        line_list[column] = syn_dict[line_list[column]]
                    else:
                        not_replaced += 1

                    out_fd.write(separator.join(line_list))
                    out_fd.write("\n")

        sys.stderr.write("Replaced: %i\nNot replaced: %i\n" % (replaced, not_replaced))

    @staticmethod
    def combine_syn_dicts(list_of_syn_dict):
        combined_dict = SynDict()
        for syn_dict in list_of_syn_dict:
            for key in syn_dict:
                if key in combined_dict:
                    combined_dict[key] += syn_dict[key]
                else:
                    combined_dict[key] = syn_dict[key]

        return combined_dict

    def combine_syn_dicts_from_file(self, list_of_syndict_files, output, key_index=0, value_index=1, separator="\t",
                                    values_separator=",", header=False):
        list_of_syn_dicts = []
        for filename in list_of_syndict_files:
            list_of_syn_dicts.append(SynDict(filename=filename, key_index=key_index, value_index=value_index,
                                             separator=separator, values_separator=values_separator, split_values=True,
                                             allow_repeats_of_key=True, header=header))

        merged_dict = self.combine_syn_dicts(list_of_syn_dicts)

        merged_dict.write(output, splited_values=True)

    def add_add_new_column_by_key_column(self, table_file, syn_dict_file, key_column,  output_file, new_column_name=None,
                                         separator='\t', absent_value="."):
        column_syn_dict = SynDict(filename=syn_dict_file, allow_repeats_of_key=True, values_separator="@")
        with open(table_file, "r") as in_fd, open(output_file, "w") as out_fd:
            if new_column_name:
                header_line = in_fd.readline().strip() + "\t%s\n" % new_column_name
                out_fd.write(header_line)
                for line in in_fd:
                    line_list = line.strip().split(separator)
                    if line_list[key_column] in column_syn_dict:
                        print (line_list[key_column])
                        print (column_syn_dict[line_list[key_column]])
                    line_list.append(absent_value if line_list[key_column] not in column_syn_dict else "|".join(column_syn_dict[line_list[key_column]]))
                    out_fd.write(separator.join(line_list) + "\n")

    @staticmethod
    def label_column_in_file(input_file, label, column_index, output_file, column_separator="\t",
                             label_position="first", label_separator="@",
                             comments_prefix="#"):
        with open(input_file, "r") as in_fd:
            with open(output_file, "w") as out_fd:
                for line in in_fd:
                    if line[0] == comments_prefix:
                        out_fd.write(line)
                        continue

                    line_list = line.strip("\n").split(column_separator)
                    if label_position == "first":
                        line_list[column_index] = "%s%s%s" % (label, label_separator, line_list[column_index])
                    elif label_position == "last":
                        line_list[column_index] = "%s%s%s" % (line_list[column_index], label_separator, label)
                    else:
                        raise ValueError("ERROR!!!Unrecognized label position")

                    out_fd.write(column_separator.join(line_list) + "\n")

    def file_line_as_list_generator(self, input_file, comments_prefix="#", separator="\t", buffering=10000000):
        comments_prefix_len = len(comments_prefix)

        with self.metaopen(input_file, "r", buffering=buffering) as in_fd:
            for line in in_fd:
                if line[:comments_prefix_len] == comments_prefix:
                    continue
                yield line.strip().split(separator)

    def file_line_generator(self, input_file, comments_prefix="#"):
        comments_prefix_len = len(comments_prefix)

        with self.metaopen(input_file, "r") as in_fd:
            for line in in_fd:
                if line[:comments_prefix_len] == comments_prefix:
                    continue
                yield line.strip()

    @staticmethod
    def combine_files_with_header(file_list, output_file, header_prefix="#", sorting_options=None):

        #list_of_files = self.make_list_of_path_to_files(file_list)
        unsorted_file = "%s.unsorted.tmp" % output_file

        string = "cat %s > %s" % (file_list[0], unsorted_file if sorting_options else output_file)
        os.system(string)

        for filename in file_list[1:]:
            string = "sed -n '/^[^%s]/p' %s >> %s" % (header_prefix,
                                                      filename,
                                                      unsorted_file if sorting_options else output_file)
            print(string)
            os.system(string)

        if sorting_options:
            sorting_string = "(sed '/^[^#]/Q' %s; sort %s %s) > %s" % (unsorted_file,
                                                                       sorting_options,
                                                                       unsorted_file,
                                                                       output_file)
            print(sorting_string)

            os.system(sorting_string)

    def combine_chunks_with_header(self, chunks_dir, chunks_prefix, output_file,
                                   starting_chunk=None, end_chunk=None, chunk_number_list=None,
                                   chunks_suffix=None, header_prefix="#", sorting_options=None,
                                   separator="_"):

        if chunk_number_list:
            file_list = []
            for chunk_n in chunk_number_list:
                file_list.append("%s/%s%s%s%s" % (chunks_dir, chunks_prefix, separator,
                                                  str(chunk_n), chunks_suffix if chunks_suffix else ""))

        elif starting_chunk or end_chunk:
            if starting_chunk and end_chunk:
                file_list = []
                for chunk_n in range(starting_chunk, end_chunk + 1):
                    file_list.append("%s/%s%s%s%s" % (chunks_dir, chunks_prefix, separator,
                                                      str(chunk_n), chunks_suffix if chunks_suffix else ""))
            else:
                raise ValueError("Either starting or end chunks was not set")
        else:
            file_list = self.make_list_of_path_to_files(chunks_dir)

        self.combine_files_with_header(file_list, output_file,
                                       header_prefix=header_prefix, sorting_options=sorting_options)

    @staticmethod
    def get_filtered_entry_list(entry_list,
                                entry_black_list=None,
                                sort_entries=False,
                                entry_ordered_list=None,
                                entry_white_list=None):
        white_set = set(entry_white_list) if entry_white_list is not None else set()
        black_set = set(entry_black_list) if entry_black_list is not None else set()
        entry_set = set(entry_list)

        if white_set:
            entry_set = entry_set & white_set
        if black_set:
            entry_set = entry_set - black_set

        filtered_entry_list = list(entry_set)
        if sort_entries:
            filtered_entry_list.sort()

        final_entry_list = []

        if entry_ordered_list:
            for entry in entry_ordered_list:
                if entry in filtered_entry_list:
                    final_entry_list.append(entry)
                    filtered_entry_list.remove(entry)
                else:
                    print("WARNING!!!Entry(%s) from order list is absent in list of entries!" % entry)
            return final_entry_list + filtered_entry_list
        else:
            return filtered_entry_list

    def merge_files_by_columns(self, file_list, column_index_list,  output_file, separator="\t", column_names_list=None,
                               comment_prefix="#", header=False):
        comment_prefix_len = len(comment_prefix)
        number_of_files = len(file_list)

        if len(column_index_list) != number_of_files:
            raise ValueError("ERROR!!! Column indexes either were not specified for some files "
                             "or specified too many of them")

        combined_header_list = []

        file_line_dict_list = []

        for file_index in range(0, number_of_files):
            file_line_dict_list.append(OrderedDict())

            with self.metaopen(file_list[file_index], "r") as file_fd:
                if header:
                    header_list = file_fd.readline().strip("\n").split()
                    for i in range(0, len(header_list)):
                        if i in column_index_list[file_index]:
                            continue
                        combined_header_list.append(header_list[i])

                for line in file_fd:
                    if line[0:comment_prefix_len] == comment_prefix:
                        continue
                    line_list = line.strip().split(separator)

                    line_combined_id = ""
                    for column_index in column_index_list[file_index]:
                        line_combined_id += line_list[column_index]
                    file_line_dict_list[-1][line_combined_id] = line_list

        file_line_dict_line_len = []
        for file_dict in file_line_dict_list:
            file_line_dict_line_len.append(len(file_dict[list(file_dict.keys())[0]]))

        with self.metaopen(output_file, "w") as out_fd:
            if header:
                out_fd.write(separator.join(column_names_list + combined_header_list) + "\n")

            for line_id in file_line_dict_list[0]:
                for j in range(1, len(file_line_dict_list)):
                    if line_id not in file_line_dict_list[j]:
                        break
                else:
                    output_line_list = []
                    for column_index in column_index_list[0]:
                        output_line_list.append(file_line_dict_list[0][line_id][column_index])

                    for file_dict_index in range(0, len(file_line_dict_list)):
                        for column_index in range(0, file_line_dict_line_len[file_dict_index]):
                            if column_index in column_index_list[file_dict_index]:
                                continue
                            output_line_list.append(file_line_dict_list[file_dict_index][line_id][column_index])

                    out_fd.write(separator.join(output_line_list) + "\n")

    def rename_ncbi_genome_files(self, genome_dir):
        species_list = os.listdir(genome_dir)

        assembly_suffix = "_genomic.fna"
        gff_suffix = "_genomic.gff"
        protein_suffix = "_protein.faa"
        rna_suffix = "_rna.fna"
        tab_suffix = ".txt"
        for species in species_list:
            species_path = "%s/%s" % (genome_dir, species)
            source_list = os.listdir(species_path)
            for source in source_list:
                species_source_path = "%s/%s" % (species_path, source)
                if os.path.isdir(species_source_path):
                    assembly_list = os.listdir(species_source_path)
                    for assembly in assembly_list:
                        species_source_assembly_path = "%s/%s" % (species_source_path, assembly)
                        if os.path.isdir(species_source_assembly_path):
                            os.system("gunzip %s/*" % species_source_assembly_path)
                            file_list = os.listdir(species_source_assembly_path)
                            for filename in file_list:
                                filename_path = "%s/%s" % (species_source_assembly_path, filename)
                                if os.path.isfile(filename_path):
                                    if assembly_suffix in filename:
                                        shutil.move(filename_path, "%s/%s.fasta" % (species_source_assembly_path, species))
                                    elif gff_suffix in filename:
                                        shutil.move(filename_path, "%s/%s.gff" % (species_source_assembly_path, species))
                                    elif protein_suffix in filename:
                                        shutil.move(filename_path, "%s/%s.pep" % (species_source_assembly_path, species))
                                    elif rna_suffix in filename:
                                        shutil.move(filename_path, "%s/%s.transcript" % (species_source_assembly_path, species))
                                    elif tab_suffix in filename:
                                        shutil.move(filename_path, "%s/%s.tab" % (species_source_assembly_path, species))



    @staticmethod
    def string2float(string):
        try:
            retval = float(string)
        except ValueError:
            retval = string
        return retval

    @staticmethod
    def string2int(string):
        return int(string) if string.isdigit() else string

    def natural_keys_int(self, string):
        """
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        """

        return [self.string2int(c) for c in re.split('(\d+)', string)]

    def natural_keys_float(self, string):
        """
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        float regex comes from https://stackoverflow.com/a/12643073/190597
        """
        return [self.string2float(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', string)]

    def get_column_value_set_from_file(self, input_file, column_number, output_file=None,
                                       separator="\t", comments_prefix="#", verbose=False):

        column_value_set = IdSet([line_list[column_number] for line_list in self.file_line_as_list_generator(input_file,
                                                                                                             separator=separator,
                                                                                                             comments_prefix=comments_prefix)])
        if output_file:
            column_value_set.write(output_file)

        if verbose:
            print("#Column %i (0-based) contains %i different values" % (column_number, len(column_value_set)))

        return column_value_set

    def count_column_values_from_file(self, input_file, column_number, output_file=None,
                                      separator="\t", comments_prefix="#", verbose=False):

        column_value_dict = SynDict()

        for line_list in self.file_line_as_list_generator(input_file,
                                                          separator=separator,
                                                          comments_prefix=comments_prefix):

            if line_list[column_number] in column_value_dict:
                column_value_dict[line_list[column_number]] += 1
            else:
                column_value_dict[line_list[column_number]] = 1

        if output_file:
            column_value_dict.write(output_file)

        return column_value_dict

filetypes_dict = {"fasta": [".fa", ".fasta", ".fa", ".pep", ".cds"],
                  "fastq": [".fastq", ".fq"],
                  "genbank": [".gb", ".genbank"],
                  "newick": [".nwk"]}


def safe_mkdir(dirname):
    try:
        os.mkdir(dirname)
    except OSError:
        pass


def detect_filetype_by_extension(filename, filetypes_dict=filetypes_dict):
    directory, prefix, extension = split_filename(filename)
    for filetype in filetypes_dict:
        if extension in filetypes_dict[filetype]:
            return filetype
    return None


def check_path(path_to_check):
    #print (path_to_check)
    #returns path with / at end or blank path
    if path_to_check != "":
        if path_to_check[-1] != "/":
            return path_to_check + "/"
    return path_to_check


def split_filename(filepath):
    directory, basename = os.path.split(filepath)
    prefix, extension = os.path.splitext(basename)
    return directory, prefix, extension


def make_list_of_path_to_files(list_of_dirs_and_files, expression=None):

    paths_list = []
    for entry in list_of_dirs_and_files:
        #print entry
        if os.path.isdir(entry):
            files_in_dir = sorted(filter(expression, os.listdir(entry)) if expression else os.listdir(entry))
            for filename in files_in_dir:
                paths_list.append("%s%s" % (check_path(entry), filename))
        elif os.path.exists(entry):
            if expression:
                if expression(entry):
                    paths_list.append(entry)
            else:
                paths_list.append(entry)
        else:
            print("%s does not exist" % entry)

    return paths_list


def read_synonyms_dict(filename, header=False, separator="\t",
                       split_values=False, values_separator=",", key_index=0, value_index=1):
    # reads synonyms from file
    synonyms_dict = OrderedDict()
    with open(filename, "r") as in_fd:
        if header:
            header_str = in_fd.readline().strip()
        for line in in_fd:
            tmp = line.strip().split(separator)
            #print line
            key, value = tmp[key_index], tmp[value_index]
            if split_values:
                value = value.split(values_separator)
            synonyms_dict[key] = value
    return synonyms_dict


def read_ids(filename, header=False, close_after_if_file_object=False):
    #reads ids from file or file object with one id per line
    id_list = []

    in_fd = filename if isinstance(filename, file) else open(filename, "r")

    if header:
        header_str = in_fd.readline().strip()
    for line in in_fd:
        id_list.append(line.strip())
    if (not isinstance(filename, file)) or close_after_if_file_object:
        in_fd.close()

    return id_list


def read_tsv_as_rows_list(filename, header=False, header_prefix="#", separator="\t"):
    tsv_list = []
    with open(filename, "r") as tsv_fd:
        if header:
            header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
        for line in tsv_fd:
            tsv_list.append(line.strip().split(separator))
    return header_list, tsv_list if header else tsv_list


def read_tsv_as_columns_dict(filename, header_prefix="#", separator="\t"):
    tsv_dict = OrderedDict()
    with open(filename, "r") as tsv_fd:
        header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
        number_of_columns = len(header_list)
        for column_name in header_list:
            tsv_dict[column_name] = []
        for line in tsv_fd:
            tmp_line = line.strip().split(separator)
            for i in range(0, number_of_columns):
                tsv_dict[header_list[i]].append(tmp_line[i])
    return tsv_dict

def tsv_extract_by_column_value(tsv_file, column_number, column_value, separator="\t",
                                header=False, outfile_prefix=None):
    # column number should start from 0
    # column_value should be string or list of strings

    splited_name = tsv_file.split(".")
    extension = splited_name[-1] if len(splited_name) > 1 else ""
    out_prefix = outfile_prefix if outfile_prefix is not None \
        else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
    out_fd = open("%s.%s" % (out_prefix, extension), "w")

    column_value_list = column_value if isinstance(column_value, Iterable) else []
    with open(tsv_file, "r") as in_fd:
        if header:
            out_fd.write(in_fd.readline())
        for line in in_fd:
            line_str = line.strip().split(separator)
            if line_str[column_number] in column_value_list:
                out_fd.write(line)
    out_fd.close()

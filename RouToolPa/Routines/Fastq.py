__author__ = 'mahajrod'
import os
import datetime
from collections import OrderedDict
from Bio.Seq import Seq
from RouToolPa.Collections.General import TwoLvlDict, IdSet
from RouToolPa.Routines.File import FileRoutines


class FastQRoutines(FileRoutines):

    def __init__(self):
        pass

    def reverse_complement(self, in_file, out_file):
        with self.metaopen(in_file, "r") as in_fd:
            with self.metaopen(out_file, "w") as out_fd:
                for line in in_fd:
                    out_fd.write(line)
                    out_fd.write(str(Seq(in_fd.next().strip()).reverse_complement()))
                    out_fd.write("\n")
                    out_fd.write(in_fd.next())
                    out_fd.write(in_fd.next().strip()[::-1])
                    out_fd.write("\n")

    @staticmethod
    def parse_illumina_name(string):
        name_list = string[1:].split(" ")
        return tuple(name_list[0].split(":"))

    def remove_tiles_from_fastq(self, forward_reads, black_list_forward_tiles_list,
                                reverse_reads, black_list_reverse_tiles_list, output_prefix):

        filtered_paired_forward_pe = "%s.ok.pe_1.fastq" % output_prefix
        filtered_forward_se = "%s.ok.forward.se.fastq" % output_prefix
        filtered_out_forward_se = "%s.bad.forward.fastq" % output_prefix

        filtered_paired_reverse_pe = "%s.ok.pe_2.fastq" % output_prefix
        filtered_reverse_se = "%s.ok.reverse.se.fastq" % output_prefix
        filtered_out_reverse_se = "%s.bad.reverse.fastq" % output_prefix

        forward_input_fd = self.metaopen(forward_reads, "r")
        reverse_input_fd = self.metaopen(reverse_reads, "r")

        filtered_paired_forward_pe_fd = self.metaopen(filtered_paired_forward_pe, "w")
        filtered_forward_se_fd = self.metaopen(filtered_forward_se, "w")
        filtered_out_forward_se_fd = self.metaopen(filtered_out_forward_se, "w")

        filtered_paired_reverse_pe_fd = self.metaopen(filtered_paired_reverse_pe, "w")
        filtered_reverse_se_fd = self.metaopen(filtered_reverse_se, "w")
        filtered_out_reverse_se_fd = self.metaopen(filtered_out_reverse_se, "w")

        for line in forward_input_fd:
            name_list = line.strip()[1:].split(":")
            #instrument_id = name_list[0]
            #run_id = name_list[1]
            #flowcell_id = name_list[2]
            #lane_id = name_list[3]
            tile_id = name_list[4]

            if (tile_id in black_list_forward_tiles_list) and (tile_id in black_list_reverse_tiles_list):
                filtered_out_forward_se_fd.write(line)
                for i in range(0, 3):
                    filtered_out_forward_se_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_out_reverse_se_fd.write(reverse_input_fd.next())

            elif (tile_id in black_list_forward_tiles_list) and (not(tile_id in black_list_reverse_tiles_list)):
                filtered_out_forward_se_fd.write(line)
                for i in range(0, 3):
                    filtered_out_forward_se_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_reverse_se_fd.write(reverse_input_fd.next())

            elif (not (tile_id in black_list_forward_tiles_list)) and (tile_id in black_list_reverse_tiles_list):
                filtered_forward_se_fd.write(line)
                for i in range(0, 3):
                    filtered_forward_se_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_out_reverse_se_fd.write(reverse_input_fd.next())

            else:
                filtered_paired_forward_pe_fd.write(line)
                for i in range(0, 3):
                    filtered_paired_forward_pe_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_paired_reverse_pe_fd.write(reverse_input_fd.next())

        filtered_paired_forward_pe_fd.close()
        filtered_forward_se_fd.close()
        filtered_out_forward_se_fd.close()

        filtered_paired_reverse_pe_fd.close()
        filtered_reverse_se_fd.close()
        filtered_out_reverse_se_fd.close()

    def combine_fastq_files(self, samples_directory, sample, output_directory,
                            use_links_if_merge_not_necessary=True, input_is_se=False):
        sample_dir = "%s/%s/" % (samples_directory, sample)

        filetypes, forward_files, reverse_files, se_files = self.make_lists_forward_and_reverse_files(sample_dir,
                                                                                                      input_is_se=input_is_se)

        uncompresed = True
        if len(filetypes) == 1:
            if ("fq.gz" in filetypes) or ("fastq.gz" in filetypes):
                command = "zcat"
                uncompresed = False
            elif ("fq.bz2" in filetypes) or ("fastq.bz2" in filetypes):
                command = "bzcat"
                uncompresed = False
            else:
                command = "cat"

            merged_forward = "%s/%s_1.fq" % (output_directory, sample)
            merged_reverse = "%s/%s_2.fq" % (output_directory, sample)
            merged_se = "%s/%s.se.fq" % (output_directory, sample)

            if use_links_if_merge_not_necessary and (len(forward_files) == 1) and (len(reverse_files) == 1) and (uncompresed == True):
                forward_string = "ln -s %s %s" % (forward_files[0], merged_forward)
                reverse_string = "ln -s %s %s" % (reverse_files[0], merged_reverse)

                print("Executing: %s" % forward_string)
                os.system(forward_string)
                print("Executing: %s" % reverse_string)
                os.system(reverse_string)

                if len(se_files) == 1:
                    os.system("ln -s %s %s" % (se_files[0], merged_se))
                    return merged_forward, merged_reverse, merged_se
                elif len(se_files) > 0:
                    os.system("%s %s > %s" % (command, " ".join(se_files), merged_se))
                    return merged_forward, merged_reverse, merged_se
                else:
                    return merged_forward, merged_reverse, None
            else:
                if (len(forward_files) > 0) and (len(reverse_files) > 0):
                    os.system("%s %s > %s" % (command, " ".join(forward_files), merged_forward))
                    os.system("%s %s > %s" % (command, " ".join(reverse_files), merged_reverse))
                    if len(se_files) > 0:
                        os.system("%s %s > %s" % (command, " ".join(se_files), merged_se))
                        return merged_forward, merged_reverse, merged_se
                    else:
                        return merged_forward, merged_reverse, None
                if len(se_files) > 0:
                    os.system("%s %s > %s" % (command, " ".join(se_files), merged_se))
                    return None, None, merged_se
                else:
                    raise IOError("No input files were found")
        else:
            raise IOError("Extracting from mix of archives in not implemented yet")

    @staticmethod
    def filter_se_by_length(input_file, filtered_file, filtered_out_file, min_len=None, max_len=None):

        if min_len and max_len:
            def expression(line):
                return min_len <= len(line) <= max_len
        elif min_len:
            def expression(line):
                return min_len <= len(line)
        elif max_len:
            def expression(line):
                return len(line) <= max_len
        else:
            raise ValueError("Both minimum and maximum thresholds for read length were not set")

        with open(input_file, "r") as in_fd:
            with open(filtered_file, "w") as filtered_fd:
                with open(filtered_out_file, "w") as filtered_out_fd:
                    for read_name in in_fd:
                        read = in_fd.next()
                        #print len(read.strip())
                        #print expression(read.strip())
                        if expression(read.strip()):
                            filtered_fd.write(read_name)
                            filtered_fd.write(read)
                            filtered_fd.write(in_fd.next())
                            filtered_fd.write(in_fd.next())
                        else:
                            filtered_out_fd.write(read_name)
                            filtered_out_fd.write(read)
                            filtered_out_fd.write(in_fd.next())
                            filtered_out_fd.write(in_fd.next())

    def split_illumina_fastq_by_lanes(self, input_fastq, output_dir, output_prefix=None, output_suffix=".fastq"):
        out_fd_dict = OrderedDict()

        with self.metaopen(input_fastq, "r") as in_fd:
            for line in in_fd:
                tmp = line.split(":")
                lane_id = ".".join(tmp[:4])[1:]

                if lane_id not in out_fd_dict:
                    out_fd_dict[lane_id] = open("%s/%s%s" % (output_dir,
                                                             ((output_prefix + ".") if output_prefix else "") + lane_id,
                                                             output_suffix), "w")

                out_fd_dict[lane_id].write(line)

                for i in range(0, 3):
                    out_fd_dict[lane_id].write(in_fd.next())

        for fd in out_fd_dict:
            out_fd_dict[fd].close()

    def find_tiles(self, fastq_file, output_file):
        with self.metaopen(fastq_file, "r") as fastq_fd, (output_file if isinstance(output_file, file) else open(output_file, "w")) as out_fd:
            out_fd.write("#Machine\tRun\tFlowcellID\tLane\tTile\n")
            tile_set = set()

            for line in fastq_fd:
                read_name_tuple = self.parse_illumina_name(line)[:5]
                tile_set.add(read_name_tuple)
                for i in 0, 1, 2:
                    fastq_fd.next()

            tile_list = list(tile_set)
            tile_list.sort()

            for tile in tile_list:
                out_fd.write("\t".join(tile))
                out_fd.write("\n")

    def count_reads_in_tiles(self, fastq_file, output_file):
        with self.metaopen(fastq_file, "r") as fastq_fd, (output_file if isinstance(output_file, file) else open(output_file, "w")) as out_fd:
            out_fd.write("#Machine\tRun\tFlowcellID\tLane\tTile\tReadNumber\n")
            tile_dict = OrderedDict

            for line in fastq_fd:
                tile_name = "\t".join(self.parse_illumina_name(line)[:5])
                if tile_name in tile_dict:
                    tile_dict[tile_name] += 1
                else:
                    tile_dict[tile_name] = 1

                for i in 0, 1, 2:
                    fastq_fd.next()

            for tile_name in sorted(tile_dict.keys()):
                out_fd.write(tile_name)
                out_fd.write("\t")
                out_fd.write(str(tile_dict[tile_name]))
                out_fd.write("\n")

    def count_reads_and_bases(self, fastq_file_list, stat_file=None):

        fastq_list = [fastq_file_list] if isinstance(fastq_file_list, str) else fastq_file_list

        counts = TwoLvlDict()

        for fastq_file in fastq_list:
            counts[fastq_file] = OrderedDict()
            counts[fastq_file]["Reads"] = 0
            counts[fastq_file]["Bases"] = 0

        for fastq_file in fastq_list:
            with self.metaopen(fastq_file, "r") as fastq_fd:
                for line in fastq_fd:
                    counts[fastq_file]["Bases"] += len(fastq_fd.next())
                    counts[fastq_file]["Reads"] += 1
                    fastq_fd.next()
                    fastq_fd.next()

                # to take into account "\n" at the end of each line
                counts[fastq_file]["Bases"] = counts[fastq_file]["Bases"] - counts[fastq_file]["Reads"]

        counts.write()

        if stat_file:
            counts.write(stat_file)

    def extract_10x_barcodes(self, forward_file_list, reverse_file_list, index_file_list, barcode_file, output_prefix,
                             buffering=1000000000, read_index_length=16, linker_length=6, min_forward_read_len=50):
        index_end = linker_start = read_index_length
        linker_end = read_start = read_index_length + linker_length

        service_seq_length = read_index_length + linker_length
        min_forward_seq_length = service_seq_length + min_forward_read_len

        with self.metaopen(barcode_file, "r", buffering=buffering) as in_fd:
            barcode_set = set(map(lambda s: s.strip(), in_fd.readlines()))

        output_dict = {
                       "good": {
                                "forward": "%s.good_1.fastq" % output_prefix,
                                "reverse": "%s.good_2.fastq" % output_prefix,
                                "index": "%s.good_I.fastq" % output_prefix,
                                "linker": "%s.good_L.fastq" % output_prefix,
                                },
                       "bad":  {
                                "forward": "%s.bad_1.fastq" % output_prefix,
                                "reverse": "%s.bad_2.fastq" % output_prefix,
                                "index": "%s.bad_I.fastq" % output_prefix,
                                "linker": "%s.bad_L.fastq" % output_prefix,
                                },
                       "short": {
                                 "forward": "%s.short_1.fastq" % output_prefix,
                                 "reverse": "%s.short_2.fastq" % output_prefix,
                                 "index": "%s.short_I.fastq" % output_prefix,
                                }
                       }
        output_dict_fd = {}
        input_dict_fd = {}

        for key, filelist in zip(["forward", "reverse", "index"], [forward_file_list, reverse_file_list, index_file_list]):
            input_dict_fd[key] = list(map(lambda s: self.metaopen(s, "r", buffering=buffering), filelist))

        for quality in "good", "bad", "short":
            output_dict_fd[quality] = {}
            for seq_type in ["forward", "reverse", "index", "linker"]:
                output_dict_fd[quality][seq_type] = self.metaopen(output_dict[quality][seq_type], "w", buffering=buffering)

        counter_dict = OrderedDict({"handled": 0,
                                    "good": 0,
                                    "bad": 0})
        print("%s\tStarting...\n" % str(datetime.datetime.now()))
        for forward_fd, reverse_fd, index_fd in zip(input_dict_fd["forward"], input_dict_fd["reverse"], input_dict_fd["index"]):

            for line in forward_fd:
                counter_dict["handled"] += 1
                if counter_dict["handled"] % 1000000 == 0:
                    print("%s\tHandled %i read pairs\n" % (str(datetime.datetime.now()), counter_dict["handled"]))

                read_seq = forward_fd.readline()
                delimiter = forward_fd.readline()
                read_qual = forward_fd.readline()

                if len(read_seq) - 1 < min_forward_seq_length:
                    label = "short"

                    for i in (0, 1, 2, 3):
                        output_dict_fd[label]["index"].write(index_fd.readline())
                else:
                    index = read_seq[:index_end]
                    label = "good" if index in barcode_set else "bad"

                    output_dict_fd[label]["linker"].write(line)
                    output_dict_fd[label]["linker"].write(read_seq[linker_start:linker_end] + "\n")
                    output_dict_fd[label]["linker"].write(delimiter)
                    output_dict_fd[label]["linker"].write(read_qual[linker_start:linker_end] + "\n")

                    output_dict_fd[label]["index"].write(index_fd.readline())
                    output_dict_fd[label]["index"].write(index_fd.readline()[:-1] + read_seq[:index_end] + "\n")
                    output_dict_fd[label]["index"].write(index_fd.readline())
                    output_dict_fd[label]["index"].write(index_fd.readline()[:-1] + read_qual[:index_end] + "\n")

                counter_dict[label] += 1

                output_dict_fd[label]["forward"].write(line)
                output_dict_fd[label]["forward"].write(read_seq[read_start:])
                output_dict_fd[label]["forward"].write(delimiter)
                output_dict_fd[label]["forward"].write(read_qual[read_start:])

                for i in (0, 1, 2, 3):
                    output_dict_fd[label]["reverse"].write(reverse_fd.readline())

        for seq_type in input_dict_fd:
            for fd in input_dict_fd[seq_type]:
                fd.close()

        for quality in output_dict_fd:
            for seq_type in output_dict_fd[quality]:
                output_dict_fd[quality][seq_type].close()

        for key in counter_dict:
            print("%s read pairs: %i\n" % (key, counter_dict[key]))

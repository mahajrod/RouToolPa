__author__ = 'mahajrod'

from pathlib import Path
from RouToolPa.Routines.Sequence import SequenceRoutines
from RouToolPa.Routines.Fastq import FastQRoutines


class FilteringRoutines(SequenceRoutines, FastQRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)
        FastQRoutines.__init__(self)

    def extract_reads_by_kraken_report(self, forward_fastq, kraken_output, taxon_id_list, output_prefix,
                                       reverse_fastq=None, gzip=True, check_read_id=False):

        if reverse_fastq is not None:
            forward_output_fd = self.metaopen(output_prefix + "_1.fastq" + (".gz" if gzip else ""), "w")
            reverse_output_fd = self.metaopen(output_prefix + "_2.fastq" + (".gz" if gzip else ""), "w")

        else:
            forward_output_fd = self.metaopen(output_prefix + ".fastq" + (".gz" if gzip else ""), "w")

        with self.metaopen(forward_fastq, "r") as forward_input_fd, \
             (self.metaopen(reverse_fastq, "r") if reverse_fastq else None) as reverse_input_fd, \
             self.metaopen(kraken_output, "r") as kraken_input_fd:
            print(taxon_id_list)
            for kraken_line in kraken_input_fd:
                print(kraken_line.split()[2])
                if kraken_line.split()[2] in taxon_id_list:

                    if check_read_id:
                        forward_first_line = forward_input_fd.readline()
                        forward_id = forward_first_line.split()[0][1:]
                        if forward_id != kraken_line[1]:
                            print(kraken_line)
                            print(forward_first_line)
                            raise ValueError("ERROR!!! Read ids in KRAKEN output and read files doesn't match! ")
                        else:
                            print("AAAA")
                            forward_output_fd.write(forward_first_line)
                            for i in 0, 1, 2:
                                forward_output_fd.write(forward_input_fd.readline())
                            if reverse_fastq is not None:
                                for i in 0, 1, 2, 3:
                                    reverse_output_fd.write(reverse_input_fd.readline())
                    else:
                        for i in 0, 1, 2, 3:
                            forward_output_fd.write(forward_input_fd.readline())
                        if reverse_fastq is not None:
                            for i in 0, 1, 2, 3:
                                reverse_output_fd.write(reverse_input_fd.readline())
                else:
                    print("BBBB")
                    for i in 0, 1, 2, 3:
                        forward_input_fd.readline()
                    if reverse_fastq is not None:
                        for i in 0, 1, 2, 3:
                            reverse_input_fd.readline()

        if reverse_fastq is not None:
            forward_output_fd.close()
            reverse_output_fd.close()
        else:
            forward_output_fd.close()

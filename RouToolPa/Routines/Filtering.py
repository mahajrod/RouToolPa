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
            forward_output_fd = self.metaopen(output_prefix + "_1.fastq" + (".gz" if gzip else ""), "w", buffering=100000000)
            reverse_output_fd = self.metaopen(output_prefix + "_2.fastq" + (".gz" if gzip else ""), "w", buffering=100000000)

        else:
            forward_output_fd = self.metaopen(output_prefix + ".fastq" + (".gz" if gzip else ""), "w", buffering=100000000)

        with self.metaopen(forward_fastq, "r", buffering=100000000) as forward_input_fd, \
             (self.metaopen(reverse_fastq, "r", buffering=100000000) if reverse_fastq else None) as reverse_input_fd, \
             self.metaopen(kraken_output, "r", buffering=100000000) as kraken_input_fd:
            counter = 0
            for kraken_line in kraken_input_fd:
                kraken_line_list = kraken_line.split()
                if kraken_line_list[2] in taxon_id_list:
                    if check_read_id:
                        forward_first_line = forward_input_fd.readline()
                        forward_id = forward_first_line.split()[0][1:]
                        if forward_id != kraken_line_list[1]:
                            raise ValueError("ERROR!!! Read ids in KRAKEN output and read files doesn't match! ")
                        else:
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
                    for i in 0, 1, 2, 3:
                        forward_input_fd.readline()
                    if reverse_fastq is not None:
                        for i in 0, 1, 2, 3:
                            reverse_input_fd.readline()

                counter += 1
                if counter % 1000000 == 0:
                    print("Processed {} reads".format(counter))

        if reverse_fastq is not None:
            forward_output_fd.close()
            reverse_output_fd.close()
        else:
            forward_output_fd.close()

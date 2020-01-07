
from RouToolPa.Tools.Abstract import Tool


class EMA(Tool):
    def __init__(self, path="", max_threads=4, max_memory="100G", max_per_thread_memory="5G"):
        Tool.__init__(self, "ema", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

    def convert_ema_read_format_to_fastq(self, input, output_prefix, bufferig=1000000000):
        forward_output = "%s_1.fastq" % output_prefix
        reverse_output = "%s_2.fastq" % output_prefix
        index_output = "%s.index" % output_prefix
        with self.metaopen(input, "r", buffering=bufferig) as in_fd, \
             self.metaopen(forward_output, "w", buffering=bufferig) as f_fd, \
             self.metaopen(reverse_output, "w", buffering=bufferig) as r_fd, \
             self.metaopen(index_output, "w", buffering=bufferig) as i_fd:
            for line in in_fd:
                line_list = line.strip().split()

                f_fd.write(line_list[1])
                f_fd.write("\n")
                f_fd.write(line_list[2])
                f_fd.write("\n")
                f_fd.write("+\n")
                f_fd.write(line_list[3])
                f_fd.write("\n")

                r_fd.write(line_list[1])
                r_fd.write("\n")
                r_fd.write(line_list[4])
                r_fd.write("\n")
                r_fd.write("+\n")
                r_fd.write(line_list[5])
                r_fd.write("\n")

                i_fd.write(line_list[0])
                i_fd.write("\n")



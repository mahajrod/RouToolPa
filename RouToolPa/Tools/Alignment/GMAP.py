
from RouToolPa.Tools.Abstract import Tool


class GMAP(Tool):
    def __init__(self, path="", max_threads=4, max_memory="10G", max_per_thread_memory="5G"):
        Tool.__init__(self, "gmap", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

    def index(self, fasta_file, prefix=None):
        pass

    def align(self):
        pass

    def filter_gff(self, in_file, out_file, min_coverage=None, min_identity=None, gene_feature="gene",
                   transcript_feature="mRNA", exon_feature="exon", cds_feature="CDS", verbose=False):
        if min_coverage is not None and min_identity is not None:
            def check(attributes_dict):
                if (float(attributes_dict["coverage"]) >= min_coverage) and (float(attributes_dict["identity"]) >= min_identity):
                    return True
                return False
        elif min_coverage is not None:
            def check(attributes_dict):
                if float(attributes_dict["coverage"]) >= min_coverage:
                    return True
                return False
        elif min_identity is not None:
            def check(attributes_dict):
                if float(attributes_dict["identity"]) >= min_identity:
                    return True
                return False
        else:
            raise ValueError("ERROR!!! Neither minimum coverage nor minimum identity was set!")

        line_list_gen = self.file_line_as_list_generator(in_file)
        total_transcripts = 0
        retained_transcripts = 0
        with open(out_file, "w") as out_fd:
            write_flag = False
            for line_list in line_list_gen:
                if line_list[2] == transcript_feature:
                    total_transcripts += 1
                    attributtes_dict = dict(map(lambda s: s.split("="), line_list[8].split(";")))
                    if check(attributtes_dict):
                        retained_transcripts += 1
                        out_fd.write("\t".join(line_list))
                        out_fd.write("\n")
                        write_flag = True
                    else:
                        write_flag = False

                elif ((line_list[2] == cds_feature) or (line_list[2] == exon_feature)) and write_flag:
                    out_fd.write("\t".join(line_list))
                    out_fd.write("\n")

        if verbose:
            print("Total transcripts: %i\n Retained transcripts: %i\n" % (total_transcripts, retained_transcripts))

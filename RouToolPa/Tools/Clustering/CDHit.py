#!/usr/bin/env python
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Parsers.CDHIT import CollectionCDHIT


class CDHit(Tool):
    """
    T
    """

    def __init__(self, path="", max_threads=4, max_memmory=0):
        Tool.__init__(self, "cd-hit", path=path, max_threads=max_threads, max_memory=max_memmory)

    def parse_options(self, input_file, output_prefix, sequence_identity=None,
                      shorter_sequence_fraction_cutoff=None,
                      shorter_sequence_length_cutoff=None,
                      longer_sequence_alignment_coverage_fraction_cutoff=None,
                      longer_sequence_alignment_coverage_length_cutoff=None,
                      shorter_sequence_alignment_coverage_fraction_cutoff=None,
                      shorter_sequence_alignment_coverage_length_cutoff=None,
                      longer_sequence_unmatched_fraction_cutoff=None,
                      shorter_sequence_unmatched_fraction_cutoff=None,
                      sort_by_cluster_size=True):
        """
        NOT ALL OPTIONS WERE IMPLEMENTED
        :param input_file:
        :param output_prefix:
        :param sequence_identity:
        :param shorter_sequence_fraction_cutoff:
        :param shorter_sequence_length_cutoff:
        :param longer_sequence_alignment_coverage_fraction_cutoff:
        :param longer_sequence_alignment_coverage_length_cutoff:
        :param shorter_sequence_alignment_coverage_fraction_cutoff:
        :param shorter_sequence_alignment_coverage_length_cutoff:
        :param longer_sequence_unmatched_fraction_cutoff:
        :param shorter_sequence_unmatched_fraction_cutoff:
        :param sort_by_cluster_size:
        :return:
        """

        options = "-T %i" % self.threads
        options += " -M %s" % str(self.max_memory)
        options += " -i %s" % input_file
        options += " -o %s" % output_prefix
        options += " -sc 1" if sort_by_cluster_size else ""
        options += " -c %f" % sequence_identity if sequence_identity else "" # default 0.9

        options += " -s %f" % shorter_sequence_fraction_cutoff if shorter_sequence_fraction_cutoff else "" # default 0.0, i.e no limit
        options += " -S %i" % shorter_sequence_length_cutoff if shorter_sequence_length_cutoff else ""  # default 999999, i.e no limit

        options += " -aL %f" % longer_sequence_alignment_coverage_fraction_cutoff if longer_sequence_alignment_coverage_fraction_cutoff else ""  # default 0.0, i.e no limit
        options += " -AL %i" % longer_sequence_alignment_coverage_length_cutoff if longer_sequence_alignment_coverage_length_cutoff else ""  # default 99999999, i.e no limit

        options += " -aS %f" % shorter_sequence_alignment_coverage_fraction_cutoff if shorter_sequence_alignment_coverage_fraction_cutoff else ""  # default 0.0, i.e no limit
        options += " -AS %i" % shorter_sequence_alignment_coverage_length_cutoff if shorter_sequence_alignment_coverage_length_cutoff else ""  # default 99999999, i.e no limit

        options += " -uS %f" % shorter_sequence_unmatched_fraction_cutoff if shorter_sequence_unmatched_fraction_cutoff else ""  # default 1.0, i.e no limit
        options += " -uL %i" % longer_sequence_unmatched_fraction_cutoff if longer_sequence_unmatched_fraction_cutoff else ""  # default 1.0, i.e no limit

        return options

    def find_duplicates(self, input_file, output_prefix):
        cluster_file = "%s.clstr" % output_prefix

        options = self.parse_options(input_file, output_prefix,
                                     sequence_identity=1.0,
                                     shorter_sequence_fraction_cutoff=1.0,
                                     shorter_sequence_length_cutoff=0,
                                     longer_sequence_alignment_coverage_fraction_cutoff=1.0,
                                     longer_sequence_alignment_coverage_length_cutoff=0,
                                     shorter_sequence_alignment_coverage_fraction_cutoff=1.0,
                                     shorter_sequence_alignment_coverage_length_cutoff=0,
                                     longer_sequence_unmatched_fraction_cutoff=0.0,
                                     shorter_sequence_unmatched_fraction_cutoff=0.0,
                                     sort_by_cluster_size=True)

        self.execute(options=options)

        cdhit_collection = CollectionCDHIT(input_file=cluster_file)

        for ext in "fam", "tab":
            cdhit_collection.write("%s.%s" % (output_prefix, ext), format=ext)



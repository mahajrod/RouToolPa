#!/usr/bin/env python
import os
from Bio import SeqIO
from RouToolPa.Collections.General import IdList, SynDict
from RouToolPa.Routines.SequenceCluster import SequenceClusterRoutines
from RouToolPa.GeneralRoutines import FileRoutines


class TreeFamRoutines(SequenceClusterRoutines):
    def __init__(self):

        pass

    @staticmethod
    def extract_proteins_from_selected_families(families_id_file, fam_file, pep_file,
                                                output_dir="./", pep_format="fasta",
                                                out_prefix=None, create_dir_for_each_family=False):
        from RouToolPa.Routines import SequenceRoutines

        fam_id_list = IdList()
        fam_dict = SynDict()
        #print(pep_file)
        FileRoutines.safe_mkdir(output_dir)
        out_dir = FileRoutines.check_path(output_dir)
        create_directory_for_each_family = True if out_prefix else create_dir_for_each_family
        if families_id_file:
            fam_id_list.read(families_id_file)
        fam_dict.read(fam_file, split_values=True, values_separator=",")
        protein_dict = SeqIO.index_db("tmp.idx", pep_file, format=pep_format)

        for fam_id in fam_id_list if families_id_file else fam_dict:
            if fam_id in fam_dict:
                if create_directory_for_each_family:
                    fam_dir = "%s%s/" % (out_dir, fam_id)
                    FileRoutines.safe_mkdir(fam_dir)
                    out_file = "%s%s.pep" % (fam_dir, out_prefix if out_prefix else fam_id)
                else:
                    out_file = "%s/%s.pep" % (out_dir, out_prefix if out_prefix else fam_id)

                SeqIO.write(SequenceRoutines.record_by_id_generator(protein_dict, fam_dict[fam_id], verbose=True),
                            out_file, format=pep_format)
            else:
                print("%s was not found" % fam_id)

        os.remove("tmp.idx")

    @staticmethod
    def add_length_to_fam_file(fam_file, len_file, out_file, close_after_if_file_object=False):
        fam_dict = SynDict()
        fam_dict.read(fam_file, split_values=True, comments_prefix="#")
        len_dict = SynDict()
        len_dict.read(len_file, comments_prefix="#")

        out_fd = out_file if isinstance(out_file, file) else open(out_file, "r")

        for family in fam_dict:
            len_list = []
            for member in fam_dict[family]:
                len_list.append(None if member not in len_dict else len_dict[member])

            out_fd.write("%s\t%s\t%s\n" % (family, ",".join(fam_dict[family]), ",".join(len_list)))

        if close_after_if_file_object:
            out_fd.close()


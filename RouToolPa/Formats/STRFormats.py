#!/usr/bin/env python

from collections import OrderedDict


STR_FMT_COLS = OrderedDict({
    "str": OrderedDict({"amplicon_len": 0,
                        "HitName": 1,
                        "FP_ID": 2,
                        "FP_seq": 3,
                        "FP_degeneracies": 4,
                        "FP_mismatches": 5,
                        "RP_ID": 6,
                        "RP_seq": 7,
                        "RP_degeneracies": 8,
                        "RP_mismatches": 9,
                        "RevcompRP": 10,
                        "ProbeID": 11,
                        "Probe_seq": 12,
                        "RevcompProbe": 13,
                        "Probe_degeneracies": 14,
                        "Probe_mismatches": 15,
                        "Probe_startOnAmplicon": 16,
                        "ProbeStrand": 17,
                        "Start": 18,
                        "End": 19,
                        "Full_Hit_ID": 20,
                        "SubjectFullLength": 21,
                        "AmpliconSeq": 22,
                        "Primary_tag": 23,
                        "Gene": 24,
                        "Product": 25,
                        "Protein_id": 26,
                        "Note": 27
                        }),

    "filtered_str": OrderedDict({"primer_pair": 0,
                                 "amplicon_len": 1,
                                 "HitName": 2,
                                 "FP_ID": 3,
                                 "FP_mismatches": 4,
                                 "RP_ID": 5,
                                 "RP_mismatches": 6,
                                 "Start": 7,
                                 "End": 8,
                                 "AmpliconSeq": 9,
                                 "max_mismatches": 10,
                                 "total_mismatches": 11,
                                 "max_mis_min_dist": 12,
                                 "tot_mis_min_dist": 13,
                                 "min_len": 14,
                                 "max_len": 15,
                                 "len_in_expected_interval": 16,
                                 "color": 17,
                                 })

    })

STR_FMT_COLUMN_NAMES_SYN = OrderedDict({
    "str": OrderedDict({
                        "primer_pair": None,
                        "scaffold_id": "HitName",
                        "start": "Start",
                        "end": "End",
                        "forward_primer_ID": "FP_ID",
                        "forward_primer_seq": "FP_seq",
                        "reverse_primer_ID": "RP_ID",
                        "reverse_primer_seq": "RP_seq",
                        "forward_primer_mismatches": "FP_mismatches",
                        "reverse_primer_mismatches": "RP_mismatches",
                        "amplicon_seq": "AmpliconSeq",
                        "amplicon_len": "amplicon_len",
                        "color": None,
                        }),

    "filtered_str": OrderedDict({
                                 "primer_pair": "primer_pair",
                                 "scaffold_id": "HitName",
                                 "start": "Start",
                                 "end": "End",
                                 "forward_primer_ID": "FP_ID",
                                 "reverse_primer_ID": "RP_ID",
                                 "forward_primer_mismatches": "FP_mismatches",
                                 "reverse_primer_mismatches": "RP_mismatches",
                                 "amplicon_seq": "AmpliconSeq",
                                 "amplicon_len": "amplicon_len",
                                 "color": "color",
                                 })

    })

#!/usr/bin/env python

from collections import OrderedDict


ALN_FMT_COLS = OrderedDict({
    "tab": OrderedDict({
                             "score": 0,
                             "target_id":   1,
                             "target_start": 2,
                             "target_hit_len": 3,
                             "target_strand": 4,
                             "target_len": 5,
                             "query_id": 6,
                             "query_start": 7,
                             "query_hit_len": 8,
                             "query_strand": 9,
                             "query_len": 10,
                             "alignment": 11,
                             "EG2": 12,
                             "E": 13,
                              }),

    "psl": OrderedDict({
                        "matches": 0,
                        "misMatches": 1,
                        "repMatches": 2,
                        "nCount": 3,
                        "qNumInsert": 4,
                        "qBaseInsert": 5,
                        "tNumInsert": 6,
                        "tBaseInsert": 7,
                        "strand": 8,
                        "qName": 9,
                        "qSize": 10,
                        "qStart": 11,
                        "qEnd": 12,
                        "tName": 13,
                        "tSize": 14,
                        "tStart": 15,
                        "tEnd": 16,
                        "blockCount": 17,
                        "blockSizes": 18,
                        "qStarts": 19,
                        "tStarts": 20,
                         }),

    })

ALN_FMT_COLUMN_NAMES_SYN = OrderedDict({
    "tab": OrderedDict({
        "target_id": "target_id",
        "target_start": "target_start",
        "target_end": None,
        "target_hit_len": "target_hit_len",
        "target_strand": "target_strand",
        "target_len": "target_len",
        "query_id": "query_id",
        "query_start": "query_start",
        "query_end": None,
        "query_hit_len": "query_hit_len",
        "query_strand": "query_strand",
        "query_len": "query_len",
        }),

    "psl": OrderedDict({
        "target_id": "tName",
        "target_start": "tStart",
        "target_end": "tEnd",
        "target_hit_len": "tHitLen",  # is calculated during parsing
        "target_strand": None,
        "target_len": "tSize",
        "query_id": "qName",
        "query_start": "qStart",
        "query_end": "qEnd",
        "query_hit_len": "qHitLen",  # is calculated during parsing
        "query_strand": "strand",
        "query_len": "qSize",
        })

    })

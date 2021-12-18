#!/usr/bin/env python

from collections import OrderedDict


SCAF_FMT_COLS = OrderedDict({
    "agp": OrderedDict({
                        "scaffold": 0,
                        "start": 1,
                        "end": 2,
                        "part_number": 3,
                        "part_type": 4,
                        "part_id/gap_length": 5,
                        "part_start/gap_type": 6,
                        "part_end/linkage": 7,
                        "orientation/evidence": 8
                         })
    })

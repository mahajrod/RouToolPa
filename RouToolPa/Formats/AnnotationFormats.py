#!/usr/bin/env python

from collections import OrderedDict

GFF_COLS = OrderedDict({
                        "scaffold": 0,
                        "source":   1,
                        "featuretype": 2,
                        "start": 3,
                        "end": 4,
                        "score": 5,
                        "strand": 6,
                        "phase": 7,
                        "attributes": 8
                        })

BED_COLS = OrderedDict({
                        "scaffold": 0,
                        "start": 1,
                        "end": 2
                        })
BUSCO_TABLE_COLS = OrderedDict({
                                "id": 0,
                                "status": 1,
                                "contig": 2,
                                "start": 3,
                                "end": 4,
                                "score": 5,
                                "length": 6})

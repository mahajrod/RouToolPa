#!/usr/bin/env python

from collections import OrderedDict

VCF_COLS = OrderedDict({"CHROM":  0,
                        "POS":    1,
                        "ID":     2,
                        "REF":    3,
                        "ALT":    4,
                        "QUAL":   5,
                        "FILTER": 6,
                        "INFO":   7,
                        "FORMAT": 8})

#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from RouToolPa.Routines import FileRoutines
from RouToolPa.Data import Tables

try:
    import importlib.resources as resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as resources


class VEPTable:
    def __init__(self, vep_tab_file, metadata=None):

        self.impact_type_list = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
        self.so_description = pd.read_csv(resources.open_text(Tables, "SO_terms.description.tsv"), sep="\t", header=True)
        self.so_description.columns = self.so_description.columns.str.lstrip("#")
        self.metadata = metadata if metadata else []
        self.records = self.read(vep_tab_file)

    def read(self, vep_tab_file):

        with FileRoutines.metaopen(vep_tab_file, "r") as in_fd:
            metadata = []
            while True:
                line = in_fd.readline()
                if line[:19] == "#Uploaded_variation":
                    header_list = line.strip()[1:].split("\t")
                    break
                metadata.append(line)

            return pd.read_csv(in_fd, sep='\t', names=header_list)




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
    def __init__(self, vep_tab_file, metadata=None, get_stats=False):
        self.so_description = pd.read_csv(resources.open_text(Tables, "SO_terms.description.tsv"),
                                          sep="\t", index="#SO term")
        self.so_description.index.name = "#SO term" # = columns = self.so_description.columns.str.lstrip("#")

        self.allowed_impact_type_list = self.so_description["IMPACT"].unique()
        self.allowed_so_terms = self.so_description["SO term"]

        self.metadata = metadata if metadata else []
        self.records = self.read(vep_tab_file)

        # stats
        self.impact_counts = None
        self.consequence_entry_set = None
        self.count_consequence = None

        if get_stats:
            self.get_stats()

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

    # ----------------------------------Stats---------------------------
    def get_stats(self):
        self.impact_counts = self.records[["Uploaded_variation", "IMPACT"]].groupby("IMPACT").count()
        self.consequence_entry_set = self.get_consequence_entry_set()
        self.count_consequence = self.count_consequence().join(self.so_description[["IMPACT"]])

    def get_consequence_entry_set(self):
        consequence_list = []
        for entry in list(map(lambda s: s.split(","), self.records["Consequence"].unique())):
            consequence_list += entry

        return set(consequence_list)

    def count_consequence(self):
        count_dict = dict([(element, 0) for element in self.consequence_entry_set])
        for entry in self.records["Consequence"]:
            for cons in entry.split(","):  # row can have multiple consequence terms
                count_dict[cons] += 1
        return pd.DataFrame.from_dict(count_dict, orient='index', columns=["Counts"]).sort_values(by="Counts",
                                                                                                  ascending=False)


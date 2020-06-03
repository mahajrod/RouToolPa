#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from RouToolPa.Routines import FileRoutines


class VEPTable:
    def __init__(self, vep_tab_file, metadata=None):

        self.metadata = metadata if metadata else []
        self.records = self.read(vep_tab_file)

    def read(self, vep_tab_file):

        with FileRoutines.metaopen(vep_tab_file, "r") as in_fd:
            metadata = []
            for line in in_fd:
                if line[:19] == "#Uploaded_variation":
                    header_list = line.strip()[1:].split("\t")
                    break
                metadata.append(line)

            return pd.read_csv(in_fd, sep='\t', names=header_list)




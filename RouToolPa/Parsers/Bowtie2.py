#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class Bowtie2Table:
    def __init__(self, files, samples=None):

        filelist = [files] if isinstance(files, str) else files

        self.sample_number = len(filelist)

        if samples:
            self.samplelist = [samples] if isinstance(samples, str) else samples
        else:
            self.samplelist = ["S%i" % i for i in range(1, self.sample_number + 1)]

        self.records = self.read(files, samples)

    def read(self, files, samples):
        records_df = []

        index = ["Input read pairs",
                 "Not aligned concordantly pairs",
                 "Uniquely and concordantly mapped pairs",
                 "Multi and concordantly mapped pairs",
                 "Uniquely and discordantly mapped pairs",
                 "Not aligned pairs",
                 ]

        for filename in files:
            with open(filename, "r") as in_fd:
                linelist = list(map(lambda s: s.strip().split()[0], list(in_fd)))
                linelist = list(map(int, linelist[1:5] + [linelist[7]] + [linelist[9]]))
                records_df.append(linelist)

        records_df = pd.DataFrame.from_records(records_df, index=samples, columns=index).transpose()

        return records_df

    def write(self, output_file):
        self.records.to_csv(output_file, sep="\t", index=True, header=True)

    def write_xlsx(self, output_file, sheet_name="bowtie2_stat"):
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')

        self.records.to_excel(writer, sheet_name=sheet_name, header=True, index=True, freeze_panes=(1, 1))
        #worksheet = writer.sheets[sheet_name]

        writer.save()

    def draw(self, output_prefix, width=0.35, figsize=(4, 4), dpi=300, extensions=("png",)):

        fig = plt.figure(1, figsize=figsize, dpi=dpi)

        ind = np.arange(self.sample_number)  # the x locations for the groups
        bottom = np.zeros(self.sample_number)

        bar_list = []
        for stat in ("Uniquely and concordantly mapped pairs",
                     "Multi and concordantly mapped pairs",
                     "Uniquely and discordantly mapped pairs",
                     "Not aligned pairs"):

            data = self.records.loc[stat] / self.records.loc["Input read pairs"]
            bar_list.append(plt.bar(ind, data, width, bottom=bottom))

            bottom += data

        plt.ylabel('Read pairs, %')
        plt.title('Alignment statistics')
        plt.xticks(ind, self.samplelist)
        plt.yticks(np.arange(0, 100, 10))

        plt.legend(map(lambda s: s[0], bar_list), ('Uniq and concordsnt', 'Multi and concordant', "Uniq and discordant", "Not aligned"))

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))




#!/usr/bin/env python
import pandas
__author__ = 'mahajrod'
from collections import OrderedDict
from RouToolPa.Collections.General import TwoLvlDict


class PLINKReport():
    def __init__(self, plink_report_file, report_type="ROH", separator="\s+"):
        #OrderedDict.__init__(self)
        self.report_type = report_type
        self.families = None
        self.samples = None
        self.roh_lengths = None

        index_col = None
        if report_type == "ROH":
            index_col = (0, 1)

        self.data = pandas.read_table(plink_report_file, sep=separator, index_col=index_col)

        if report_type == "ROH":
            self.families = self.data.index.levels[0]
            self.samples = self.data.index.levels[1]
        """
        with open(plink_report_file, "r") as in_fd:
            if report_type == "ROH":


                header_list = in_fd.readline().strip().split()

            for line in in_fd:
                tmp = line.strip().split("\t")
                if tmp[0] not in self:
                    self[tmp[0]] = OrderedDict()
                if tmp[1] not in self[tmp[0]]:
                    self[tmp[0]][tmp[1]] = OrderedDict()
                for i in range(2, len(header_list)):
                    self[tmp[0]][tmp[1]][header_list[i]] = tmp[i]
                for j in 2, 8, 10, 11, 12:
                    self[tmp[0]][tmp[1]][header_list[j]] = float(self[tmp[0]][tmp[1]][header_list[j]])
                for j in 6, 7, 9:
                    self[tmp[0]][tmp[1]][header_list[j]] = int(self[tmp[0]][tmp[1]][header_list[j]])
            """
    def __len__(self):
        return len(self.data)

    def get_roh_lengths(self, save=False):
        """
        return: pandas Series
        """
        if save:
            self.roh_lengths = self.data['POS2'] - self.data['POS1'] + 1
        return self.roh_lengths if self.roh_lengths is not None else self.data['POS2'] - self.data['POS1'] + 1


class PLINKReportCollection(OrderedDict):
    def __init__(self, report_dict=None, file_dict=None, report_type=None):
        OrderedDict.__init__(self)
        self.report_type = report_type
        if report_dict:
            for report_id in report_dict:
                self[report_id] = report_dict
        if file_dict:
            for file_id in file_dict:
                if file_id in self:
                    raise KeyError("Report with same id(%s) is already present in PLINKReportCollection" % str(file_id))
                self[file_id] = PLINKReport(file_dict[file_id])

    def get_general_stats(self):
        pass

    def write_general_stats(self, output_file):
        stat_dict = self.get_general_stats()
        stat_dict.write(output_file, absent_symbol=".", sort=False)

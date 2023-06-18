#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from collections import OrderedDict
from collections.abc import Iterable
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from RouToolPa.Parsers.Abstract import Record, Collection, Metadata, Header
from RouToolPa.Parsers.VCF import CollectionVCF, MetadataVCF, HeaderVCF

from RouToolPa.General.GeneralCollections import WDict


class RecordCCF(Record, Iterable):

    def __init__(self, size=None, start=None, end=None, info_dict=None, flags=None,
                 collection_vcf=None, bad_vcf_records=0, from_records=True, subclusters=None, features=None,
                 ):
        # possible flags:
        # IP - indel(s) are present in record
        # BR - record is located in bad region
        self.records = collection_vcf
        self.subclusters = subclusters
        if from_records:
            self.size = len(collection_vcf)
            self.start = collection_vcf[0].pos
            self.end = collection_vcf[-1].pos - 1 + \
                       max(map(lambda x: len(x), collection_vcf[-1].alt_list + [collection_vcf[-1].ref]))
            self.flags = set(flags) if flags else set()
            for record in self.records:
                if record.check_indel():
                    self.flags.add("IP")
                    break
            self.info_dict = info_dict if info_dict else OrderedDict({})
            self.features = features

            self.bad_records = 0
            for record in self.records:
                if "BR" in record.flags:
                    self.bad_records += 1
        else:
            self.size = size
            self.start = start
            self.end = end
            self.features = features
            self.info_dict = info_dict
            self.mean_dist = None
            self.flags = set(flags) if flags else set()
            self.bad_records = bad_vcf_records
        self.len = self.end - self.start + 1
        if bad_vcf_records > 0:
            self.flags.add("BR")
        else:
            self.flags - set(["BR"])

        distances = self.distances()
        if self.info_dict is None:
            self.info_dict = OrderedDict({})

        if self.size > 1:
            self.info_dict["Mean"] = np.mean(distances)
            self.info_dict["Median"] = np.median(distances)
            self.info_dict["Power"] = self.size / np.median(distances)

    def __len__(self):
        return self.size

    def __iter__(self):
        for record in self.records:
            yield record

    def __str__(self):
        attributes_string = "Size=%i;Bad_records=%i" % (self.size, self.bad_records)
        if self.flags:
            attributes_string += ";" + ";".join(self.flags)
        if self.info_dict:
            attributes_string += ";" + ";".join(["%s=%s" % (key, ",".join(map(lambda x: str(x), self.info_dict[key])) if isinstance(self.info_dict[key], list) or isinstance(self.info_dict[key], set)  else str(self.info_dict[key]))
                                                 for key in self.info_dict])
        if self.subclusters is not None:
            #print(self.subclusters)
            attributes_string += ";Subclusters=" + ",".join(map(lambda x: str(x), self.subclusters))

        cluster_string = "%i\t%i\t%s" % (self.start, self.end, attributes_string)
        return cluster_string + "\n\t" + "\n\t".join([str(record) for record in self.records])

    def distances(self):
        #returns distances between variants in cluster
        positions = np.array([record.pos for record in self])
        return np.ediff1d(positions)

    def reset_flags(self, flags_to_reset='all'):
        # flags_to_reset can take values 'all', set or list, or None
        # if all - all flags will be removed
        # if None - nothing will be done
        # if set/list - only flags in set/list will be removed
        if flags_to_reset == 'all':
            self.flags = set([])
        elif flags_to_reset:
            flags_to_remove = set(flags_to_reset)
            self.flags = self.flags - flags_to_remove

    def check_flags(self, flag_list, mismatch_list=[], expression_list=[], remove_mismatch_list=[],
                    flags_to_reset=None, mode="all", min_cluster_size=[]):
        # TODO: fix changes in  loc in info_dict after removal
        # at moment subclustering resets after check if some record was removed

        # allow_flag_mismatch - number of allowed mismaches. if zero or None no mismatches are allowed
        self.reset_flags(flags_to_reset)
        if self.flags is None:
            self.flags = set([])


        mismatches = mismatch_list if mismatch_list else [0 for flag in flag_list]
        min_cluster = min_cluster_size if min_cluster_size else [3 for flag in flag_list]
        remove_list = remove_mismatch_list if remove_mismatch_list else[False for flag in flag_list]
        # possible from_records_flag_mode:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list

        # option flags_to_consider works only for remove_mismatch_records
        if mode == "one":
            for record in self:
                self.flags |= record.flags
        elif mode == "all":
            record_to_remove_dict = dict([(flag, set([])) for flag in flag_list])
            mismatch_count_dict = dict([(flag, 0) for flag in flag_list])
            index = 0
            for record in self:
                expressions = [eval(expression) for expression in expression_list] \
                    if expression_list else [True for flag in flag_list]
                for flag, expression, min_size in zip(flag_list, expressions, min_cluster):
                    #if self.size < min_size:
                    #    continue
                    if flag not in record.flags:
                        if mismatch_count_dict[flag] is not None:
                            mismatch_count_dict[flag] += 1
                        if expression and self.size >= min_size:# expression to reconsider cluster with mismatch as good
                            record_to_remove_dict[flag].add(index)
                        else:
                            # if expression is not True flag will not be set, None is used to indicate it
                            mismatch_count_dict[flag] = None
                index += 1
            for flag, mismatch, min_size in zip(flag_list, mismatches, min_cluster):
                if (mismatch_count_dict[flag] is not None) and mismatch_count_dict[flag] <= mismatch:

                    self.flags.add(flag)
                    if mismatch_count_dict[flag] > 0:
                        self.info_dict["N%sR" % flag] = [mismatch_count_dict[flag]]
                else:
                    record_to_remove_dict[flag].clear()
            records_to_remove = set([])
            for flag, remove in zip(flag_list, remove_list):
                if remove and (mismatch_count_dict[flag] is not None):
                    records_to_remove |= record_to_remove_dict[flag]
            records_to_remove = sorted(records_to_remove, reverse=True)
            #print(records_to_remove)
            record_scaffold = self.records.scaffold_list[0]
            if records_to_remove:
                #print(self)
                for index in records_to_remove:
                    self.records.records[record_scaffold].pop(index)
                    self.subclusters = np.delete(self.subclusters, index)
                self.__init__(collection_vcf=CollectionVCF(metadata=self.records.metadata,
                                                           records_dict=self.records.records,
                                                           header=self.records.header,
                                                           samples=self.records.samples, from_file=False),
                              info_dict=self.info_dict,
                              flags=self.flags, subclusters=self.subclusters,
                              from_records=True)

    def check_location(self):
        self.bad_records = 0
        for variant in self:
            if "BR" in variant.flags:
                self.bad_records += 1
        if self.bad_records > 0:
            self.flags.add("BR")

    def get_location(self, record_scaffold, record_dict, key="Loc", use_synonym=False, synonym_dict=None,
                     feature_type_black_list=[], strand_key=None):
        # function is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        # id of one SeqRecord in record_dict must be equal to record.pos
        if not self.info_dict:
            self.info_dict = {}
        if not self.features:
            self.features = []
        if key not in self.info_dict:
            self.info_dict[key] = set([])
        for variant in self:
            if key in variant.info_dict:
                self.info_dict[key] |= set(variant.info_dict[key])
            for feature in record_dict[record_scaffold].features:
                if (variant.pos - 1) in feature:
                    self.features.append(feature)
                    self.info_dict[key].add(self.get_synonym(feature.type, use_synonym=use_synonym,
                                                               synonym_dict=synonym_dict))
                for sub_feature in feature.sub_features:
                    if (variant.pos - 1) in sub_feature:
                        self.info_dict[key].add(self.get_synonym(sub_feature.type, use_synonym=use_synonym,
                                                                   synonym_dict=synonym_dict))

    def subclustering(self,
                      method="inconsistent",
                      threshold=0.8,
                      cluster_distance='average'):
        tmp = self.records.get_clusters(extracting_method=method,
                                        threshold=threshold,
                                        cluster_distance=cluster_distance,
                                        draw_dendrogramm=False,
                                        return_collection=False,
                                        write_inconsistent=False,
                                        write_correlation=False)
        self.subclusters = tmp[tmp.keys()[0]]

    def adjust(self, border_limit=None, min_size_to_adjust=2, remove_border_subclusters=False, remove_size_limit=1):
        # adjusts cluster borders, returns list of new cluster records
        # skip adjustment for clusters with 3 or less mutations
        if (self.size < min_size_to_adjust) or (self.subclusters is None):
            #return -1
            return [self]
        limit = border_limit if border_limit else len(self.subclusters)
        for i in range(0, limit):
            if self.subclusters[i] == self.subclusters[0]:
                left_subcluster_end = i
            else:
                break
        # exit if cluster doesnt have subclusters
        if left_subcluster_end == len(self.subclusters) - 1:
            #return 1
            return [self]

        for i in range(-1, -limit - 1, -1):
            if self.subclusters[i] == self.subclusters[-1]:
                right_subcluster_start = i
            else:
                break

        if remove_border_subclusters:
            start = left_subcluster_end + 1 if left_subcluster_end < remove_size_limit else 0
            end = right_subcluster_start if right_subcluster_start >= -remove_size_limit else len(self.subclusters)

            new_left_cluster, new_right_cluster = None, None
            record_scaffold = self.records.scaffold_list[0]

            if start > 0:
                new_left_cluster = RecordCCF(collection_vcf=CollectionVCF(records_dict=dict([(record_scaffold, self.records.records[record_scaffold][:start])]),
                                                                          from_file=False),
                                             subclusters=self.subclusters[:start], from_records=True)

            if end < len(self.subclusters):
                new_right_cluster = RecordCCF(collection_vcf=CollectionVCF(records_dict=dict([(record_scaffold, self.records.records[record_scaffold][end:])]),
                                                                           from_file=False),
                                              subclusters=self.subclusters[end:], from_records=True)

            new_middle_cluster = RecordCCF(collection_vcf=CollectionVCF(records_dict=dict([(record_scaffold, self.records.records[record_scaffold][start:end])]),
                                                                        from_file=False),
                                           subclusters=self.subclusters[start:end], from_records=True)

            cluster_list = [new_left_cluster] if new_left_cluster else []
            cluster_list += [new_middle_cluster]
            cluster_list += [new_right_cluster] if new_right_cluster else []
            return cluster_list

    def check_strandness(self):
        #for desaminases only
        count_C = 0.
        for record in self:
            if record.ref == "C":
                count_C += 1.
        homogeneity = count_C / self.size
        self.info_dict["Homogeneity"] = homogeneity if homogeneity >= 0.5 else 1.0 - homogeneity
        self.info_dict["Strand"] = "P" if homogeneity >= 0.5 else "N"

    def gff_string(self):
        attributes_string = "Size=%i;Bad_records=%i" % (self.size, self.bad_records)
        if self.info_dict:
            for key in self.info_dict:
                attributes_string += ";%s=%s" % (key, self.info_dict[key])
        if self.flags:
            attributes_string += ";" + ";".join(self.flags)
        if self.subclusters != None:
            #print(self.subclusters)
            attributes_string += ";Subclusters=" + ",".join(map(lambda x: str(x), self.subclusters))
        return "%s\t%s\t%i\t%i\t.\t.\t.\t%s" % ("custom", "cluster", self.start, self.end, attributes_string)


class MetadataCCF(Metadata):

    def __init__(self, samples, vcf_metadata=None, vcf_header=None, metadata={}):
        self.samples = samples      #list
        self.metadata = metadata
        self.vcf_metadata = vcf_metadata
        self.vcf_header = vcf_header

    def __str__(self):
        metadata_string = None
        if self.vcf_metadata:
            metadata_string = "#VCF_METADATA START\n" + str(self.vcf_metadata) + "\n#VCF_METADATA END"
        if self.vcf_header:
            metadata_string = metadata_string + "\n#VCF_HEADER\n" + str(self.vcf_header) \
                if metadata_string else "#VCF_HEADER\n" + str(self.vcf_header)
        if self.metadata:
            metadata_string += "\n##" + "\n##".join(["%s=%s" % (key, self.metadata[key]) for key in self.metadata])
        return metadata_string


class HeaderCCF(list, Header):
    def __str__(self):
        return "#CCF_HEADER\n#" + "\t".join(self)


class CollectionCCF(Collection):
    """
    def read(self, input_file):
        # TODO: write read from ccf file
        with open(input_file, "r") as in_fd:

                stripped_line = in_fd.readline().strip()
                if stripped_line == "#VCF_METADATA START":
                    vcf_metadata = MetadataVCF()
                    stripped_line = in_fd.readline().strip()
                    while (stripped_line != "#VCF_METADATA END"):
                        vcf_metadata.add_metadata(stripped_line)
                        stripped_line = in_fd.readline().strip()
                stripped_line = in_fd.readline().strip()
                if stripped_line == "#VCF_HEADER":
                    header_line = in_fd.readline().strip()
                    vcf_header = HeaderVCF(header_line[1:].split("\t"))
                    #print("a\na\na\na\na\n")
                    #print(vcf_header)
                    self.metadata = MetadataCCF(vcf_header[9:], vcf_metadata=vcf_metadata, vcf_header=vcf_header)
                stripped_line = in_fd.readline().strip()
                if stripped_line == "#CCF_HEADER":
                    header_line = in_fd.readline().strip()
                    self.header = HeaderCCF(header_line[1:].split("\t"))
                flag = 0
                self.records = []
                while True:
                    data_line = in_fd.readline()

                    if data_line == "" or data_line == "\n":
                        break
                    stripped_line = data_line.strip()
                    if data_line[0] == "\t":
                        #stripped_line = stripped_line[1:]
                        #print(collection_vcf)
                        collection_vcf.records.append(collection_vcf.add_record(stripped_line, external_metadata=self.metadata.vcf_metadata))
                        flag = 1
                        #print("aaaa")
                        continue

                    if flag != 0:
                        self.records.append(RecordCCF(id=cluster_id, chrom=chrom, size=size, start=start, end=end,
                                                      info_dict=info_dict, flags=flags,
                                                      collection_vcf=collection_vcf, bad_vcf_records=bad_records,
                                                      from_records=False, subclusters=subclusters))
                        #collection_vcf = None

                    if stripped_line[0] == ">":
                        flag = 0
                        cluster_id, chrom, start, end, info_dict_and_flags = stripped_line[1:].split("\t")
                        start = int(start)
                        end = int(end)
                        info_dict_and_flags = info_dict_and_flags.split(";")
                        info_dict = OrderedDict({})
                        flags = set([])
                        subclusters = None
                        for descr_entry in info_dict_and_flags:
                            descr_entry_splited = descr_entry.split("=")
                            if len(descr_entry_splited) == 1:
                                flags.add(descr_entry_splited[0])
                                continue
                            if descr_entry_splited[0] == "Size":
                                size = int(descr_entry_splited[1])
                            elif descr_entry_splited[0] == "Bad_records":
                                bad_records = int(descr_entry_splited[1])
                            elif descr_entry_splited[0] == "Mean" or descr_entry_splited[0] == "Median" or descr_entry_splited[0] == "Power" or descr_entry_splited[0] == "Homogeneity":
                                info_dict[descr_entry_splited[0]] = float(descr_entry_splited[1])
                            elif descr_entry_splited[0] == "Loc":
                                info_dict[descr_entry_splited[0]] = descr_entry_splited[1].split(",")
                            elif descr_entry_splited[0] == "Subclusters":
                                subclusters = [int(x) for x in descr_entry_splited[1].split(",")]
                            else:
                                info_dict[descr_entry_splited[0]] = descr_entry_splited[1].split(",")
                                if len(info_dict[descr_entry_splited[0]]) == 1:
                                    info_dict[descr_entry_splited[0]] = info_dict[descr_entry_splited[0]][0]
                        collection_vcf = CollectionVCF(metadata=None, record_list=None, header=None, vcf_file=None, samples=None, from_file=False, external_metadata=None)
                        continue
                self.records.append(RecordCCF(id=cluster_id, chrom=chrom, size=size, start=start, end=end,
                                              info_dict=info_dict, flags=flags,
                                              collection_vcf=collection_vcf, bad_vcf_records=bad_records,
                                              from_records=False, subclusters=subclusters))
    """
    def filter(self, expression):
        # expression should be a function with one argument - record entry
        filtered_records, filtered_out_records = self.filter_records(expression)
        return CollectionCCF(metadata=self.metadata, records_dict=filtered_records,
                             from_file=False,
                             header=self.header), \
               CollectionCCF(metadata=self.metadata, records_dict=filtered_out_records,
                             from_file=False,
                             header=self.header)

    def filter_by_size(self, min_size=3):
        return self.filter(lambda record: record.size >= min_size)

    def count(self):
        sizes = []
        for record in self:
            sizes.append(record.size)
        return sizes

    def statistics(self, filename="cluster_size_distribution.svg", title="Distribution of sizes of clusters",
                   dpi=150, figsize=(8, 8), facecolor="green", draw_hist=True):
        plt.figure(1, dpi=dpi, figsize=figsize)
        plt.subplot(1, 1, 1)
        plt.suptitle(title)
        counts = self.count()
        if not counts:
            return -1
        maximum = max(counts)
        #print(maximum)
        bins = np.linspace(1, maximum+1, maximum+1)
        #print(bins)
        n, bins, patches = plt.hist(counts, bins, facecolor=facecolor, align='left')
        if draw_hist:
            plt.xticks(np.arange(1, maximum+1, 1))
            plt.xlim(xmin=0.5, xmax=maximum+0.5)
            plt.xlabel("Size")
            plt.ylabel("Number")
            plt.savefig(filename)
        plt.close()
        return WDict([map(int, pair) for pair in zip(bins, n)])

    def filter_by_flags(self, white_flag_list=[], black_flag_list=[]):
        white_list = set(white_flag_list)
        black_list = set(black_flag_list)

        if white_list and not black_list:
            expression = lambda record: white_list & record.flags
        elif black_list and not white_list:
            expression = lambda record: not (black_list & record.flags)
        elif black_list and not white_list:
            expression = lambda record: (white_list & record.flags) and not (black_list & record.flags)
        elif not black_list and (not white_list):
            raise ValueError("Both white and black lists were not set")
        return self.filter(expression)

    def subclustering(self,
                      method="inconsistent",
                      threshold=0.8,
                      cluster_distance='average'):
        for record in self:
            if len(record) < 3:
                continue
            record.subclustering(method=method,
                                 threshold=threshold,
                                 cluster_distance=cluster_distance)

    def check_location(self):
        for record in self:
            record.check_location()

    def check_flags(self, flag_list, mismatch_list=[], expression_list=[], remove_mismatch_list=[],
                    flags_to_reset=None, mode="all", min_cluster_size=[]):
        for record in self:
            record.check_flags(flag_list, mismatch_list=mismatch_list, expression_list=expression_list,
                               remove_mismatch_list=remove_mismatch_list, flags_to_reset=flags_to_reset, mode=mode,
                               min_cluster_size=min_cluster_size)

    def check_strandness(self):
        for record in self:
            record.check_strandness()

    def get_data_for_stat(self, additional_data=("Median", "Mean", "Power")):
        data = []
        for record in self:
            if record.size == 1:
                continue
            data.append([record.len, record.size] + ([record.info_dict[add_data] for add_data in additional_data] if additional_data else []))
        return np.array(data)

    def heatmap_statistics(self, filename="heatmap_statistics.svg", suptitle="Heatmap_statistics",
                           dpi=150, facecolor="green", n_bins_default=20,
                           additional_data=("Median", "Mean", "Power")):
        labels_dict = {"Median": "Median distance",
                       "Mean": "Mean distance",
                       "Power": "Power(Size/median dist.)",
                       "Homogeneity": "Strand specificity"}
        names = ("Length", "Size")
        if additional_data is not None:
            names += additional_data
        data = [self.get_data_for_stat(additional_data=additional_data)[:, i]
                for i in range(0, 2 if additional_data is None else 2 + len(additional_data))]
        data.append(data[1] / data[2])
        data.append((data[1]**2) / data[2])
        size = 6 * len(names)
        plt.figure(1, dpi=dpi, figsize=(size, size))
        plt.suptitle(suptitle)

        if len(data[0]) == 0:
            return 0

        for i in range(0, len(names)):
            for j in range(0, len(names)):
                if i == j:
                    continue
                plt.subplot(len(names), len(names), i * len(names) + j + 1)
                plt.xlabel(labels_dict[names[i]] if names[i] in labels_dict else names[i])
                plt.ylabel(labels_dict[names[j]] if names[j] in labels_dict else names[j])
                n_x_bins = 10 if names[i] == "Homogeneity" else int(max(data[i])) if names[i] == "Size" else \
                    int(max(data[i]) * 20) if names[i] == "Power" else n_bins_default
                n_y_bins = 10 if names[j] == "Homogeneity" else int(max(data[j])) if names[j] == "Size" else \
                    int(max(data[j]) * 20) if names[j] == "Power" else n_bins_default
                #print(names[i])
                #print(data[i])
                #print(names[j])
                #print(data[j])
                #print(n_x_bins, n_y_bins)

                #cmap = colors.ListedColormap(['white', 'red'])
                #bounds = [0, 5, 10]
                #norm = colors.BoundaryNorm(bounds, cmap.N)

                plt.hist2d(data[i], data[j], (n_x_bins, n_y_bins), cmin=1)  # normed=True)
                plt.colorbar()
        plt.savefig(filename)
        plt.close()

    def get_data_for_strand_stat(self):
        num_of_P_clusters = 0
        num_of_N_clusters = 0
        strandness_P_clusters = []
        strandness_N_clusters = []
        size_P_clusters = []
        size_N_clusters = []
        length_P_clusters = []
        length_N_clusters = []

        for record in self:
            if "Strand" not in record.info_dict:
                print(record)
            if record.info_dict["Strand"] == "P":
                num_of_P_clusters += 1
                strandness_P_clusters.append(record.info_dict["SHom"])
                size_P_clusters.append(record.size)
                length_P_clusters.append(record.len)
            else:
                num_of_N_clusters += 1
                strandness_N_clusters.append(record.info_dict["SHom"])
                size_N_clusters.append(record.size)
                length_N_clusters.append(record.len)

        return num_of_P_clusters, strandness_P_clusters, size_P_clusters, length_P_clusters, \
               num_of_N_clusters,  strandness_N_clusters, size_N_clusters, length_N_clusters,

    def strandness_statistics(self, filename="Strandness_statistics.svg", suptitle="Cluster strandness_statistics",
                   dpi=150, figsize=(20, 20), facecolor="green"):
        num_of_P_clusters, strandness_P_clusters, size_P_clusters, length_P_clusters, \
        num_of_N_clusters,  strandness_N_clusters, size_N_clusters, length_N_clusters = \
                                                                                     self.get_data_for_strand_stat()

        plt.figure(1, dpi=dpi, figsize=figsize)
        plt.subplot(3, 4, 1)
        points = np.arange(2)
        bar_width = 0.55

        rects1 = plt.bar(points, [num_of_P_clusters, num_of_N_clusters], bar_width,
                         color='b')

        plt.xlabel('Strandness')
        plt.ylabel('Counts')
        plt.title("Distribution of clusters in strands")
        plt.xticks(points + bar_width, ('P', 'M'))

        for subplot_index, counts, title in zip([5, 9], [strandness_P_clusters, strandness_N_clusters], ["P", "N"]):
            plt.subplot(3, 4, subplot_index)
            if len(counts) == 0:
                continue
            bins = np.linspace(0.5, 1.0, 11)
            plt.hist(counts, bins, facecolor=facecolor)
            plt.xticks(np.arange(0.5, 1.0, 0.1))
            plt.title("Strandness coefficient of clusters in %s strand" % title)

        for subplot_index, coeff, size, title in zip([6, 10], [strandness_P_clusters, strandness_N_clusters],
                                                     [size_P_clusters, size_N_clusters], ["P", "N"]):
            plt.subplot(3, 4, subplot_index)
            if len(coeff) == 0:
                continue
            #plt.plot(size, coeff, "b.")
            #heatmap, xedges, yedges = np.histogram2d(size, coeff, bins=(10, max(size)))
            #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            #plt.imshow(heatmap)#, extent=extent)
            plt.hist2d(size, coeff, (max(size), 10))
            #plt.yticks(np.arange(0.5, 1.0, 0.1))
            plt.title("Coefficient and size of clusters in %s strand" % title)
            plt.xlabel("Size of cluster")
            plt.ylabel("Strandness")

        for subplot_index, coeff, length, title in zip([7, 11], [strandness_P_clusters, strandness_N_clusters],
                                                     [length_P_clusters, length_N_clusters], ["P", "N"]):
            plt.subplot(3, 4, subplot_index)
            if len(coeff) == 0:
                continue
            plt.plot(length, coeff, "b.")
            plt.title("Coefficient and length of clusters in %s strand" % title)
            plt.xlabel("Length of cluster")
            plt.ylabel("Strandness")

        for subplot_index, size, length, title in zip([8, 12], [size_P_clusters, size_N_clusters],
                                                     [length_P_clusters, length_N_clusters], ["P", "N"]):
            plt.subplot(3, 4, subplot_index)
            if len(size) == 0:
                continue
            plt.plot(length, size, "b.")
            plt.title("Length and size of clusters in %s strand" % title)
            plt.xlabel("Length of cluster")
            plt.ylabel("Size of cluster")

        plt.suptitle(suptitle)
        plt.savefig(filename)
        plt.close()

    def adjust(self, border_limit=None, min_size_to_adjust=2, remove_border_subclusters=False, remove_size_limit=1):
        new_records_dict = OrderedDict()
        for scaffold in self.scaffold_list:
            new_records_dict[scaffold] = []
        for scaffold in self.scaffold_list:
            for record in self.records[scaffold]:
                new_records_dict[scaffold] += record.adjust(border_limit=border_limit,
                                                            min_size_to_adjust=min_size_to_adjust,
                                                            remove_border_subclusters=remove_border_subclusters,
                                                            remove_size_limit=remove_size_limit)
        self.records = new_records_dict

    def extract_vcf(self):
        vcf_records_dict = OrderedDict()
        for scaffold in self.records:
            vcf_records_dict[scaffold] = []
            for record_ccf in self.records[scaffold]:
                vcf_records_dict[scaffold] += record_ccf.records.records[scaffold]
        return CollectionVCF(metadata=self.metadata.vcf_metadata,
                             records_dict=vcf_records_dict,
                             header=self.metadata.vcf_header,
                             samples=self.metadata.samples, from_file=False)

    def write_gff(self, outfile):
        with open(outfile, "w") as out_fd:
            for scaffold in self.scaffold_list:
                for record in self.records[scaffold]:
                    out_fd.write("%s\t%s\n" % (scaffold, record.gff_string()))
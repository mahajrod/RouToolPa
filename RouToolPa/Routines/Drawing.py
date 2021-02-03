__author__ = 'mahajrod'
import sys
import datetime
from collections import OrderedDict
import numpy as np
#import matplotlib

import pandas as pd

#matplotlib.use('Agg')
#os.environ['MPLCONFIGDIR'] = '/tmp/'

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

if sys.version_info[0] == 2:
    sys.stderr.write("WARNING!!! Package venn doesn't work in Python 2, related functionality is disabled!"
                     "Run scripts with Python 3...")
else:
    import venn

import scipy.stats as stats

from BCBio import GFF

from Bio import AlignIO

from RouToolPa.Collections.General import SynDict
from RouToolPa.Routines.Matplotlib import MatplotlibRoutines
from RouToolPa.Routines.Sequence import SequenceRoutines
from RouToolPa.Parsers.DESeq2 import CollectionPWC
from RouToolPa.Parsers.LAST import CollectionLast
from RouToolPa.Parsers.PSL import CollectionPSL


class DrawingRoutines(MatplotlibRoutines, SequenceRoutines):
    def __init__(self):
        MatplotlibRoutines.__init__(self)

    @staticmethod
    def millions(x, pos):
        return '%1.1fMbp' % (x*1e-6)

    @staticmethod
    def billions(x, pos):
        return '%1.1fGbp' % (x*1e-9)

    def draw_chromosomes_with_features_simple(self, chromosomes_gff, genes_gff, output_prefix, figsize=(10, 10),
                                              sense_feature_color="green", antisense_feature_color="red",
                                              chromosome_color="black", label_fontsize=15,
                                              chromosome_name_field=None,
                                              ext_list=("png",), dpi=None,
                                              alias_fields_in_gff=("ID", "Name", "gene", "Alias"),
                                              id_field_in_gff="ID", deseq2_pwc_file=None, upregulated_color="green",
                                              downregulated_color="red", absent_expression_data_color="black",
                                              coloring_mode="strand"):
        if deseq2_pwc_file and (coloring_mode == "expression"):
            deseq2_pwc_collection = CollectionPWC(from_file=True, pwc_file=deseq2_pwc_file)

        figure = plt.figure(figsize=figsize, dpi=dpi)
        """
        if dpi:
            figure = plt.figure(figsize=figsize, dpi=dpi)
        else:
            figure = plt.figure(figsize=figsize)
        """
        subplot = plt.subplot(1, 1, 1)

        subplot.get_yaxis().set_visible(False)
        #subplot.get_xaxis().set_visible(False)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')
        subplot.spines['right'].set_color('none')
        subplot.spines['left'].set_color('none')
        subplot.spines['top'].set_color('none')

        chr_dict = OrderedDict()

        max_chr_len = 0
        chr_number = 0
        with open(chromosomes_gff, "r") as chr_fd:
            for line in chr_fd:
                if line[0] == "#":
                    continue
                line_list = line.split("\t")
                feature_scaffold = line_list[0]
                feature_type = line_list[2]
                feature_start = line_list[3]
                feature_end = line_list[4]
                feature_description_list = line_list[8].split(";")

                chr_name = None
                if chromosome_name_field:

                    for entry in feature_description_list:
                        entry_list = entry.split("=")
                        if entry_list[0] == chromosome_name_field:
                            chr_name = entry_list[1]

                if feature_type == "chromosome" or feature_type == "plasmid" or feature_type == "region":
                    chr_dict[feature_scaffold] = OrderedDict()
                    chr_dict[feature_scaffold]["Name"] = chr_name if chr_name else feature_scaffold
                    chr_dict[feature_scaffold]["Genes"] = OrderedDict()
                    chr_dict[feature_scaffold]["Length"] = int(feature_end)
                    if chr_dict[feature_scaffold]["Length"] > max_chr_len:
                        max_chr_len = chr_dict[feature_scaffold]["Length"]
                        chr_number += 1
                elif feature_type == "centromere":
                    # add check for case when centromere line is first in gff file

                    chr_dict[feature_scaffold]["Centromere"] = [int(feature_start), int(feature_end)]
        #print chr_dict
        if genes_gff:
            with open(genes_gff, "r") as gene_fd:
                for line in gene_fd:
                    if line[0] == "#":
                        continue
                    line_list = line.strip().split("\t")
                    feature_scaffold = line_list[0]
                    feature_type = line_list[2]
                    feature_start = line_list[3]
                    feature_end = line_list[4]
                    feature_strand = line_list[6]
                    description_list = line_list[-1].split(";")

                    entry_id = None
                    alias_list = []
                    #print "AAAAAAAAAAA"
                    #print line
                    for entry in description_list:
                        if entry[:len(id_field_in_gff)+1] == ("%s=" % id_field_in_gff):
                            entry_id = entry[len(id_field_in_gff)+1:]
                            alias_list.append(entry_id)
                        for alias_id in alias_fields_in_gff:
                            if entry[:len(alias_id)+1] == ("%s=" % alias_id):
                                alias_list += entry.split("=")[1].split(",")
                    #print entry_id
                    if not entry_id:
                        continue

                    chr_dict[feature_scaffold]["Genes"][entry_id] = [int(feature_start), int(feature_end), feature_strand, alias_list]

        #print chr_dict
        centromere_radius = int(max_chr_len/100)
        distance_between_chromosome_and_gene = centromere_radius * 2
        distance_between_chromosomes = centromere_radius * 8
        chromosome_width = int(centromere_radius / 2)
        gene_width = chromosome_width
        chromosome_position = - int(distance_between_chromosomes / 3)

        text_x_offset = -max_chr_len/15
        for chromosome in chr_dict:
            print("Drawing chromosome %s" % chr_dict[chromosome]["Name"])
            chromosome_position += distance_between_chromosomes
            chromosome_fragments_list = []
            if "Centromere" not in chr_dict[chromosome]:
                chromosome_fragments_list.append(Rectangle((1, chromosome_position), chr_dict[chromosome]["Length"],
                                                           chromosome_width, fill=True, edgecolor=chromosome_color,
                                                           facecolor=chromosome_color))
            else:
                # left arm of chromosome
                #print "Centromere"
                #print "Left arm"
                #print (1, chromosome_position), chr_dict[chromosome]["Centromere"][0], chromosome_width


                chromosome_fragments_list.append(Rectangle((1, chromosome_position), chr_dict[chromosome]["Centromere"][0],
                                                           chromosome_width, fill=True, edgecolor=chromosome_color,
                                                           facecolor=chromosome_color))
                chromosome_fragments_list.append(Circle((chr_dict[chromosome]["Centromere"][0] + centromere_radius,
                                                         chromosome_position + chromosome_width/2), centromere_radius,
                                                        fill=False, edgecolor=chromosome_color,
                                                        facecolor=chromosome_color))
                chromosome_fragments_list.append(Rectangle((chr_dict[chromosome]["Centromere"][0] + 2 * centromere_radius,
                                                            chromosome_position), chr_dict[chromosome]["Length"] - chr_dict[chromosome]["Centromere"][1],
                                                           chromosome_width, fill=True, edgecolor=chromosome_color,
                                                           facecolor=chromosome_color))
            for patch in chromosome_fragments_list:
                subplot.add_patch(patch)

            subplot.annotate(chr_dict[chromosome]["Name"], xy=(text_x_offset, chromosome_position),
                             xycoords='data', fontsize=label_fontsize,
                             ha='center', va='center')

            if chr_dict[chromosome]["Genes"]:
                for gene in chr_dict[chromosome]["Genes"]:
                    print("Adding feature %s: %i-%i, %s. Aliases: %s" % (gene, chr_dict[chromosome]["Genes"][gene][0],
                                                                         chr_dict[chromosome]["Genes"][gene][1],
                                                                         chr_dict[chromosome]["Genes"][gene][2],
                                                                         ",".join(chr_dict[chromosome]["Genes"][gene][3])))
                    gene_start = chr_dict[chromosome]["Genes"][gene][0]
                    if "Centromere" in chr_dict[chromosome]:
                        if gene_start >= chr_dict[chromosome]["Centromere"][1]:
                            gene_start += centromere_radius * 2

                    if deseq2_pwc_file and (coloring_mode == "expression"):
                        for alias in chr_dict[chromosome]["Genes"][gene][3]:
                            if alias in deseq2_pwc_collection:
                                gene_log2foldchange = deseq2_pwc_collection[alias].log2foldchange
                                break
                        else:
                            gene_log2foldchange = None
                        gene_color = absent_expression_data_color if gene_log2foldchange is None else upregulated_color if gene_log2foldchange > 0 else downregulated_color
                    else:
                        gene_color = sense_feature_color if chr_dict[chromosome]["Genes"][gene][2] == "+" else antisense_feature_color

                    gene_patch = Rectangle((gene_start, chromosome_position + (1 if chr_dict[chromosome]["Genes"][gene][2] == "+" else -1) * distance_between_chromosome_and_gene),
                                           chr_dict[chromosome]["Genes"][gene][1] - chr_dict[chromosome]["Genes"][gene][0] + 1,
                                           gene_width, fill=True, edgecolor=gene_color,
                                           facecolor=gene_color)
                    subplot.add_patch(gene_patch)
        plt.xlim(xmax=max_chr_len, xmin=text_x_offset)
        plt.ylim(ymax=chromosome_position+distance_between_chromosomes)
        plt.subplots_adjust(right=0.95)#bottom=0.1, right=0.8, top=0.9)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension), dpi=dpi)



    def draw_alignment(self, alignment, features, output_prefix, record_style=None, ext_list=["svg", "png"],
                       label_fontsize=13, left_offset=0.2, figure_width=8, id_synonym_dict=None,
                       id_replacement_mode="partial", domain_style="vlines"):
        """
        id_replacement_mode have to be either partial or exact
        """
        #from Routines import SequenceRoutines
        sequence_number = len(alignment)
        alignment_length = len(alignment[0].seq)

        figure = plt.figure(figsize=(figure_width, sequence_number))
        subplot = plt.subplot(1, 1, 1)

        subplot.get_yaxis().set_visible(False)
        #subplot.get_xaxis().set_visible(False)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')
        subplot.spines['right'].set_color('none')
        subplot.spines['left'].set_color('none')
        subplot.spines['top'].set_color('none')

        protein_height = 10

        dist_between_proteins = 10
        start_x = 0
        start_y = - dist_between_proteins

        gap_line_y_shift = int(protein_height/2)
        gap_line_y_jump = int(protein_height/2)

        domen_colors = []
        for feature in features:
            if (feature.type == "domen") or (feature.type == "domain"):
                domen_colors.append(subplot._get_lines.color_cycle.readline())

        for record in alignment:

            gap_coords_list, gap_len_list = self.find_homopolymers(record.seq, "-", min_size=1,
                                                                   search_type="perfect")
            #print gap_coords_list, gap_len_list

            start_y += protein_height + dist_between_proteins
            gap_y_start = gap_line_y_shift + start_y
            gap_y_jump = gap_y_start + gap_line_y_jump
            prev_x = 0
            """
            figure.text(0, start_y, record.id, rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                         horizontalalignment='center',
                         verticalalignment='center')
            """
            if id_synonym_dict:
                if id_replacement_mode == "exact":
                    if record.id in id_synonym_dict:
                        record_label = id_synonym_dict[record.id]
                    else:
                        record_label = record.id
                        print("WARNING!!! Synonym for %s was not found" % record.id)
                elif id_replacement_mode == "partial":

                    partial_syn_list = []
                    for partial_syn in id_synonym_dict:
                        if partial_syn in record.id:
                            partial_syn_list.append(partial_syn)

                    if len(partial_syn_list) > 1:
                        print("WARNING!!! More than one possible replacement for %s was found: %s. No replacement then." % (record.id, ",".join(partial_syn_list)))
                        record_label = record.id
                    elif not partial_syn_list:
                        record_label = record.id
                        print("WARNING!!! Synonym for %s was not found" % record.id)
                    else:
                        record_label = id_synonym_dict[partial_syn_list[0]]
                else:
                    raise ValueError("Unknown id replacement mode")

            else:
                record_label = record.id

            subplot.annotate(record_label, xy=(0, gap_y_start), xycoords='data', fontsize=16,
                             xytext=(-15, 1.5 * gap_line_y_shift), textcoords='offset points',
                             ha='right', va='top')

            for gap_coords, gap_len in zip(gap_coords_list, gap_len_list):

                if gap_coords[0] != 0:

                    fragment = Rectangle((prev_x, start_y), gap_coords[0] - prev_x, protein_height, fill=False,
                                         edgecolor="black", facecolor="grey")
                    #print prev_x
                    #print gap_coords[0] - prev_x

                    subplot.add_patch(fragment)
                prev_x = gap_coords[1]
                #print [gap_coords[0], gap_coords[0] + int(gap_len/2) + 1, gap_coords[1]]
                plt.plot([gap_coords[0], gap_coords[0] + int(gap_len/2) + 1, gap_coords[1]], #plt.plot([gap_coords[0] + 2, gap_coords[0] + int(gap_len/2) + 1, gap_coords[1] - 1],
                         [gap_y_start, gap_y_jump, gap_y_start], color="black", linewidth=1)

            if not gap_coords_list:
                fragment = Rectangle((prev_x, start_y), alignment_length, protein_height, fill=False,
                                     edgecolor="black", facecolor="grey")
                subplot.add_patch(fragment)
            else:
                if gap_coords_list[-1][-1] != alignment_length:
                    fragment = Rectangle((prev_x, start_y), alignment_length - prev_x, protein_height, fill=False,
                                         edgecolor="black", facecolor="grey")
                    #print prev_x, alignment_length - prev_x
                    subplot.add_patch(fragment)
            i = 0
            for feature in features:
                if (feature.type == "domen") or (feature.type == "domain"):
                    #print feature.id, feature.location
                    if domain_style == "rectangle":
                        fragment = Rectangle((feature.location.start, start_y), len(feature)-1, protein_height, fill=False,
                                             facecolor="grey", edgecolor=domen_colors[i]) #edgecolor="green",
                        subplot.add_patch(fragment)
                    elif domain_style == "vlines":
                        plt.vlines(feature.location.start, protein_height, start_y + protein_height, colors=domen_colors[i])
                        plt.vlines(feature.location.end - 1, protein_height, start_y + protein_height, colors=domen_colors[i])
                    i += 1

        for feature in features:
            if feature.type == "domen":
                #print feature.id, feature.location
                subplot.annotate(feature.id, xy=(feature.location.start + len(feature)/2, gap_y_start + protein_height),
                                 xycoords='data', fontsize=label_fontsize,
                                 xytext=(0, 1.5 * gap_line_y_shift), textcoords='offset points', ha='center', va='top')

        plt.xlim(xmin=0, xmax=alignment_length + 10)
        plt.ylim(ymin=0, ymax=start_y + 2 * protein_height)
        #plt.tight_layout()
        plt.subplots_adjust(left=left_offset, right=0.95)#bottom=0.1, right=0.8, top=0.9)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension))

    def draw_alignment_from_file(self, alignment_file, feature_gff, output_prefix, alignment_style=None,
                                 alignment_format="fasta", ext_list=["svg", "png"], label_fontsize=13,
                                 left_offset=0.2, figure_width=8, id_synonym_dict=None,
                                 id_replacement_mode="partial"):

        alignment = AlignIO.read(alignment_file, format=alignment_format)
        if feature_gff:
            with open(feature_gff, "r") as gff_fd:
                record = list(GFF.parse(gff_fd))[0]
                features = record.features
                record_id = record.id
                gap_coords_list, gap_len_list = self.find_homopolymers(record.seq, "-", min_size=1,
                                                                       search_type="perfect")
        else:
            features = []

        self.draw_alignment(alignment, features, output_prefix, ext_list=ext_list, label_fontsize=label_fontsize,
                            left_offset=left_offset, figure_width=figure_width, id_synonym_dict=id_synonym_dict,
                            id_replacement_mode=id_replacement_mode)

    @staticmethod
    def draw_length_histogram(sequence_dict, output_prefix, number_of_bins=None, width_of_bins=None,
                              min_length=1, max_length=None, extensions=("png", "svg"),
                              legend_location='best'):
        length_dict = SynDict()

        for record in sequence_dict:
            length_dict[record] = len(sequence_dict[record].seq)

        length_dict.write("%s.len" % output_prefix)

        lengths = length_dict.values()

        max_len = max(lengths)
        min_len = min(lengths)
        median = np.median(lengths)
        mean = np.mean(lengths)

        if max_length is None:
            maximum_length = max_len
        else:
            maximum_length = max_length

        filtered = []

        if (maximum_length < max_len) and (min_length > 1):
            for entry in lengths:
                if min_length <= entry <= maximum_length:
                    filtered.append(entry)
        elif min_length > 1:
            for entry in lengths:
                if min_length <= entry:
                    filtered.append(entry)
        elif maximum_length < max_len:
            for entry in lengths:
                if entry <= maximum_length:
                    filtered.append(entry)
        else:
            filtered = lengths

        plt.figure(1, figsize=(6, 6))
        plt.subplot(1, 1, 1)

        if number_of_bins:
            bins = number_of_bins
        elif width_of_bins:
            bins = np.arange(min_length - 1, maximum_length, width_of_bins, dtype=np.int32)
            bins[0] += 1
            bins = np.append(bins, [maximum_length])
        else:
            bins = 30
        plt.hist(filtered, bins=bins)
        plt.xlim(xmin=min_length, xmax=maximum_length)
        plt.xlabel("Length")
        plt.ylabel("N")
        plt.title("Distribution of sequence lengths")
        plt.legend(("Min: %i\nMax: %i\nMean: %i\nMedian: %i" % (min_len, max_len, mean, median),), loc=legend_location)
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    def draw_heatmap_and_three_percent_histograms(self, first_histo_values, second_histo_values,
                                                  third_histo_values, output_prefix, figsize=(12, 12),
                                                  extensions=("png", "svg"), stats_as_legend=True):
        """
        second_histo_values and third_histo_values are used to build heatmap
        """

        plt.figure(1, figsize=figsize)

        for (index, histo_values, title) in zip([1, 2, 3],
                                                [first_histo_values, second_histo_values, third_histo_values],
                                                ["Total support", "CDS support", "Intron support"]):
            subplot = plt.subplot(2, 2, index)

            self.percent_histogram(histo_values, output_prefix=None, n_bins=20, title=title, xlabel="%",
                                   ylabel="Number of transcripts", label=None, extensions=("png", "svg"),
                                   legend=None, legend_location="best", input_mode="percent", xmax=None,
                                   xmin=None, stats_as_legend=stats_as_legend)

        bins = np.linspace(0, 100, 21)

        subplot = plt.subplot(2, 2, 4)
        #print bins
        counts, xedges, yedges, image = plt.hist2d(second_histo_values,
                                                   third_histo_values,
                                                   bins=(bins, bins),
                                                   range=[[0, 100], [0, 100]])
        max_counts = int(np.nanmax(counts))

        cmap = plt.get_cmap('jet', max_counts)
        #cmap.set_under('gray')
        mappable = plt.cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(0.00001, max_counts)
        #mappable.set_array([])
        #mappable.set_clim(-0.5, ncolors+0.5)
        colorbar = plt.colorbar(mappable)
        plt.xlabel("CDS support")
        plt.ylabel("Intron support")
        plt.title("Transcript support")
        plt.tight_layout()
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    @staticmethod
    def draw_plot(input_file, output_prefix, x_column_index=0, y_column_index=1, separator="\t", min_x=None,
                  max_x=None, min_y=None, max_y=None, extensions=["png", "svg"], xlabel=None, ylabel=None,
                  title=None, width=6, height=6, markersize=2, ylogbase=10, type="plot", grid=False, 
                  correlation=False, close_plot=True):

        data = np.loadtxt(input_file, comments="#", usecols=(x_column_index, y_column_index), delimiter=separator)
        plt.figure(1, figsize=(width, height), dpi=300)
        plt.subplot(1, 1, 1)
        if type == "plot":
            plt.plot(data[:, 0], data[:, 1], markersize=markersize)
        elif type == "scatter":
            plt.scatter(data[:, 0], data[:, 1], s=markersize)
        plt.xlim(xmin=min_x, xmax=max_x)
        plt.ylim(ymin=min_y, ymax=max_y)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)
        if grid:
            plt.grid()
        if correlation:
            print("Kendal's tau")
            print(stats.kendalltau(data[:, 0], data[:, 1]))

            print("Pearson's r")
            print(stats.pearsonr(data[:, 0], data[:, 1]))
        for ext in extensions:
            plt.savefig("%s.%s.%s" % (output_prefix, type, ext))

        plt.yscale("log")

        for ext in extensions:
            plt.savefig("%s.%s.ylog%i.%s" % (output_prefix, type, ylogbase, ext))

        if close_plot:
            plt.close()

    @staticmethod
    def draw_precalculated_heatmap(heatmap_array, output_prefix=None, figsize=(5, 5), extensions=("png", "svg")):

        if output_prefix:
            figure = plt.figure(1, figsize=figsize)

        heatmap = plt.imshow(heatmap_array, origin='low', interpolation='none')

        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

    def draw_histogram_from_multiple_files(self, filelist, output_prefix,  filelabels=None, nbins=30,
                                           figsize=(5, 5), title=None, xlabel=None, ylabel=None,
                                           extensions=("png", "svg"), separator="\n", xmin=None, xmax=None,
                                           histtype="stepfilled"):

        colors = ["aqua", "blue", "black", "brown",
                  "red", "fuchsia", "maroon", "orange",
                  "purple", "yellow", "sienna", "green",
                  "olive"]

        if filelabels:
            labels = filelabels
        else:
            labels = []
            for filename in filelist:
                path, basename, extention = self.split_filename(filename)
                labels.append(basename)

        if (xmin is not None) and (xmax is not None):
            bins = np.linspace(xmin, xmax, nbins + 1)
            bin_middles = []
            for i in range(0, len(bins) - 1):
                bin_middles.append((bins[i] + bins[i+1])/2)
        else:
            bins = nbins

        data = []
        for filename, label in zip(filelist, labels):
            filedata = np.fromfile(filename, sep=separator)
            data.append(filedata)

            print("%s\tmin %f\t max %f" % (label, np.min(filedata), np.max(filedata)))

            # a = np.histogram(filedata, bins=bins)
            #print a

        #print xmin
        #print xmax
        #print bins

        figure = plt.figure(1, figsize=figsize)
        ax = plt.axes()
        ax.set_color_cycle(colors)

        #colors = plt.cm.jet(np.linspace(0.2, 1, len(filelist)))

        plt.hist(data, label=labels, bins=bins, histtype=histtype) #, color=colors)

        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)

        plt.xlim(xmin=xmin, xmax=xmax)
        plt.legend(loc="best")

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    def test_colormap(self, colormap_name, color_number=10, output_prefix=None, extension_list=("png",), side=10):
        cmap = plt.get_cmap(colormap_name, color_number)
        figure = plt.figure(figsize=(4, 4))
        subplot = plt.subplot(1, 1, 1)
        subplot.get_yaxis().set_visible(False)
        subplot.get_xaxis().set_visible(False)

        posx = 0
        posy = 0
        for i in range(0, color_number):
            fragment = Rectangle((posx, posy), side, side, fill=True, edgecolor="black",
                                 facecolor=cmap(i), linewidth=0.5)
            posy += side
            subplot.add_patch(fragment)
        plt.ylim(ymin=0, ymax=posy)
        plt.title("Colormap %s" % colormap_name)

        if output_prefix:
            for ext in extension_list:
                plt.savefig("%s.%s" % (output_prefix, ext))
        else:
            plt.show()

    def draw_dot_plot_per_scaffold_from_last_self_alignment(self, last_collection,
                                                            output_prefix=None,
                                                            extension_list=("png", ),
                                                            scaffold_black_list=(), scaffold_white_list=(),
                                                            scaffold_reverse_list=(),
                                                            figsize=(16, 16), dpi=300,
                                                            grid_color='black',
                                                            bar_color='grey',
                                                            same_strand_color='red',
                                                            diff_strand_color='blue',
                                                            title=None,
                                                            label=None,
                                                            gridwidth=1,
                                                            show_grid=True,
                                                            linewidth=0.01,
                                                            scaffold_label_fontsize=13,
                                                            axes_label_fontstyle="italic",
                                                            axes_label_weight="normal",
                                                            axes_label_fontsize=15,
                                                            axes_label_distance=6,
                                                            antialiased_lines=None,
                                                            target_scaffold_labels_angle=45,
                                                            query_scaffold_labels_angle=0,
                                                            show_labels=True,
                                                            bottom_offset=0.1,
                                                            top_offset=0.9,
                                                            left_offset=0.1,
                                                            right_offset=0.9,
                                                            x_axis_visible=False,
                                                            y_axis_visible=False,
                                                            x_labelpad=None,
                                                            y_labelpad=None
                                                            ):

        for scaffold in scaffold_white_list:
            self.draw_dot_plot_from_last_alignment(last_collection,
                                              output_prefix="%s.%s" % (output_prefix, scaffold),
                                              extension_list=extension_list,
                                              target_black_list=scaffold_black_list, target_white_list=[scaffold,],
                                              target_ordered_list=[scaffold,], target_reverse_list=scaffold_reverse_list,
                                              query_black_list=scaffold_black_list, query_white_list=[scaffold,],
                                              query_ordered_list=[scaffold,], query_reverse_list=scaffold_reverse_list,
                                              figsize=figsize, dpi=dpi,
                                              grid_color=grid_color,
                                              bar_color=bar_color,
                                              same_strand_color=same_strand_color,
                                              diff_strand_color=diff_strand_color,
                                              title=title,
                                              target_label=label,
                                              query_label=label,
                                              gridwidth=gridwidth,
                                              show_grid=show_grid,
                                              linewidth=linewidth,
                                              scaffold_label_fontsize=scaffold_label_fontsize,
                                              axes_label_fontstyle=axes_label_fontstyle,
                                              axes_label_weight=axes_label_weight,
                                              axes_label_fontsize=axes_label_fontsize,
                                              axes_label_distance=axes_label_distance,
                                              antialiased_lines=antialiased_lines,
                                              target_scaffold_labels_angle=target_scaffold_labels_angle,
                                              query_scaffold_labels_angle=query_scaffold_labels_angle,
                                              show_target_labels=show_labels,
                                              show_query_labels=show_labels,
                                              bottom_offset=bottom_offset,
                                              top_offset=top_offset,
                                              left_offset=left_offset,
                                              right_offset=right_offset,
                                              x_axis_visible=x_axis_visible,
                                              y_axis_visible=y_axis_visible,
                                              x_labelpad=x_labelpad,
                                              y_labelpad=y_labelpad
                                                   )

    def draw_dot_plot_per_scaffold_vs_all_from_last_self_alignment(self, last_collection,
                                                            output_prefix=None,
                                                            extension_list=("png", ),
                                                            scaffold_black_list=(), scaffold_white_list=(),
                                                            scaffold_reverse_list=(),
                                                            ordered_list=(),
                                                            figsize=(16, 16), dpi=300,
                                                            grid_color='black',
                                                            bar_color='grey',
                                                            same_strand_color='red',
                                                            diff_strand_color='blue',
                                                            title=None,
                                                            label=None,
                                                            gridwidth=1,
                                                            show_grid=True,
                                                            linewidth=0.01,
                                                            scaffold_label_fontsize=13,
                                                            axes_label_fontstyle="italic",
                                                            axes_label_weight="normal",
                                                            axes_label_fontsize=15,
                                                            axes_label_distance=6,
                                                            antialiased_lines=None,
                                                            target_scaffold_labels_angle=45,
                                                            query_scaffold_labels_angle=0,
                                                            show_labels=True,
                                                            bottom_offset=0.1,
                                                            top_offset=0.9,
                                                            left_offset=0.1,
                                                            right_offset=0.9,
                                                            x_axis_visible=False,
                                                            y_axis_visible=False,
                                                            x_labelpad=None,
                                                            y_labelpad=None
                                                           ):

        for scaffold in scaffold_white_list:
            self.draw_dot_plot_from_last_alignment(last_collection,
                                                   output_prefix="%s.%s" % (output_prefix, scaffold),
                                                   extension_list=extension_list,
                                                   target_black_list=scaffold_black_list, target_white_list=scaffold_white_list,
                                                   target_ordered_list=ordered_list, target_reverse_list=scaffold_reverse_list,
                                                   query_black_list=scaffold_black_list, query_white_list=[scaffold,],
                                                   query_ordered_list=[scaffold,], query_reverse_list=scaffold_reverse_list,
                                                   figsize=figsize, dpi=dpi,
                                                   grid_color=grid_color,
                                                   bar_color=bar_color,
                                                   same_strand_color=same_strand_color,
                                                   diff_strand_color=diff_strand_color,
                                                   title=title,
                                                   target_label=label,
                                                   query_label=label,
                                                   gridwidth=gridwidth,
                                                   show_grid=show_grid,
                                                   linewidth=linewidth,
                                                   scaffold_label_fontsize=scaffold_label_fontsize,
                                                   axes_label_fontstyle=axes_label_fontstyle,
                                                   axes_label_weight=axes_label_weight,
                                                   axes_label_fontsize=axes_label_fontsize,
                                                   axes_label_distance=axes_label_distance,
                                                   antialiased_lines=antialiased_lines,
                                                   target_scaffold_labels_angle=target_scaffold_labels_angle,
                                                   query_scaffold_labels_angle=query_scaffold_labels_angle,
                                                   show_target_labels=show_labels,
                                                   show_query_labels=show_labels,
                                                   bottom_offset=bottom_offset,
                                                   top_offset=top_offset,
                                                   left_offset=left_offset,
                                                   right_offset=right_offset,
                                                   x_axis_visible=x_axis_visible,
                                                   y_axis_visible=y_axis_visible,
                                                   x_labelpad=x_labelpad,
                                                   y_labelpad=y_labelpad
                                                   )

    def draw_dot_plot_from_last_alignment(self, last_collection,
                                          output_prefix=None,
                                          extension_list=("png", ),
                                          target_black_list=(), target_white_list=(),
                                          target_ordered_list=(), target_reverse_list=(),
                                          query_black_list=(), query_white_list=(),
                                          query_ordered_list=(), query_reverse_list=(),
                                          remove_scaffolds_absent_in_target_ordered_list=False,
                                          remove_scaffolds_absent_in_query_ordered_list=False,
                                          figsize=(16, 16), dpi=300,
                                          grid_color='black',
                                          bar_color='grey',
                                          same_strand_color='red',
                                          diff_strand_color='blue',
                                          title=None,
                                          target_label=None,
                                          query_label=None,
                                          gridwidth=1,
                                          show_grid=True,
                                          linewidth=0.01,
                                          scaffold_label_fontsize=13,
                                          axes_label_fontstyle="italic",
                                          axes_label_weight="normal",
                                          axes_label_fontsize=15,
                                          axes_label_distance=6,
                                          antialiased_lines=None,
                                          target_scaffold_labels_angle=45,
                                          query_scaffold_labels_angle=0,
                                          show_target_labels=True,
                                          show_query_labels=True,
                                          bottom_offset=0.1,
                                          top_offset=0.9,
                                          left_offset=0.1,
                                          right_offset=0.9,
                                          x_axis_visible=False,
                                          y_axis_visible=False,
                                          x_labelpad=None,
                                          y_labelpad=None,
                                          show_length_ticks=False,
                                          tick_step=10000000,
                                          tick_unit=1000000):

        target_scaffold_list = self.get_filtered_scaffold_list(last_collection.target_scaffold_list,
                                                               scaffold_black_list=target_black_list,
                                                               sort_scaffolds=False,
                                                               scaffold_ordered_list=target_ordered_list,
                                                               scaffold_white_list=target_white_list,
                                                               remove_scaffolds_absent_in_ordered_list=remove_scaffolds_absent_in_target_ordered_list)

        query_scaffold_list = self.get_filtered_scaffold_list(last_collection.query_scaffold_list,
                                                              scaffold_black_list=query_black_list,
                                                              sort_scaffolds=False,
                                                              scaffold_ordered_list=query_ordered_list,
                                                              scaffold_white_list=query_white_list,
                                                              remove_scaffolds_absent_in_ordered_list=remove_scaffolds_absent_in_query_ordered_list)

        target_length_df = last_collection.target_scaffold_lengths.loc[target_scaffold_list]
        target_length_df["cum_end"] = target_length_df["length"].cumsum()
        target_length_df["cum_start"] = target_length_df["cum_end"] - target_length_df["length"]

        total_target_len = target_length_df["length"].sum()

        query_length_df = last_collection.query_scaffold_lengths.loc[query_scaffold_list]
        query_length_df["cum_end"] = query_length_df["length"].cumsum()
        query_length_df["cum_start"] = query_length_df["cum_end"] - query_length_df["length"]

        total_query_len = query_length_df["length"].sum()

        bar_width_fraction = 0.02
        bar_width = int(max(total_query_len, total_target_len) * bar_width_fraction)

        print("%s\tDrawing..." % str(datetime.datetime.now()))
        print("%s\t\tInitializing figure..." % str(datetime.datetime.now()))

        figure = plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.subplot(1, 1, 1)

        print("%s\t\tInitializing figure finished..." % str(datetime.datetime.now()))
        print("%s\t\tDrawing grid..." % str(datetime.datetime.now()))

        ax.add_patch(Rectangle((0, total_query_len), total_target_len, bar_width, color=bar_color))  # top bar
        ax.add_patch(Rectangle((0, -bar_width), total_target_len, bar_width, color=bar_color))       # bottom bar
        ax.add_patch(Rectangle((-bar_width, 0), bar_width, total_query_len, color=bar_color))        # left bar
        ax.add_patch(Rectangle((total_target_len, 0), bar_width, total_query_len, color=bar_color))  # right bar

        """
        @staticmethod
        def create_tick_formatter_function(max_value, tick_type="nucleotide"):
            max_val = max_value * 1.1
            if tick_type == "nucleotide":
                if max_val // (10 ** 9) > 2:
                    def tick_formater(x, pos):
                        return '%1.1f Gbp' % (x * 1e-9)
                elif max_val // (10 ** 6) > 200:
                    def tick_formater(x, pos):
                        return '%.0f Mbp' % (x * 1e-6)
                elif max_val // (10 ** 6) > 2:
                    def tick_formater(x, pos):
                        return '%.1f Mbp' % (x * 1e-6)
                elif max_val // (10 ** 3) > 2:
                    def tick_formater(x, pos):
                        return '%.1f kbp' % (x * 1e-3)
                else:
                    def tick_formater(x, pos):
                        return '%i bp' % (int(x))

                return FuncFormatter(tick_formater)

            else:
                raise ValueError("ERROR!!! Tick formter for %s is not implemented yet!" % tick_type)
        """

        if show_grid:
            for query_cum_start in query_length_df["cum_start"]:
                ax.add_line(Line2D((-bar_width, total_target_len+bar_width), (query_cum_start, query_cum_start),
                                   color=grid_color, linewidth=gridwidth))

            ax.add_line(Line2D((-bar_width, total_target_len+bar_width), (total_query_len, total_query_len),
                               color=grid_color, linewidth=gridwidth))

            for target_cum_start in target_length_df["cum_start"]:
                ax.add_line(Line2D((target_cum_start, target_cum_start), (-bar_width, total_query_len + bar_width),
                                   color=grid_color, linewidth=gridwidth))
            ax.add_line(Line2D((total_target_len, total_target_len), (-bar_width, total_query_len + bar_width),
                               color=grid_color, linewidth=gridwidth))
        else:
            for query_cum_start in query_length_df["cum_start"]:
                ax.add_line(Line2D((-bar_width, 0), (query_cum_start, query_cum_start),
                                   color=grid_color, linewidth=gridwidth))
                ax.add_line(Line2D((total_target_len, total_target_len + bar_width), (query_cum_start, query_cum_start),
                                   color=grid_color, linewidth=gridwidth))

            ax.add_line(Line2D((-bar_width, 0), (total_query_len, total_query_len),
                               color=grid_color, linewidth=gridwidth))
            ax.add_line(Line2D((total_target_len, total_target_len + bar_width), (total_query_len, total_query_len),
                               color=grid_color, linewidth=gridwidth))

            for target_cum_start in target_length_df["cum_start"]:
                ax.add_line(Line2D((target_cum_start, target_cum_start), (-bar_width, 0),
                                   color=grid_color, linewidth=gridwidth))
                ax.add_line(Line2D((target_cum_start, target_cum_start), (total_query_len, total_query_len + bar_width),
                                   color=grid_color, linewidth=gridwidth))

            ax.add_line(Line2D((total_target_len, total_target_len), (-bar_width, 0),
                               color=grid_color, linewidth=gridwidth))
            ax.add_line(Line2D((total_target_len, total_target_len), (total_query_len, total_query_len + bar_width),
                               color=grid_color, linewidth=gridwidth))

        if show_length_ticks:
            for target_scaffold in target_length_df.index:
                tick_labels = np.arange(0, target_length_df.loc[target_scaffold, "length"], tick_step)
                tick_list = list(tick_labels + target_length_df.loc[target_scaffold, "cum_start"])
                tick_labels = list(map(str, tick_labels // tick_unit))

                if len(tick_list) > 1:
                    for tick in tick_list[1:]:
                        ax.add_line(Line2D((tick, tick), (-bar_width/4, 0),
                                           color="black", linewidth=gridwidth / 4))
                    if len(tick_list) >= 5:
                        for tick, tick_label in zip(tick_list[::5][1:], tick_labels[::5][1:]):
                            ax.add_line(Line2D((tick, tick), (-bar_width/2, 0),
                                               color="red", linewidth=gridwidth / 2))
                            ax.text(tick, -bar_width * 0.6, tick_label,
                                    fontsize=5, #scaffold_label_fontsize/3,
                                    horizontalalignment='center',
                                    verticalalignment='top', )
            for query_scaffold in query_length_df.index:
                tick_labels = np.arange(0, query_length_df.loc[query_scaffold, "length"], tick_step)
                tick_list = list(tick_labels + query_length_df.loc[query_scaffold, "cum_start"])
                tick_labels = list(map(str, tick_labels // tick_unit))

                if len(tick_list) > 1:
                    for tick in tick_list[1:]:
                        ax.add_line(Line2D((-bar_width/4, 0), (tick, tick),
                                           color="black", linewidth=gridwidth / 4))
                    if len(tick_list) >= 5:
                        for tick, tick_label in zip(tick_list[::5][1:], tick_labels[::5][1:]):
                            ax.add_line(Line2D((-bar_width/2, 0), (tick, tick),
                                               color="red", linewidth=gridwidth / 2))
                            ax.text(-bar_width * 0.6, tick, tick_label,
                                    fontsize=5, #scaffold_label_fontsize/3,
                                    horizontalalignment='right',
                                    verticalalignment='center', )

        print("%s\t\tDrawing grid finished..." % str(datetime.datetime.now()))
        print("%s\t\tAdding labels..." % str(datetime.datetime.now()))

        if show_target_labels:
            for target_scaffold_id in target_scaffold_list:
                ax.text((target_length_df.loc[target_scaffold_id]["cum_start"] + target_length_df.loc[target_scaffold_id]["cum_end"])/2,
                        total_query_len + 1.5 * bar_width, target_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=target_scaffold_labels_angle,
                        horizontalalignment='left',
                        verticalalignment='bottom',)
                ax.text((target_length_df.loc[target_scaffold_id]["cum_start"] + target_length_df.loc[target_scaffold_id]["cum_end"])/2,
                        -1.5 * bar_width, target_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=target_scaffold_labels_angle,
                        horizontalalignment='right',
                        verticalalignment='top',)
        if show_query_labels:
            for query_scaffold_id in query_scaffold_list:
                ax.text(total_target_len + 1.5 * bar_width,
                        (query_length_df.loc[query_scaffold_id]["cum_start"] + query_length_df.loc[query_scaffold_id]["cum_end"])/2,
                         query_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=query_scaffold_labels_angle,
                        horizontalalignment='left',
                        verticalalignment='center' if query_scaffold_labels_angle == 0 else 'bottom')
                ax.text(-1.5 * bar_width,
                        (query_length_df.loc[query_scaffold_id]["cum_start"] + query_length_df.loc[query_scaffold_id]["cum_end"])/2,
                         query_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=query_scaffold_labels_angle,
                        horizontalalignment='right',
                        verticalalignment='center' if query_scaffold_labels_angle == 0 else 'top')

        if title:
            plt.title(title)
        if target_label:
            ax.text(total_target_len/2,
                    -bar_width*axes_label_distance,
                    target_label,
                    fontsize=axes_label_fontsize,
                    fontstyle=axes_label_fontstyle,
                    fontweight=axes_label_weight,
                    horizontalalignment='center',
                    verticalalignment='center')
        if query_label:
            ax.text(-bar_width*axes_label_distance,
                    total_query_len/2,
                    query_label,
                    fontsize=axes_label_fontsize,
                    fontstyle=axes_label_fontstyle,
                    fontweight=axes_label_weight,
                    rotation=90,
                    horizontalalignment='center',
                    verticalalignment='center')

        plt.xlim(xmin=-bar_width * 2, xmax=total_target_len + 2 * bar_width)
        plt.ylim(ymin=-bar_width * 2, ymax=total_query_len + 2 * bar_width)

        if not x_axis_visible:
            ax.spines['bottom'].set_color('none')
            ax.spines['top'].set_color('none')

        if not y_axis_visible:
            ax.spines['right'].set_color('none')
            ax.spines['left'].set_color('none')

        ax.get_yaxis().set_visible(y_axis_visible)
        ax.get_xaxis().set_visible(x_axis_visible)

        print("%s\t\tAdding labels finished..." % str(datetime.datetime.now()))
        print("%s\t\tDrawing alignments..." % str(datetime.datetime.now()))

        if isinstance(last_collection, CollectionLast):
            def get_strand_specific_records(collection, target_scaffold_id, query_scaffold_id):
                return (collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                            & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                            & (collection.records[collection.query_strand_syn] == collection.records[collection.target_strand_syn])],
                        collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                            & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                            & (collection.records[collection.query_strand_syn] != collection.records[collection.target_strand_syn])])

        elif isinstance(last_collection, CollectionPSL):
            def get_strand_specific_records(collection, target_scaffold_id, query_scaffold_id):
                return (collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                           & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                           & (collection.records[collection.query_strand_syn] == "+")],
                        collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                           & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                           & (collection.records[collection.query_strand_syn] == "-")])
        else:
            raise ValueError("ERROR!!! Unknown collection type (neither CollectionPSL nor CollectionLast)")

        def line_segments_generator(dataframe):
            for row_tuple in dataframe.itertuples(index=False):
                yield (row_tuple[:2], row_tuple[2:])

        for query_scaffold_id in query_scaffold_list:
            for target_scaffold_id in target_scaffold_list:

                same_strand_records_df, diff_strand_records_df = get_strand_specific_records(last_collection, target_scaffold_id, query_scaffold_id)
                """
                same_strand_records_df = \
                    last_collection.records[last_collection.records["target_id"].isin([target_scaffold_id])
                                            & last_collection.records["query_id"].isin([query_scaffold_id])
                                            & (last_collection.records["query_strand"] == last_collection.records["target_strand"])]


                
                diff_strand_records_df = \
                    last_collection.records[last_collection.records["target_id"].isin([target_scaffold_id])
                                            & last_collection.records["query_id"].isin([query_scaffold_id])
                                            & (last_collection.records["query_strand"] != last_collection.records["target_strand"])]
                """
                if not same_strand_records_df.empty:
                    data = pd.DataFrame()
                    """
                    data["x1"] = same_strand_records["target_start"] + target_length_df.loc[target_scaffold_id]["cum_start"]
                    data["y1"] = same_strand_records["query_start"] + query_length_df.loc[query_scaffold_id]["cum_start"]

                    data["x2"] = data["x1"] + same_strand_records["target_hit_len"] - 1
                    data["y2"] = data["y1"] + same_strand_records["query_hit_len"] - 1
                    """
                    data["x1"] = same_strand_records_df[last_collection.target_start_syn] + target_length_df.loc[target_scaffold_id]["cum_start"]
                    data["y1"] = same_strand_records_df[last_collection.query_start_syn] + query_length_df.loc[query_scaffold_id]["cum_start"]
                    if isinstance(last_collection, CollectionLast):
                        data["x2"] = data["x1"] + same_strand_records_df[last_collection.target_hit_len_syn] - 1
                        data["y2"] = data["y1"] + same_strand_records_df[last_collection.query_hit_len_syn] - 1
                    elif isinstance(last_collection, CollectionPSL):
                        data["x2"] = same_strand_records_df[last_collection.target_end_syn] + target_length_df.loc[target_scaffold_id]["cum_start"] - 1
                        data["y2"] = same_strand_records_df[last_collection.query_end_syn] + query_length_df.loc[query_scaffold_id]["cum_start"] - 1

                    lines = LineCollection(line_segments_generator(data), colors=same_strand_color, linestyle='solid',
                                           linewidths=linewidth, antialiased=antialiased_lines)

                    ax.add_collection(lines)

                if not diff_strand_records_df.empty:
                    data = pd.DataFrame()
                    data["x1"] = diff_strand_records_df[last_collection.target_start_syn] + target_length_df.loc[target_scaffold_id]["cum_start"]
                    if isinstance(last_collection, CollectionLast):
                        # in last tab format coordinates for minus strand start from the end of sequence
                        data["y1"] = query_length_df.loc[query_scaffold_id]["length"] - diff_strand_records_df[last_collection.query_start_syn] + query_length_df.loc[query_scaffold_id]["cum_start"]
                        data["x2"] = data["x1"] + diff_strand_records_df[last_collection.target_hit_len_syn] - 1
                        data["y2"] = data["y1"] - diff_strand_records_df[last_collection.query_hit_len_syn] + 1

                    elif isinstance(last_collection, CollectionPSL):
                        # in PSL format coordinates for minus strand start from the start of sequence
                        data["y1"] = diff_strand_records_df[last_collection.query_end_syn] + query_length_df.loc[query_scaffold_id]["cum_start"] - 1
                        data["x2"] = diff_strand_records_df[last_collection.target_end_syn] + target_length_df.loc[target_scaffold_id]["cum_start"] - 1
                        data["y2"] = diff_strand_records_df[last_collection.query_start_syn] + query_length_df.loc[query_scaffold_id]["cum_start"] - 1
                    #print(data)
                    lines = LineCollection(line_segments_generator(data), colors=diff_strand_color, linestyle='solid',
                                           linewidths=linewidth, antialiased=antialiased_lines)
                    #print(data)
                    #print(diff_strand_records_df[["target_start", "target_hit_len", "query_start", "query_hit_len"]])
                    ax.add_collection(lines)

        print("%s\t\tDrawng alignments finished..." % str(datetime.datetime.now()))
        print("%s\tDrawing finished..." % str(datetime.datetime.now()))

        plt.subplots_adjust(left=left_offset, bottom=bottom_offset, right=right_offset, top=top_offset)

        if output_prefix:
            print("%s\tWriting to file..." % str(datetime.datetime.now()))
            for extension in extension_list:
                plt.savefig("%s.%s" % (output_prefix, extension))

            print("%s\tWriting to file finished..." % str(datetime.datetime.now()))

    def draw_venn(self, set_list, label_list, output_prefix, extensions=("png", "svg"), title=None):
        if sys.version_info[0] == 2:
            raise ImportError("Package venn doesn't work in Python 2")
        else:
            #try:
            count_dict = OrderedDict(dict([(sample, set) for sample, set in zip(label_list, set_list)]))
            subplot = venn.venn(count_dict)
            if title:
                plt.title(title)
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

            return subplot
            #except NameError:
            #    pass

    def draw_venn_from_files(self, file_list, label_list, output_prefix, extensions=("png", "svg"), title=None):
        if sys.version_info[0] == 2:
            raise ImportError("Package venn doesn't work in Python 2")
        else:
            set_list = list(map(set, [pd.read_csv(filename, sep="\t", squeeze=True) for filename in file_list]))

            return self.draw_venn(set_list, label_list, output_prefix, extensions=extensions, title=title)

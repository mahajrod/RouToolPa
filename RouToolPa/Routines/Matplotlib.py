__author__ = 'mahajrod'
import os
import sys

from itertools import cycle, islice
from collections.abc import Iterable
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, BboxConnectorPatch
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
if sys.version_info[0] == 3:
    from venn import venn
elif sys.version_info[0] == 2:
    import pyvenn as venn


import numpy as np

from RouToolPa.Collections.General import IdList, SynDict, IdSet


class MatplotlibRoutines:
    def __init__(self):
        pass

    @staticmethod
    def connect_bbox(bbox1, bbox2,
                     loc1a, loc2a, loc1b, loc2b,
                     prop_lines, prop_patches=None):
        if prop_patches is None:
            prop_patches = prop_lines.copy()
            prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

        c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
        c1.set_clip_on(False)
        c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
        c2.set_clip_on(False)

        bbox_patch1 = BboxPatch(bbox1, **prop_patches)
        bbox_patch2 = BboxPatch(bbox2, **prop_patches)

        p = BboxConnectorPatch(bbox1, bbox2,
                               loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                               **prop_patches)
        p.set_clip_on(False)

        return c1, c2, bbox_patch1, bbox_patch2, p

    def zoom_effect(self, ax1, ax2, xmin, xmax, alpha=0.22, color="gray", **kwargs):
        """
        ax1 : the main axes
        ax2 : the zoomed axes
        (xmin,xmax) : the limits of the colored area in both plot axes.

        connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
        be marked.  The keywords parameters will be used ti create
        patches.

        """

        trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
        trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

        bbox = Bbox.from_extents(xmin, 0, xmax, 1)

        mybbox1 = TransformedBbox(bbox, trans1)
        mybbox2 = TransformedBbox(bbox, trans2)

        prop_patches = kwargs.copy()
        prop_patches["ec"] = "none"
        prop_patches["alpha"] = alpha
        prop_patches["color"] = color

        c1, c2, bbox_patch1, bbox_patch2, p = self.connect_bbox(mybbox1, mybbox2,
                                                                loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                                                                prop_lines=kwargs, prop_patches=prop_patches)

        ax1.add_patch(bbox_patch1)
        ax2.add_patch(bbox_patch2)
        ax2.add_patch(c1)
        ax2.add_patch(c2)
        ax2.add_patch(p)

        return c1, c2, bbox_patch1, bbox_patch2, p

    @staticmethod
    def annotated_heatmap(data, row_labels, column_labels, subplot=None, colorbar_kw={}, colorbar_label="",
                          title=None, xlabel=None, ylabel=None,
                          horizontal_ticks_angle=0, **kwargs):
        """
        Create a heatmap from a numpy array and two lists of labels.

        Arguments:
            data       : A 2D numpy array of shape (N,M)
            row_labels : A list or array of length N with the labels
                         for the rows
            col_labels : A list or array of length M with the labels
                         for the columns
        Optional arguments:
            ax         : A matplotlib.axes.Axes instance to which the heatmap
                         is plotted. If not provided, use current axes or
                         create a new one.
            cbar_kw    : A dictionary with arguments to
                         :meth:`matplotlib.Figure.colorbar`.
            cbarlabel  : The label for the colorbar
        All other arguments are directly passed on to the imshow call.
        """

        if not subplot:
            subplot = plt.gca()

        # Plot the heatmap
        image = subplot.imshow(data, **kwargs)

        # Create colorbar
        colorbar = subplot.figure.colorbar(image, ax=subplot, **colorbar_kw)
        colorbar.ax.set_ylabel(colorbar_label, rotation=-90, va="bottom")

        # We want to show all ticks...
        subplot.set_xticks(np.arange(data.shape[1]))
        subplot.set_yticks(np.arange(data.shape[0]))
        # ... and label them with the respective list entries.
        subplot.set_xticklabels(column_labels)
        subplot.set_yticklabels(row_labels)

        # Let the horizontal axes labeling appear on top.
        #subplot.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

        # Rotate the tick labels and set their alignment.
        plt.setp(subplot.get_xticklabels(), rotation=horizontal_ticks_angle, ha="right",
                 rotation_mode="anchor")

        # Turn spines off and create white grid.
        for edge, spine in subplot.spines.items():
            spine.set_visible(False)

        subplot.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        subplot.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        subplot.grid(which="minor", color="w", linestyle='-', linewidth=3)
        subplot.tick_params(which="minor", bottom=False, left=False)

        if title is not None:
            plt.title(title)
        if xlabel is not None:
            plt.xlabel(xlabel)
        if ylabel is not None:
            plt.ylabel(ylabel)

        return image, colorbar

    @staticmethod
    def percent_histogram(data, output_prefix=None, n_bins=20, title="", xlabel="%", ylabel="Number", label=None,
                          extensions=("png", "svg"), legend=None, legend_location="best", input_mode="percent", xmax=None,
                          xmin=None, stats_as_legend=True):
        if output_prefix:
            figure = plt.figure()
            subplot = plt.subplot(1, 1, 1)

        if isinstance(n_bins, Iterable):
            number_of_bins = n_bins
        else:
            if input_mode == "percent":
                number_of_bins = np.linspace(0, 100, n_bins+1)
            elif input_mode == "fraction":
                number_of_bins = np.linspace(0, 1.0, n_bins+1)

        n, bins, patches = plt.hist(data, bins=number_of_bins, label=label)

        #print n
        #print bins
        #print patches
        if input_mode == "percent":
            plt.xlim(xmin=0, xmax=100)
        elif input_mode == "fraction":
            plt.xlim(xmin=0, xmax=1)
        else:
            raise ValueError("Unrecognized type of input data(neither percents nor fractions)")

        if title:
            plt.title(title)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)

        if stats_as_legend:
            legenda = "Total: %i\nMean: %.2f %%\nMedian: %.2f %%" if input_mode == "percent" else "Total: %i\nMean: %.2f\nMedian: %.2f"
            legenda = legenda % (len(data), np.mean(data), np.median(data)) if stats_as_legend else legend
            plt.legend((legenda,), loc=legend_location)
        else:
            if legend:
                plt.legend((legend,), loc=legend_location)
        if label:
            plt.legend(loc=legend_location)

        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

    def percent_histogram_from_file(self, data_file, output_prefix, data_type=float, column_list=None, separator=None,
                                    comments="#", n_bins=20, title="", xlabel="%", ylabel="Number",
                                    extensions=("png", "svg"), legend=None, legend_location="best",
                                    stats_as_legend=False, input_mode="percent", label=None):
        #print column_list
        data = np.loadtxt(data_file, dtype=data_type, comments=comments, delimiter=separator, usecols=column_list)
        if input_mode == "percent":
            n_bins = np.linspace(0, 100, n_bins+1)
        elif input_mode == "fraction":
            n_bins = np.linspace(0, 1.0, n_bins+1)
        else:
            raise ValueError("Unrecognized type of input data(neither percents nor fractions)")
        #legenda = "Total: %i\nMean: %.2f %%\nMedian: %.2f %%" if input_mode == "percent" else "Total: %i\nMean: %.2f\nMedian: %.2f"
        #legenda = legenda % (len(data), np.mean(data), np.median(data)) if stats_as_legend else legend
        
        self.percent_histogram(data, output_prefix=output_prefix, n_bins=n_bins, title=title, xlabel=xlabel,
                               ylabel=ylabel, extensions=extensions, legend=legend, legend_location=legend_location,
                               input_mode=input_mode, label=label, stats_as_legend=stats_as_legend)

    @staticmethod
    def add_line(axes, start, end, color):
        line = Line2D([start[0], end[0]], [start[1], end[1]], color=color)
        return axes.add_line(line)

    @staticmethod
    def int_histogram(data, output_prefix, n_bins=None, title="", xlabel="%", ylabel="Number",
                      extensions=("png", "svg"), legend=None, legend_location="best"):



        figure = plt.figure()
        subplot = plt.subplot(1, 1, 1)

        maximum = np.max(data)
        minimum = np.min(data)
        number_of_bins = n_bins if n_bins else maximum - minimum + 1
        plt.hist(data, bins=number_of_bins)

        plt.xlim(xmin=minimum - 1, xmax=maximum)

        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if legend:
            plt.legend((legend,), loc=legend_location)
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    @staticmethod
    def extended_percent_histogram(data, output_prefix, n_bins=20, title="", xlabel="%", ylabel="Number", label=None,
                                   extensions=("png", "svg"), legend=None, legend_location="best", input_mode="percent",
                                   xmax=None, xmin=None):

        figure = plt.figure(figsize=(16, 16))


        maximum = np.max(data)

        bins = np.linspace(0, 100, 21)
        bins = np.append(bins, maximum + 1)

        histogram_list = [np.histogram(dataset, bins=bins)[0] for dataset in data]
        bins = bins[:-1]
        #print bins
        #print histogram_list

        number_of_histograms = len(histogram_list)

        width = 5.0 / float((number_of_histograms + 2))
        #print width

        color_cycle = cycle(['b', 'r', 'c', 'g', 'indigo', 'y', 'k', 'olive', 'violet', 'darkgrey', 'gold'])

        color_list = list(islice(color_cycle, None, number_of_histograms))

        for subplot_index in (1, 2):
            subplot = plt.subplot(2, 1, subplot_index)
            for i in range(0, number_of_histograms):
                left = bins + (i + 1) * width
                #print left
                #print histogram_list[i]
                subplot.bar(left, histogram_list[i], width, label=list(label)[i], color=color_list[i],)

            plt.xlim(xmin=0, xmax=105)

            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

            if subplot_index == 1:
                subplot.set_yscale('log', basey=10)
            elif subplot_index == 2:
                if legend:
                    plt.legend((legend,), loc=legend_location)
                if label:
                    plt.legend(loc=legend_location)
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    @staticmethod
    def draw_histogram(data_array, output_prefix=None, number_of_bins=None, width_of_bins=None, bins_list=None,
                       max_threshold=None, min_threshold=None, xlabel=None, ylabel=None,
                       title=None, extensions=("png",), ylogbase=None, xlogbase=None, subplot=None, suptitle=None,
                       close_figure=False, save_histovalues_only=False, figsize=(6, 6),
                       header=None):
        if (number_of_bins is not None) and (width_of_bins is not None):
            raise AttributeError("Options -w/--width_of_bins and -b/--number_of_bins mustn't be set simultaneously")

        if max_threshold and min_threshold:
            if max_threshold < min_threshold:
                raise ValueError("Maximum threshold (%s) is lower than minimum threshold(%s)" % (str(max_threshold),
                                                                                                 str(min_threshold)))

        max_lenn = max(data_array) if (data_array.size > 0) else 0
        min_lenn = min(data_array) if (data_array.size > 0)  else 0

        max_len = max_threshold if (max_threshold is not None) else max_lenn
        min_len = min_threshold if (min_threshold is not None) else min_lenn


        #max_len = max_threshold if (max_threshold is not None) and (max_threshold < max_lenn) else max_lenn
        #min_len = min_threshold if (min_threshold is not None) and (min_lenn < min_threshold) else min_lenn

        filtered = []

        if (max_len < max_lenn) and (min_len > min_lenn):
            for entry in data_array:
                if min_len <= entry <= max_len:
                    filtered.append(entry)
        elif max_len < max_lenn:
            for entry in data_array:
                if entry <= max_len:
                    filtered.append(entry)
        elif min_len > min_lenn:
            for entry in data_array:
                if min_len <= entry:
                    filtered.append(entry)
        else:
            filtered = data_array
        if subplot is None:
            #print "aaaaaaaaaa"
            figure = plt.figure(1, figsize=figsize,)
            subplot = figure.add_subplot(1, 1, 1)
        else:
            plt.axes(subplot)

        if bins_list:
            bins = bins_list
        elif number_of_bins:
            bins = number_of_bins
        elif width_of_bins:
            bins = np.arange(min_len, max_len, width_of_bins)
            #print bins
            #bins[0] += 1
            bins = np.append(bins, [max_len])
        else:
            bins = 30

        n, bins, patches = plt.hist(filtered, bins=bins) # , log=False if ylogbase is None else True)

        bin_centers = (bins + ((bins[1] - bins[0])/2))[:-1]
        #print bin_centers
        #print len(n)
        #print len(bin_centers)
        print ("Minimum x in histogram: %f" % float(min_len))
        print ("Maximum x in histogram: %f" % float(max_len))
        print ("Bins:")
        print (bins)

        plt.xlim(xmin=min_len, xmax=max_len)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)
        if suptitle:
            plt.suptitle(suptitle)

        if ylogbase:
            subplot.set_yscale('log', basey=ylogbase)
            plt.ylim(ymin=1)

        if xlogbase:
            subplot.set_xscale('log', basex=xlogbase)
            plt.xlim(xmin=1)

        if output_prefix:
            if not save_histovalues_only:
                for ext in extensions:
                    plt.savefig("%s.%s" % (output_prefix, ext))

            # save histo values
            np.savetxt("%s.histo" % output_prefix, np.column_stack(bins[:-1], n), fmt="%i\t%i")
            np.savetxt("%s.bins" % output_prefix, bins, fmt="%i")
        if close_figure:
            plt.close(figure)

        return n, bins

        #return zip(bin_centers, n)

    def draw_multi_histogram_picture(self, list_of_data_arrays, subplot_tuple, output_prefix=None,
                                     figsize=(10, 10), number_of_bins_list=None, width_of_bins_list=None,
                                     max_threshold_list=None, min_threshold_list=None, xlabel_list=None, ylabel_list=None,
                                     title_list=None, ylogbase_list=None, label_list=None,
                                     extensions=("png",), suptitle=None, share_y_axis=False,
                                     share_x_axis=False):
        figure = plt.figure(1, figsize=figsize)
        if suptitle:
            plt.suptitle(suptitle)
        if len(subplot_tuple) != 2:
            raise ValueError("Subplot tuple should contain exactly two values, not %i!" % len(subplot_tuple))
        if not (isinstance(subplot_tuple[0], int) and isinstance(subplot_tuple[1], int)):
            raise ValueError("Subplot tuple should contain two values, not (%s, %s)!" % (str(type(subplot_tuple[0])),
                                                                                         str(type(subplot_tuple[1]))))

        number_of_subplots = subplot_tuple[0] * subplot_tuple[1]
        number_of_datasets = len(list_of_data_arrays)

        parameters_list = [number_of_bins_list, width_of_bins_list, max_threshold_list, min_threshold_list,
                           xlabel_list, ylabel_list, title_list, ylogbase_list, label_list]

        subplot_list = []
        for dataset_index in range(0, number_of_datasets):
            parameters = [None, None, None, None, None, None, None, None, None]
            for parameter_index in range(0, 9):
                if parameters_list[parameter_index]:
                    if dataset_index < len(parameters_list[parameter_index]):
                        parameters[parameter_index] = parameters_list[parameter_index][dataset_index]
            if dataset_index > 0:
                if share_x_axis and share_y_axis:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0],
                                                           sharey=subplot_list[0]))
                elif share_x_axis:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0]))
                elif share_y_axis:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0],
                                                           sharey=subplot_list[0]))
                else:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1))
            else:
                subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                       subplot_tuple[1],
                                                       dataset_index + 1))

            histo = self.draw_histogram(list_of_data_arrays[dataset_index],  number_of_bins=parameters[0],
                                        width_of_bins=parameters[1], max_threshold=parameters[2],
                                        min_threshold=parameters[3], xlabel=parameters[4], ylabel=parameters[5],
                                        title=parameters[6], extensions=("png",), ylogbase=parameters[7],
                                        subplot=subplot_list[dataset_index],
                                        suptitle=None)
            #print histo
            if output_prefix and (list_of_data_arrays[dataset_index].size > 0):
                output_histo_file = "%s.%s.%shisto" % (output_prefix,
                                                       dataset_index if parameters[8] is None else parameters[9],
                                                       ("log%i." % parameters[7]) if parameters[7] else "")
                print(histo)
                np.savetxt(output_histo_file, np.column_stack((histo[1][:-1], histo[0])), fmt="%f\t%f")

        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

        return figure

    def draw_tetra_histogram_with_two_logscaled(self, list_of_data_arrays,  output_prefix=None,
                                                figsize=(10, 10), number_of_bins_list=None, width_of_bins_list=None,
                                                max_threshold_list=None, min_threshold_list=None, xlabel=None, ylabel=None,
                                                title_list=None, logbase=10, label_list=None,
                                                extensions=("png",), suptitle=None, share_y_axis=False,
                                                share_x_axis=False):
        subplot_tuple = (2, 2)

        none_list = [None, None, None, None]

        list_of_data = [list_of_data_arrays[0], list_of_data_arrays[1],
                        list_of_data_arrays[0], list_of_data_arrays[1]]

        number_of_bins_listtt = [number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0] if len(number_of_bins_list) == 1 else number_of_bins_list[1],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0] if len(number_of_bins_list) == 1 else number_of_bins_list[1]] if number_of_bins_list else none_list

        width_of_bins_listtt = [width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0] if len(width_of_bins_list) == 1 else width_of_bins_list[1],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0] if len(width_of_bins_list) == 1 else width_of_bins_list[1]] if width_of_bins_list else none_list

        max_threshold_listtt = [max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0] if len(max_threshold_list) == 1 else max_threshold_list[1],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0] if len(max_threshold_list) == 1 else max_threshold_list[1]] if max_threshold_list else none_list

        min_threshold_listtt = [min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0] if len(min_threshold_list) == 1 else min_threshold_list[1],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0] if len(min_threshold_list) == 1 else min_threshold_list[1]] if min_threshold_list else none_list

        title_listtt = [title_list[0], title_list[1],
                        "%s(logscaled)" % title_list[0], "%s(logscaled)" % title_list[1]] if title_list else none_list

        ylogbase_list = [None, None, logbase, logbase]

        xlabel_list = [None,
                       None,
                       xlabel if isinstance(xlabel, str) else xlabel[0] if len(xlabel) == 1 else xlabel[0],
                       xlabel if isinstance(xlabel, str) else xlabel[0] if len(xlabel) == 1 else xlabel[1]]

        ylabel_list = [ylabel if isinstance(ylabel, str) else ylabel[0] if len(ylabel) == 1 else ylabel[0],
                       None,
                       ylabel if isinstance(ylabel, str) else ylabel[0] if len(ylabel) == 1 else ylabel[1],
                       None]

        self.draw_multi_histogram_picture(list_of_data, subplot_tuple, output_prefix=output_prefix,
                                          figsize=figsize, number_of_bins_list=number_of_bins_listtt,
                                          width_of_bins_list=width_of_bins_listtt,
                                          max_threshold_list=max_threshold_listtt,
                                          min_threshold_list=min_threshold_listtt, xlabel_list=xlabel_list, ylabel_list=ylabel_list,
                                          title_list=title_listtt, ylogbase_list=ylogbase_list, label_list=None,
                                          extensions=extensions, suptitle=suptitle,
                                          share_y_axis=share_y_axis, share_x_axis=share_x_axis)

    def draw_hexa_histogram_with_three_logscaled(self, list_of_data_arrays,  output_prefix=None,
                                                 figsize=(10, 15), number_of_bins_list=None, width_of_bins_list=None,
                                                 max_threshold_list=None, min_threshold_list=None, xlabel=None, ylabel=None,
                                                 title_list=None, logbase=10, label_list=None,
                                                 extensions=("png",), suptitle=None, share_y_axis=False,
                                                 share_x_axis=False):
        subplot_tuple = (2, 3)

        none_list = [None, None, None, None, None, None]

        list_of_data = [list_of_data_arrays[0], list_of_data_arrays[1], list_of_data_arrays[2],
                        list_of_data_arrays[0], list_of_data_arrays[1], list_of_data_arrays[2]]

        number_of_bins_listtt = [number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0] if len(number_of_bins_list) == 1 else number_of_bins_list[1],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0] if len(number_of_bins_list) == 1 else number_of_bins_list[2],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0] if len(number_of_bins_list) == 1 else number_of_bins_list[1],
                                 number_of_bins_list if isinstance(number_of_bins_list, str) else number_of_bins_list[0] if len(number_of_bins_list) == 1 else number_of_bins_list[2],] if number_of_bins_list else none_list

        width_of_bins_listtt = [width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0] if len(width_of_bins_list) == 1 else width_of_bins_list[1],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0] if len(width_of_bins_list) == 1 else width_of_bins_list[2],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0] if len(width_of_bins_list) == 1 else width_of_bins_list[1],
                                width_of_bins_list if isinstance(width_of_bins_list, str) else width_of_bins_list[0] if len(width_of_bins_list) == 1 else width_of_bins_list[2]] if width_of_bins_list else none_list

        max_threshold_listtt = [max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0] if len(max_threshold_list) == 1 else max_threshold_list[1],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0] if len(max_threshold_list) == 1 else max_threshold_list[2],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0] if len(max_threshold_list) == 1 else max_threshold_list[1],
                                max_threshold_list if isinstance(max_threshold_list, str) else max_threshold_list[0] if len(max_threshold_list) == 1 else max_threshold_list[2]] if max_threshold_list else none_list

        min_threshold_listtt = [min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0] if len(min_threshold_list) == 1 else min_threshold_list[1],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0] if len(min_threshold_list) == 1 else min_threshold_list[2],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0] if len(min_threshold_list) == 1 else min_threshold_list[1],
                                min_threshold_list if isinstance(min_threshold_list, str) else min_threshold_list[0] if len(min_threshold_list) == 1 else min_threshold_list[2],] if min_threshold_list else none_list

        """


        number_of_bins_listtt = [number_of_bins_list[0], number_of_bins_list[1], number_of_bins_list[2],
                                 number_of_bins_list[0], number_of_bins_list[1], number_of_bins_list[2]] if number_of_bins_list else none_list
        width_of_bins_listtt = [width_of_bins_list[0], width_of_bins_list[1], width_of_bins_list[2],
                                width_of_bins_list[0], width_of_bins_list[1], width_of_bins_list[2]] if width_of_bins_list else none_list
        max_threshold_listtt = [max_threshold_list[0], max_threshold_list[1], max_threshold_list[2],
                                max_threshold_list[0], max_threshold_list[1], max_threshold_list[2]] if max_threshold_list else none_list

        min_threshold_listtt = [min_threshold_list[0], min_threshold_list[1], min_threshold_list[2],
                                min_threshold_list[0], min_threshold_list[1], min_threshold_list[2]] if min_threshold_list else none_list

        """

        title_listtt = [title_list[0], title_list[1], title_list[2],
                        "%s, logscaled" % title_list[0], "%s, logscaled" % title_list[1], "%s, logscaled" % title_list[2]] if title_list else none_list

        ylogbase_list = [None, None, None, logbase, logbase, logbase]

        xlabel_list = [None,
                       None,
                       None,
                       xlabel if isinstance(xlabel, str) else xlabel[0] if len(xlabel) == 1 else xlabel[0],
                       xlabel if isinstance(xlabel, str) else xlabel[0] if len(xlabel) == 1 else xlabel[1],
                       xlabel if isinstance(xlabel, str) else xlabel[0] if len(xlabel) == 1 else xlabel[2]]

        ylabel_list = [ylabel if isinstance(ylabel, str) else ylabel[0] if len(ylabel) == 1 else ylabel[0],
                       None,
                       None,
                       ylabel if isinstance(ylabel, str) else ylabel[0] if len(ylabel) == 1 else ylabel[1],
                       None,
                       None]

        """
        xlabel_list = [None,
                       None,
                       None,
                       xlabel if isinstance(xlabel, str) else xlabel[0],
                       xlabel if isinstance(xlabel, str) else xlabel[1],
                       xlabel if isinstance(xlabel, str) else xlabel[2]]
        ylabel_list = [ylabel if isinstance(ylabel, str) else ylabel[0],
                       None,
                       None,
                       ylabel if isinstance(ylabel, str) else ylabel[1],
                       None,
                       None]
        """

        self.draw_multi_histogram_picture(list_of_data, subplot_tuple, output_prefix=output_prefix,
                                          figsize=figsize, number_of_bins_list=number_of_bins_listtt,
                                          width_of_bins_list=width_of_bins_listtt,
                                          max_threshold_list=max_threshold_listtt,
                                          min_threshold_list=min_threshold_listtt, xlabel_list=xlabel_list, ylabel_list=ylabel_list,
                                          title_list=title_listtt, ylogbase_list=ylogbase_list, label_list=None,
                                          extensions=extensions, suptitle=suptitle, share_y_axis=share_y_axis,
                                          share_x_axis=share_x_axis)

    def draw_tetra_histogram_with_two_logscaled_from_file(self, file_list,  column_idx_list, output_prefix,
                                                          figsize=(10, 10), number_of_bins_list=None, width_of_bins_list=None,
                                                          max_threshold_list=None, min_threshold_list=None, xlabel=None, ylabel=None,
                                                          title_list=None, logbase=10, label_list=None,
                                                          extensions=("png",), suptitle=None, separator=None,
                                                          share_y_axis=False, share_x_axis=False,
                                                          comments="#"):

        list_of_data_arrays = []
        for filename, column_idx in zip(file_list, column_idx_list):
            #print filename, column_idx
            list_of_data_arrays.append(np.loadtxt(filename, usecols=(column_idx,),
                                                  delimiter=separator, comments=comments))

        self.draw_tetra_histogram_with_two_logscaled(list_of_data_arrays, output_prefix=output_prefix,
                                                     figsize=figsize, number_of_bins_list=number_of_bins_list,
                                                     width_of_bins_list=width_of_bins_list,
                                                     max_threshold_list=max_threshold_list,
                                                     min_threshold_list=min_threshold_list, xlabel=xlabel,
                                                     ylabel=ylabel, title_list=title_list, logbase=logbase,
                                                     label_list=label_list, extensions=extensions, suptitle=suptitle,
                                                     share_y_axis=share_y_axis, share_x_axis=share_x_axis)

    def draw_hexa_histogram_with_three_logscaled_from_file(self, list_of_files, column_idx_list, output_prefix,
                                                           figsize=(15, 10), number_of_bins_list=None,
                                                           width_of_bins_list=None,
                                                           max_threshold_list=None, min_threshold_list=None,
                                                           xlabel=None, ylabel=None,
                                                           title_list=None, logbase=10, label_list=None,
                                                           extensions=("png",), suptitle=None, separator=None,
                                                           share_y_axis=False, share_x_axis=False,
                                                           comments="#"):

        list_of_data_arrays = []
        for filename, column_idx in zip(list_of_files, column_idx_list if column_idx_list is not None else [0] * len(list_of_files)):
            list_of_data_arrays.append(np.loadtxt(filename, usecols=(column_idx,),
                                                  delimiter=separator, comments=comments))

        self.draw_hexa_histogram_with_three_logscaled(list_of_data_arrays,  output_prefix=output_prefix,
                                                      figsize=figsize, number_of_bins_list=number_of_bins_list,
                                                      width_of_bins_list=width_of_bins_list,
                                                      max_threshold_list=max_threshold_list,
                                                      min_threshold_list=min_threshold_list, xlabel=xlabel,
                                                      ylabel=ylabel, title_list=title_list, logbase=logbase,
                                                      label_list=label_list, extensions=extensions, suptitle=suptitle,
                                                      share_y_axis=share_y_axis, share_x_axis=share_x_axis)

    @staticmethod
    def draw_histogram_from_file(input_file, output_prefix, number_of_bins=None, width_of_bins=None,
                                 separator="\n", max_length=None, min_length=1, xlabel=None, ylabel=None,
                                 title=None, extensions=("png",), logbase=10):
        if (number_of_bins is not None) and (width_of_bins is not None):
            raise AttributeError("Options -w/--width_of_bins and -b/--number_of_bins mustn't be set simultaneously")

        if min_length < 0:
            raise ValueError("Minimum length can't be negative")
        if (max_length is not None) and (max_length < 0):
            raise ValueError("Maximum length can't be negative")

        lengths = np.fromfile(input_file, sep=separator)

        max_lenn = max(lengths)
        min_lenn = min(lengths)

        #max_len = max_length if (max_length is not None) and (max_length < max_lenn) else max_lenn
        #min_len = min_length if (min_length is not None) and min_lenn < min_length else min_lenn

        max_len = max_length if (max_length is not None) else max_lenn
        min_len = min_length if (min_length is not None) else min_lenn

        #print(lengths)
        #print(max_len)
        #print(min_len)

        filtered = []

        if (max_len < max_lenn) and (min_len > min_lenn):
            for entry in lengths:
                if min_len <= entry <= max_len:
                    filtered.append(entry)
        elif max_len < max_lenn:
            for entry in lengths:
                if entry <= max_len:
                    filtered.append(entry)
        elif min_len > min_lenn:
            for entry in lengths:
                if min_len <= entry:
                    filtered.append(entry)
        else:
            filtered = lengths

        figure = plt.figure(1, figsize=(6, 6))
        subplot = plt.subplot(1, 1, 1)

        if number_of_bins:
            bins = number_of_bins
        elif width_of_bins:
            bins = np.arange(min_len, max_len, width_of_bins)
            #print bins
            #bins[0] += 1
            bins = np.append(bins, [max_len])
        else:
            bins = 30

        n, bins, patches = plt.hist(filtered, bins=bins)

        print("Minimum x in histogram: %f" % float(min_len))
        print("Maximum x in histogram: %f" % float(max_len))
        print("Bins:")
        print(bins)

        #bin_centers = (bins + ((bins[1] - bins[0])/2))[:-1]
        #print bin_centers
        #print len(n)
        #print len(bin_centers)

        plt.xlim(xmin=min_len, xmax=max_len)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

        subplot.set_yscale('log', base=logbase)
        #subplot.set_xscale('log', basex=args.logbase)
        for ext in extensions:
            plt.savefig("%s.logscale.%s" % (output_prefix, ext))

        # save histo values
        np.savetxt("%s.histo" % output_prefix, zip(bins[:-1], n), fmt="%i\t%i")
        np.savetxt("%s.bins" % output_prefix, bins)

    @staticmethod
    def generate_bin_array_by_width(min_value, max_value, bin_width, add_max_value=True):

        #print min_value, max_value, bin_width
        bin_array = np.arange(min_value, max_value, bin_width)
        if add_max_value:
            bin_array = np.append(bin_array, [max_value])

        return bin_array
    """
    def generate_bin_array(self, x, y, bin_number=20, bin_width=None, bin_array=None,
                           min_x_value=None, max_x_value=None, min_y_value=None, max_y_value=None, add_max_value=True):
        if (bin_width is not None) and (bin_array is not None):
            raise ValueError("Both bin width and bin array were set")

        min_x, max_x = min(x), max(x)
        min_y, max_y = min(y), max(y)

        if bin_width:
            xbins = self.generate_bin_array_by_width(min_x_value if min_x_value is not None else min_x,
                                                     max_x_value if max_x_value is not None else max_x,
                                                     bin_width if isinstance(bin_width, int) else bin_width[0],
                                                     add_max_value=add_max_value)
            ybins = self.generate_bin_array_by_width(min_y_value if min_y_value is not None else min_y,
                                                     max_y_value if max_y_value is not None else max_y,
                                                     bin_width if isinstance(bin_width, int) else bin_width[1],
                                                     add_max_value=add_max_value)
            bins = (xbins, ybins)

        elif bin_array:
            bins = bin_array
        else:
            xbins = np.linspace(min_x_value if min_x_value is not None else min_x,
                                max_x_value if max_x_value is not None else max_x,
                                bin_number if isinstance(bin_number, int) else bin_number[0])
            ybins = np.linspace(min_y_value if min_y_value is not None else min_y,
                                max_y_value if max_y_value is not None else max_y,
                                bin_number if isinstance(bin_number, int) else bin_number[1])
            bins = (xbins, ybins)

        return bins
    """

    def generate_bin_array(self, x_list, y_list=None, bin_number=20, bin_width=None, bin_array=None,
                           min_x_value=None, max_x_value=None, min_y_value=None,
                           max_y_value=None, add_max_value=True):
        if (bin_width is not None) and (bin_array is not None):
            raise ValueError("Both bin width and bin array were set")
        #print(x_list)
        min_x, max_x = min(map(min, x_list) if isinstance(x_list[0], Iterable) else x_list), \
                       max(map(max, x_list) if isinstance(x_list[0], Iterable) else x_list)
        if not(y_list is None):
            min_y, max_y = min(map(min, y_list) if isinstance(y_list[0], Iterable) else y_list), \
                           max(map(max, y_list) if isinstance(y_list[0], Iterable) else y_list)

        if bin_width:
            xbins = self.generate_bin_array_by_width(min_x_value if min_x_value is not None else min_x,
                                                     max_x_value if max_x_value is not None else max_x,
                                                     bin_width if isinstance(bin_width, int) else bin_width[0],
                                                     add_max_value=add_max_value)
            if not(y_list is None):
                ybins = self.generate_bin_array_by_width(min_y_value if min_y_value is not None else min_y,
                                                         max_y_value if max_y_value is not None else max_y,
                                                         bin_width if isinstance(bin_width, int) else bin_width[1],
                                                         add_max_value=add_max_value)
            bins = (xbins, ybins) if not(y_list is None) else xbins

        elif bin_array:
            bins = bin_array
        else:
            if bin_number is None:
                print("WARNINNG!!! No bin_number or bin_width or bin_array were set. "
                      "Instead default value(20) for bin number is used.")
            xbins = np.linspace(min_x_value if min_x_value is not None else min_x,
                                max_x_value if max_x_value is not None else max_x,
                                21 if bin_number is None else bin_number + 1 if isinstance(bin_number, int) else bin_number[0] + 1)
            if not(y_list is None):
                ybins = np.linspace(min_y_value if min_y_value is not None else min_y,
                                    max_y_value if max_y_value is not None else max_y,
                                    21 if bin_number is None else bin_number + 1 if isinstance(bin_number, int) else bin_number[1] + 1)
            bins = (xbins, ybins) if not(y_list is None) else xbins

        return bins

    def draw_heatmap(self, x, y, output_prefix=None, xlabel=None, ylabel=None, title=None,
                     figsize=(8, 8), minimum_counts_to_show=1,
                     extensions=("png", "svg"), show_colorbar=True, bin_number=20, bin_width=None, bin_array=None,
                     min_x_value=None, max_x_value=None, min_y_value=None, max_y_value=None, add_max_value=True,
                     subplot=None, save_histovalues_only=False, header=None, logscaled=False):
        """
        bin_width: numeric or tuple of two numerics
        """

        bins = self.generate_bin_array(x, y, bin_number=bin_number, bin_width=bin_width, bin_array=bin_array,
                                       min_x_value=min_x_value, max_x_value=max_x_value,
                                       min_y_value=min_y_value, max_y_value=max_y_value,
                                       add_max_value=add_max_value)
        print("X axis bins:")
        print(bins[0])
        print("Y axis bins:")
        print(bins[1])

        if subplot is None:
            #print "aaaaaaaaaa"
            #figure = plt.figure(1, figsize=figsize,)
            #subplott = figure.add_subplot(1, 1, 1)
            fig, ax = plt.subplots(figsize=figsize)
        else:
            ax = plt.axes(subplot)

        counts, xedges, yedges, image = plt.hist2d(x, y, bins, cmin=minimum_counts_to_show,
                                                   norm=LogNorm() if logscaled else None)
        print("X edges bins:")
        print(xedges)
        print("Y edges bins:")
        print(yedges)
        print(counts)
        #print x
        #print y
        #print minimum_counts_to_show
        plt.xlim(xmin=min_x_value, xmax=max_x_value)
        plt.ylim(ymin=min_y_value, ymax=max_y_value)

        if xlabel:
            plt.xlabel(xlabel)

        if ylabel:
            plt.ylabel(ylabel)

        if title:
            plt.title(title)

        #np.savetxt("tara.t", counts, delimiter='\n')
        if show_colorbar:
            max_counts = np.nanmax(counts)
            #print max_counts
            plt.colorbar(image, ax=ax)
            """
            cmap = plt.get_cmap('jet', max_counts)
            #cmap.set_under('gray')
            mappable = plt.cm.ScalarMappable(cmap=cmap)
            mappable.set_array([])
            mappable.set_clim(1, max_counts + 1)
            #mappable.set_array([])
            #mappable.set_clim(-0.5, ncolors+0.5)
            colorbar = plt.colorbar(mappable)
            #colorbar.set_ticks(np.linspace(1 + 0.5, max_counts + 0.5, max_counts), minor=True)
            decimal = int(max_counts / 10) + 1
            max_major_tick = max_counts - (max_counts % decimal)
            major_ticks = np.linspace(decimal + 0.5, max_major_tick + 0.5, int(max_major_tick / decimal))
            colorbar.set_ticks(major_ticks)
            print (decimal, max_major_tick + decimal, decimal)
            colorbar.set_ticklabels(range(decimal, int(max_major_tick + decimal), decimal))
            """
        if output_prefix:
            if not save_histovalues_only:
                for ext in extensions:
                    plt.savefig("%s.%s" % (output_prefix, ext))
            with open("%s.2d_histo" % output_prefix, "w") as out_fd:
                if header:
                    out_fd.write(header if header[-1] == "\n" else header + "\n")
                for ix in range(0, len(xedges) - 1):
                    for iy in range(0, len(yedges) - 1):
                        out_fd.write("%d\t%d\t%d\n" % (xedges[ix], yedges[iy], 0 if np.isnan(counts[ix, iy]) else counts[ix, iy]))

        return counts, xedges, yedges

    def draw_percent_heatmap(self, x, y, output_prefix=None, xlabel=None, ylabel=None, title=None,
                             figsize=(8, 8), minimum_counts_to_show=1,
                             extensions=("png", "svg"), show_colorbar=True, bin_number=20, bin_width=None, bin_array=None,
                             type="percent", add_max_value=True,
                             subplot=None, save_histovalues_only=False, header=None, logscaled=False):
        if type == "percent":
            return self.draw_heatmap(x, y, output_prefix, xlabel=xlabel, ylabel=ylabel, title=title,
                                     figsize=figsize, minimum_counts_to_show=minimum_counts_to_show,
                                     extensions=extensions, show_colorbar=show_colorbar, bin_number=bin_number,
                                     bin_width=bin_width, bin_array=bin_array, min_x_value=0,
                                     max_x_value=100, min_y_value=0, max_y_value=100,
                                     add_max_value=add_max_value, subplot=subplot,
                                     save_histovalues_only=save_histovalues_only, header=header,
                                     logscaled=logscaled)
        elif type == "fraction":
            return self.draw_heatmap(x, y, output_prefix, xlabel=xlabel, ylabel=ylabel, title=title,
                                     figsize=figsize, minimum_counts_to_show=minimum_counts_to_show,
                                     extensions=extensions, show_colorbar=show_colorbar, bin_number=bin_number,
                                     bin_width=bin_width, bin_array=bin_array, min_x_value=0.0,
                                     max_x_value=1.0, min_y_value=0.0, max_y_value=1.0,
                                     add_max_value=add_max_value, subplot=subplot,
                                     save_histovalues_only=save_histovalues_only, header=header,
                                     logscaled=logscaled)
        else:
            raise ValueError("Unrecognized type of percent histogram(neither fraction nor percent)")

    def draw_heatmap_from_file(self, input_file, output_prefix, x_column=0, y_column=1, xlabel=None, ylabel=None,
                               title=None, figsize=(8, 8), minimum_counts_to_show=1, extensions=("png", "svg"),
                               show_colorbar=True, bin_number=20, bin_width=None, bin_array=None,
                               min_x_value=None, max_x_value=None, min_y_value=None, max_y_value=None,
                               add_max_value=True):

        x = np.loadtxt(input_file, usecols=(x_column,))
        y = np.loadtxt(input_file, usecols=(y_column,))

        self.draw_heatmap(x, y, output_prefix=output_prefix, xlabel=xlabel, ylabel=ylabel, title=title,
                          figsize=figsize, minimum_counts_to_show=minimum_counts_to_show, extensions=extensions,
                          show_colorbar=show_colorbar, bin_number=bin_number, bin_width=bin_width, bin_array=bin_array,
                          min_x_value=min_x_value, max_x_value=max_x_value,
                          min_y_value=min_y_value, max_y_value=max_y_value,
                          add_max_value=add_max_value,)

    def draw_double_histo_from_file(self, file_list, column_idx_list, subplot_tuple=(2, 1), output_prefix=None,
                                    figsize=(5, 10), number_of_bins_list=None, width_of_bins_list=None,
                                    max_threshold_list=None, min_threshold_list=None, xlabel_list=None, ylabel_list=None,
                                    title_list=None, ylogbase_list=None, label_list=None,
                                    extensions=("png",), suptitle=None, separator=None, comments='#',
                                    share_y_axis=False,
                                    share_x_axis=False):
        list_of_data_arrays = []
        for filename, column_idx in zip(file_list, column_idx_list):
            #print filename, column_idx
            list_of_data_arrays.append(np.loadtxt(filename, usecols=(column_idx,),
                                                  delimiter=separator, comments=comments))

        self.draw_multi_histogram_picture(list_of_data_arrays, subplot_tuple, output_prefix=output_prefix,
                                          figsize=figsize, number_of_bins_list=number_of_bins_list,
                                          width_of_bins_list=width_of_bins_list,
                                          max_threshold_list=max_threshold_list, min_threshold_list=min_threshold_list,
                                          xlabel_list=xlabel_list, ylabel_list=ylabel_list,
                                          title_list=title_list, ylogbase_list=ylogbase_list, label_list=label_list,
                                          extensions=extensions, suptitle=suptitle,
                                          share_y_axis=share_y_axis, share_x_axis=share_x_axis)

    @staticmethod
    def venn_diagram_from_sets(set1, set2, set3=None, set_labels=None, set_colors=None,
                               output_prefix=None, extensions=("png",), title=None):
        if set3 is None:
            diagramm = venn(subsets=(set1, set2), set_labels=set_labels, set_colors=set_colors)
        else:
            diagramm = venn(subsets=(set1, set2, set3), set_labels=set_labels, set_colors=set_colors)

        if title:
            plt.title(title)

        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

        return diagramm

    def venn_diagram_from_sets_from_files(self, set1_file, set2_file, set3_file=None, set_labels=None, set_colors=None,
                                          output_prefix=None, extensions=("png",), title=None):
        set1 = IdSet(filename=set1_file)
        set2 = IdSet(filename=set2_file)
        set3 = IdSet(filename=set3_file) if set3_file else None

        return self.venn_diagram_from_sets(set1, set2, set3=set3, set_labels=set_labels, set_colors=set_colors,
                                           output_prefix=output_prefix, extensions=extensions, title=title)

    def draw_bar_plot(self, input_data, output_prefix, extentions=["png"],
                      xlabel=None, ylabel=None, title=None, min_value=None, max_value=None, new_figure=True,
                      figsize=(6, 6), close_figure=True):
        
        data = np.array(input_data)
        if min_value and max_value:
            data = data[(data <= min_value) & (data >= max_value)]
        elif max_value:
            data = data[data <= max_value]
        elif min_value:
            data = data[data >= min_value]
        if new_figure:
            plt.figure(1, figsize=figsize)
            plt.subplot(1, 1, 1)
        plt.bar(np.arange(1, len(data) + 1, 1), data)
        plt.xlim(xmin=0, xmax=len(data))
        
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)
        
        for ext in extentions:
            plt.savefig("%s.%s" % (output_prefix, ext))

        plt.close()

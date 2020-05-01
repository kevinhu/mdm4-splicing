import json
import huygens as huy
import galileo as gal
import os
import sys
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm
from textwrap import wrap


prop = fm.FontProperties(fname='../plots/arial.ttf')

plt.rcParams['ps.useafm'] = True
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


sys.path.append(os.path.relpath("../../huygens"))
sys.path.append(os.path.relpath("../../galileo"))


with open("experiments.json", "r") as f:
    exp = json.load(f)

    experiments = exp["experiments"]
    experiment_ids = exp["experiment_ids"]
    display_names = exp["display_names"]
    display_groups = exp["display_groups"]
    contexts = exp["contexts"]

kallisto_sleuth_path = "../data/processed/kallisto_sleuth_merge/"

rpl22_oe_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22_oe_genes.h5", key="sleuth_diff")
rpl22l1_oe_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22l1_oe_genes.h5", key="sleuth_diff")
rpl22l1_kd1_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22l1_kd1_genes.h5", key="sleuth_diff")
rpl22l1_kd2_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22l1_kd2_genes.h5", key="sleuth_diff")
rpl22_a_ko1_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22_a_ko1_genes.h5", key="sleuth_diff")
rpl22_a_ko2_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22_a_ko2_genes.h5", key="sleuth_diff")
rpl22_b_ko1_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22_b_ko1_genes.h5", key="sleuth_diff")
rpl22_b_ko2_genes = pd.read_hdf(kallisto_sleuth_path + "rpl22_b_ko2_genes.h5", key="sleuth_diff")

rpl22_oe_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22_oe_transcripts.h5", key="sleuth_diff")
rpl22l1_oe_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22l1_oe_transcripts.h5", key="sleuth_diff")
rpl22l1_kd1_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22l1_kd1_transcripts.h5", key="sleuth_diff")
rpl22l1_kd2_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22l1_kd2_transcripts.h5", key="sleuth_diff")
rpl22_a_ko1_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22_a_ko1_transcripts.h5", key="sleuth_diff")
rpl22_a_ko2_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22_a_ko2_transcripts.h5", key="sleuth_diff")
rpl22_b_ko1_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22_b_ko1_transcripts.h5", key="sleuth_diff")
rpl22_b_ko2_transcripts = pd.read_hdf(kallisto_sleuth_path + "rpl22_b_ko2_transcripts.h5",key="sleuth_diff")


rpl22_oe_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22_oe.txt",sep="\t",index_col=0)
rpl22l1_oe_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22l1_oe.txt",sep="\t",index_col=0)
rpl22l1_kd1_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22l1_kd1.txt",sep="\t",index_col=0)
rpl22l1_kd2_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22l1_kd2.txt",sep="\t",index_col=0)
rpl22_a_ko1_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22_a_ko1.txt",sep="\t",index_col=0)
rpl22_a_ko2_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22_a_ko2.txt",sep="\t",index_col=0)
rpl22_b_ko1_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22_b_ko1.txt",sep="\t",index_col=0)
rpl22_b_ko2_rmats = pd.read_csv("../data/processed/rmats_merge/rpl22_b_ko2.txt",sep="\t",index_col=0)

rpl22_oe_rmats = rpl22_oe_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22l1_oe_rmats = rpl22l1_oe_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22l1_kd1_rmats = rpl22l1_kd1_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22l1_kd2_rmats = rpl22l1_kd2_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22_a_ko1_rmats = rpl22_a_ko1_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22_a_ko2_rmats = rpl22_a_ko2_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22_b_ko1_rmats = rpl22_b_ko1_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)
rpl22_b_ko2_rmats = rpl22_b_ko2_rmats.rename({"PValue":"pval","FDR":"qval"},axis=1)

splice_types = ["A3SS","A5SS","MXE","RI","SE"]


def as_si(x, ndp):
    """
    Convert a number to scientific notation

    Parameters
    ----------
    x : float
        number to convert
    ndp: float
        number of decimal places

    Returns
    -------
    x_si : string
        x formatted in scientific notation
    """

    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    x_si = r'{m:s} × $10^{{{e:d}}}$'.format(m=m, e=int(e))

    return x_si


def three_bars(annotation_id, 
               experiment_id_1, 
               experiment_id_2, 
               diff_results_1, 
               diff_results_2, ax=None, xlabel=None, ylabel=None):

    if ax is None:
        ax = plt.subplot(111)

    if xlabel is None:
        xlabel = experiment_id_1

    if ylabel is None:
        ylabel = "mRNA expression, log2(TPM + 1)"

    # get experiment controls and treatments
    select_abundance_1 = diff_results_1.loc[annotation_id]
    select_abundance_2 = diff_results_2.loc[annotation_id]

    controls = experiments[experiment_id_1][0]
    treatments_1 = experiments[experiment_id_1][1]
    treatments_2 = experiments[experiment_id_2][1]

    # extract control and treatment values
    control_values = np.log2(select_abundance_1[controls].astype(np.float64)+1)
    treatment_1_values = np.log2(
        select_abundance_1[treatments_1].astype(np.float64)+1)
    treatment_2_values = np.log2(
        select_abundance_2[treatments_2].astype(np.float64)+1)

    control_mean = np.mean(control_values)
    treatment_1_mean = np.mean(treatment_1_values)
    treatment_2_mean = np.mean(treatment_2_values)

    # draw bars based on means
    offset = 0.05
    width = 0.4

    control_rect = Rectangle([offset, 0],
                             width,
                             control_mean,
                             color=control_color,
                             alpha=1,
                             linewidth=0,
                             zorder=-100
                             )
    ax.add_patch(control_rect)

    treatment_1_rect = Rectangle([0.5+offset, 0],
                                 width,
                                 treatment_1_mean,
                                 color=treatment_color,
                                 alpha=1,
                                 linewidth=0,
                                 zorder=-100
                                 )
    ax.add_patch(treatment_1_rect)

    treatment_2_rect = Rectangle([1+offset, 0],
                                 width,
                                 treatment_2_mean,
                                 color=treatment_color,
                                 alpha=1,
                                 linewidth=0,
                                 zorder=-100
                                 )
    ax.add_patch(treatment_2_rect)

    # draw the points themselves
    ax.scatter([0.25-width/4, 0.25, 0.25+width/4],
               control_values,
               color="white",
               linewidth=1,
               s=32,
               edgecolor=control_color
               )
    ax.scatter([0.75-width/4, 0.75, 0.75+width/4],
               treatment_1_values,
               color="white",
               linewidth=1,
               s=32,
               edgecolor=treatment_color
               )
    ax.scatter([1.25-width/4, 1.25, 1.25+width/4],
               treatment_2_values,
               color="white",
               linewidth=1,
               s=32,
               edgecolor=treatment_color
               )

    # figure formatting
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([0.25, 0.75, 1.25])
    ax.set_xticklabels(["", ""])
    ax.set_xlabel(xlabel, rotation=45)
    ax.set_ylabel(ylabel)

    # space out axes
    ax.spines['bottom'].set_position(('axes', 0))
    ax.spines['left'].set_position(('axes', -0.25))

    # set y minimum to 0
    plt.ylim(-0.01)

    y_max = max(list(control_values)+list(treatment_1_values) +
                list(treatment_2_values))

    qval_1 = select_abundance_1["qval"]
    qval_2 = select_abundance_2["qval"]

    treatment_1_max = max(treatment_1_values)
    treatment_2_max = max(treatment_2_values)

    if not np.isnan(qval_1):

        if qval_1 < 0.001:
            compare_text = "**"
        elif qval_1 < 0.01:
            compare_text = "*"
        else:
            compare_text = "n.s"

        ax.text(0.75,
                treatment_1_max*1.15,
                compare_text,
                ha="center",
                fontsize=12,
                )

    if not np.isnan(qval_2):

        if qval_2 < 0.001:
            compare_text = "**"
        elif qval_2 < 0.01:
            compare_text = "*"
        else:
            compare_text = "n.s"

        ax.text(1.25,
                treatment_2_max*1.15,
                compare_text,
                ha="center",
                fontsize=12,
                )

    return ax, y_max


control_color = "#393e46"
treatment_color = "#f67280"
alpha = 0.5


def bars(annotation_id, experiment_id, diff_results, ax=None, xlabel=None, ylabel=None):

    if ax is None:
        ax = plt.subplot(111)

    if xlabel is None:
        xlabel = experiment_id

    if ylabel is None:
        ylabel = "mRNA expression, log2(TPM + 1)"

    # get experiment controls and treatments
    select_abundance = diff_results.loc[annotation_id]

    controls = experiments[experiment_id][0]
    treatments = experiments[experiment_id][1]

    # extract control and treatment values
    control_values = np.log2(select_abundance[controls].astype(np.float64)+1)
    treatment_values = np.log2(
        select_abundance[treatments].astype(np.float64)+1)

    control_mean = np.mean(control_values)
    treatment_mean = np.mean(treatment_values)

    # draw bars based on means
    offset = 0.05
    width = 0.4

    control_rect = Rectangle([offset, 0],
                             width,
                             control_mean,
                             color=control_color,
                             alpha=1,
                             linewidth=0,
                             zorder=-100
                             )
    ax.add_patch(control_rect)

    treatment_rect = Rectangle([0.5+offset, 0],
                               width,
                               treatment_mean,
                               color=treatment_color,
                               alpha=1,
                               linewidth=0,
                               zorder=-100
                               )
    ax.add_patch(treatment_rect)

    # draw the points themselves
    ax.scatter([0.25-width/4, 0.25, 0.25+width/4],
               control_values,
               color="white",
               linewidth=1,
               s=32,
               edgecolor=control_color
               )
    ax.scatter([0.75-width/4, 0.75, 0.75+width/4],
               treatment_values,
               color="white",
               linewidth=1,
               s=32,
               edgecolor=treatment_color
               )

    # figure formatting
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([0.25, 0.75])
    ax.set_xticklabels(["", ""])
    ax.set_xlabel(xlabel, rotation=45)
    ax.set_ylabel(ylabel)

    # space out axes
    ax.spines['bottom'].set_position(('axes', 0))
    ax.spines['left'].set_position(('axes', -0.25))

    # set y minimum to 0
    plt.ylim(-0.01)

    y_max = max(list(control_values)+list(treatment_values))
    
    treatment_max = max(treatment_values)

    qval = select_abundance["qval"]

    if not np.isnan(qval):

        if qval < 0.001:
            compare_text = "**"
        elif qval < 0.01:
            compare_text = "*"
        else:
            compare_text = "n.s"

        ax.text(0.75,
                treatment_max*1.15,
                compare_text,
                ha="center",
                fontsize=12
                )

    return ax, y_max


def all_bars(annotation_id, annotation_type, legend=False):

    plt.figure(figsize=(4, 2.5))

    axes_widths = [2, 2, 3, 3, 3]
    total_width = sum(axes_widths)

    cumulative_widths = [sum(axes_widths[:x]) for x in range(len(axes_widths))]
    
    plot_axes_height = 4
    label_axes_height = 2

    axes = [plt.subplot2grid((plot_axes_height+label_axes_height, total_width), (0, cumulative_widths[x]),
                             colspan=axes_widths[x], rowspan=plot_axes_height) for x in range(len(axes_widths))]
    
    exp_label_axes = [plt.subplot2grid((plot_axes_height+label_axes_height, total_width), (plot_axes_height, cumulative_widths[x]),
                             colspan=axes_widths[x], rowspan=label_axes_height) for x in range(len(axes_widths))]

    maxes = []

    if annotation_type == "transcript":
        sleuth_sets = [rpl22_oe_transcripts,
                       rpl22l1_oe_transcripts,
                       [rpl22l1_kd1_transcripts,
                        rpl22l1_kd2_transcripts],
                       [rpl22_a_ko1_transcripts,
                        rpl22_a_ko2_transcripts],
                       [rpl22_b_ko1_transcripts,
                        rpl22_b_ko2_transcripts]
                       ]

    elif annotation_type == "gene":
        sleuth_sets = [rpl22_oe_genes,
                       rpl22l1_oe_genes,
                       [rpl22l1_kd1_genes,
                        rpl22l1_kd2_genes],
                       [rpl22_a_ko1_genes,
                        rpl22_a_ko2_genes],
                       [rpl22_b_ko1_genes,
                        rpl22_b_ko2_genes]
                       ]

    elif annotation_type == "splicing":
        sleuth_sets = [rpl22_oe_rmats,
                       rpl22l1_oe_rmats,
                       [rpl22l1_kd1_rmats,
                        rpl22l1_kd2_rmats],
                       [rpl22_a_ko1_rmats,
                        rpl22_a_ko2_rmats],
                       [rpl22_b_ko1_rmats,
                        rpl22_b_ko2_rmats]
                       ]

    for sleuth_idx, sleuth_set in enumerate(sleuth_sets[:2]):
        ax = axes[sleuth_idx]

        ax, y_max = bars(annotation_id,
                         experiment_ids[sleuth_idx],
                         sleuth_set,
                         ax=ax,
                         ylabel="")

        maxes.append(y_max)

        if sleuth_idx > 0:

            ax.spines["left"].set_visible(False)
            ax.tick_params(axis='y', which='both', right=False,
                           left=False, labelleft=False)

        ax.set_ylim(0)

        ax.set_xticklabels(display_groups[sleuth_idx], rotation=45, ha="right", rotation_mode="anchor", size=8)
        ax.set_xlabel("")
        
        exp_ax = exp_label_axes[sleuth_idx]
        exp_ax.patch.set_alpha(0)
        exp_ax.spines["left"].set_visible(False)
        exp_ax.spines["top"].set_visible(False)
        exp_ax.spines["right"].set_visible(False)
        exp_ax.set_xticks([])
        exp_ax.set_yticks([])
        exp_ax.set_xlabel(contexts[sleuth_idx])

    for sleuth_idx, sleuth_set in enumerate(sleuth_sets[2:]):
        ax = axes[2+sleuth_idx]

        ax, y_max = three_bars(annotation_id,
                               experiment_ids[2+sleuth_idx*2],
                               experiment_ids[2+sleuth_idx*2+1],
                               sleuth_set[0],
                               sleuth_set[1],
                               ax=ax,
                               ylabel="")

        maxes.append(y_max)

        ax.spines["left"].set_visible(False)
        ax.tick_params(axis='y', which='both', right=False,
                       left=False, labelleft=False)

        ax.set_ylim(0)

        ax.set_xticklabels(
            display_groups[2+sleuth_idx], rotation=45, ha="right", rotation_mode="anchor", size=8)
        ax.set_xlabel("")
        
        exp_ax = exp_label_axes[2+sleuth_idx]
        exp_ax.patch.set_alpha(0)
        exp_ax.spines["left"].set_visible(False)
        exp_ax.spines["top"].set_visible(False)
        exp_ax.spines["right"].set_visible(False)
        exp_ax.set_xticks([])
        exp_ax.set_yticks([])
        exp_ax.set_xlabel(contexts[2+sleuth_idx])

    if annotation_type == "transcript" or annotation_type == "gene":
        axes[0].set_ylabel("mRNA expression")

    elif annotation_type == "splicing":
        axes[0].set_ylabel("PSI")

    y_max = max(maxes)

    for ax in axes:
        ax.set_ylim(0, y_max*1.25)

    plt.subplots_adjust(wspace=0.25,hspace=1)
    
    if legend:

        legend_background = "#eaeaea"

        legend_elements = [Patch(label='Control',
                                 color=control_color, alpha=1),
                           Patch(label='Treatment',
                                 color=treatment_color, alpha=1)]

        legend = plt.legend(handles=legend_elements,
                            loc='upper left', bbox_to_anchor=(1, 1),)
        frame = legend.get_frame()
        frame.set_facecolor(legend_background)
        frame.set_edgecolor(legend_background)

    return axes
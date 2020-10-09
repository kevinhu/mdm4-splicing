import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc


def config_visuals():

    # disable SVG fonts
    plt.rcParams["ps.useafm"] = True
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    # set Seaborn style
    sns.set_style("ticks", {"axes.linewidth": 0.5})
    plt.rcParams["xtick.major.size"] = 8
    plt.rcParams["xtick.major.width"] = 1
    plt.rcParams["ytick.major.size"] = 8
    plt.rcParams["ytick.major.width"] = 1
    plt.rcParams["xtick.bottom"] = True
    plt.rcParams["ytick.left"] = True

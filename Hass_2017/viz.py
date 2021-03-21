from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from observable import observables

class Visualization(object):
    def __init__(self):
        self.cm = plt.cm.get_cmap("tab10")

        self.timecourse_options = [
            {
                "divided_by": 1,
                "xlim": (),
                "xticks": None,
                "xlabel": "Time",
                "ylim": (),
                "yticks": None,
                "ylabel": observables[i].replace("__", "\n").replace("_", " "),
                "exp_data": True,
                "legend_loc": None,
                "cmap": [self.cm.colors[j] for j in range(10)],
                "shape": Line2D.filled_markers,
                "dont_show": [],
            }
            for i, _ in enumerate(observables)
        ]

        self.multiplot_options = {
            "fig_name": "multiplot_observables",
            "observables": [],
            "condition": None,
            "xlim": (),
            "xticks": None,
            "xlabel": "Time",
            "ylim": (),
            "yticks": None,
            "ylabel": "",
            "cmap": [self.cm.colors[j] for j in range(10)],
            "shape": Line2D.filled_markers,
        }

        self.sensitivity_options = {
            "figsize": (12, 5),
            "width": 0.3,
            "legend_loc": "upper left",
            "cmap": [self.cm.colors[j] for j in range(10)],
        }
    
    def get_timecourse_options(self):
        for i, _ in enumerate(observables):
            self.timecourse_options[i]["xlim"] = (-10, 250)
            self.timecourse_options[i]["xticks"] = [50 * i for i in range(5)]
            self.timecourse_options[i]["xlabel"] = "Time (min)"

        return self.timecourse_options

    def multiplot_observables(self):

        return self.multiplot_options

    @staticmethod
    def set_timecourse_rcParams():
        """figure/simulation"""
        plt.rcParams["font.size"] = 18
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["xtick.major.width"] = 1.5
        plt.rcParams["ytick.major.width"] = 1.5
        plt.rcParams["lines.linewidth"] = 1.8
        plt.rcParams["lines.markersize"] = 12
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_param_range_rcParams():
        """figure/param_range"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def set_sensitivity_rcParams():
        """figure/sensitivity"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        """
        if name == 'MP':
            return 'Per mRNA'
        elif name == 'MC':
            return 'Cry mRNA'
        elif name == 'MB':
            return 'Bmal1 mRNA'
        """
        return name
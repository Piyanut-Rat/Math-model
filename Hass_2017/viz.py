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
import sys
import os
#from IPython.display import Image, display_png
sys.path.append("Hass_2017")

from biomass import run_simulation

import Hass_2017

model = Hass_2017.create()
#run_simulation(model, viz_type='original')
run_simulation(model, viz_type="original", save_format="png",show_all=True)
#run_simulation(model, viz_type="average", save_format="png")
#run_simulation(model, viz_type="best", save_format="png")
#run_simulation(model, viz_type="experiment", save_format="png")
'''
for observable in model.obs:
    with open(
        os.path.join(
            model.path,
            "figure",
            "simulation",
            "original",
            f"{observable}.png",
        ),
        mode="rb",
    ) as f:
        display_png(Image(f.read()))

'''
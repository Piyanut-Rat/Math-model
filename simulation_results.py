import sys
import os
sys.path.append("Hass_2017")

from biomass import run_simulation

import Hass_2017

model = Hass_2017.create()
#run_simulation(model, viz_type='original') #defult pdf format
run_simulation(model, viz_type="original", save_format="png")
#run_simulation(model, viz_type="original", save_format="png",show_all=True)
#run_simulation(model, viz_type="average", save_format="png")
#run_simulation(model, viz_type="best", save_format="png")
#run_simulation(model, viz_type="experiment", save_format="png")

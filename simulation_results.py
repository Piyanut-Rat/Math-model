from biomass import run_simulation

import Hass_2017
model = Hass_2017.create()
run_simulation(model, viz_type='original')


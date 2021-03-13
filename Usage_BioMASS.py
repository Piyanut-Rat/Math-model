#Model construction
import Hass_2017
Hass_2017.show_info()
model =  Hass_2017.create()

#Parameter estimation
from biomass import optimize
    # Estimate 10 parameter sets simultaneously
optimize(
   model=model, start=1, end=10, options={
      "popsize": 5,
      "max_generation": 1000,
      "allowable_error": 0.5,
      "local_search_method": "mutation",
      "n_children": 200
   }
)

#Visualization of simulation results
from biomass import run_simulation
run_simulation(model, viz_type="average", show_all=False, stdev=True)

#Sensitivity analysis
from biomass import run_analysis
run_analysis(model, target="reaction", metric='integral', style='barplot')





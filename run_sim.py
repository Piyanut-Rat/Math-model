import sys
import os
#from IPython.display import Image, display_png
sys.path.append("Hass_2017")

from simulation import Simulation
#import HHass2017 
import plot_func

def run_simulation():
    sim = Simulation()

    plot_func.timecourse(sim)
    
    
if __name__ == "__main__":
    run_simulation()
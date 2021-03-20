import numpy as np

from biomass.dynamics.solver import get_steady_state, solve_ode

from name2idx import parameters as C
from name2idx import species as V
from set_model import DifferentialEquation

observables = [
    "FACS_EGFR",
    "FACS_ErbB2",
    'FACS_ErbB3',
    'FACS_IGF1R',
    'FACS_Met',

    'tEGFR_au ',
    'tErbB2_au ',
    'tErbB3_au',
    'tIGF1R_au',

    'pEGFR_au',
    'pErbB2_au',
    'pErbB3_au',
    'pIGF1R_au',
    'pMet_au',

    'pMEK_au',
    'pERK_au',

    'pAKT_au',

    'pS6K1_au',
    'pS6_au',
]

class NumericalSimulation(DifferentialEquation):
    # inheritance : make new class (subclass) from old class(superclass)
    def __init__(self):
        super().__init__(pertubation = {})
        self.normalization = {}

    t = range(240) # unit; sec.

    # experimental conditions
    conditions = ['control']
    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation = {}):
        if _perturbation:
            self.pertubation = _perturbation
        #get steady state
        x[C.Ligand] = x[C.no_ligand] # No ligand binding
        
        y0 = get_steady_state(self.diffeq, y0, tuple(x))
        
        if not y0:
            return False
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == 'EGF':
                x[C.Ligand] = x[C.EGF]
            elif condition == 'HRG':
                x[C.Ligand] = x[C.HRG]
            elif condition == 'IGF1':
                x[C.Ligand] = x[C.IGF1]
            elif condition == 'HGF':
                x[C.Ligand] = x[C.HGF]
            elif condition == 'BTC':
                x[C.Ligand] = x[C.BTC]
            
            
            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[observables.index('FACS_EGFR'), :, i] = np.log10(
                    sol.y[V.EGFR, :] + sol.y[V.EGFR_EGF, :] + sol.y[V.EGFR_BTC, :] + 2*sol.y[V.EGFR_BTC, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB13, :] + sol.y[V.pMetEGFR, :]
                )
                self.simulations[observables.index('FACS_ErbB2'), :, i] = np.log10(
                    sol.y[V.ErbB2, :] + 2*sol.y[V.pErbB2, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB32, :]
                )
                self.simulations[observables.index('FACS_ErbB3'), :, i] = np.log10(
                    sol.y[V.ErbB3, :] + sol.y[V.ErbB3_HRG, :] + 2*sol.y[V.pErbB3d, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB13, :] + sol.y[V.pMetErbB3]
                )
                self.simulations[observables.index('FACS_IGF1R'), :, i] = np.log10(
                    sol.y[V.IGF1R, :] + sol.y[V.IGF1R_IGF1, :] + 2*sol.y[V.pIGF1Rd, :]
                )
                self.simulations[observables.index('FACS_Met'), :, i] = np.log10(
                    sol.y[V.Met, :] + 2*sol.y[V.pMetd, :] + sol.y[V.pMetEGFR, :] + sol.y[V.pMetErbB3, :]
                )

                self.simulations[observables.index('tEGFR_au'), :, i] = np.log10(
                    offset_tEGFR_Cellline + scale_tEGFR_Cellline* (sol.y[V.EGFR, :] + sol.y[V.EGFR_EGF, :] + sol.y[V.EGFR_BTC, :] + sol.y[V.EGFRi, :] + 
                    2*sol.y[V.pEGFRd, :] + 2*sol.y[V.pEGFRi, :] + 2*sol.y[V.pEGFRi_ph, :] + sol.y[V.pErbB12i, :] + sol.y[V.pErbB13i, :] + 
                    sol.y[V.pErbB12, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB12i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetEGFR, :] + 
                    sol.y[V.pMetEGFRi, :] + sol.y[V.pMetEGFRi_ph, :])  
                )
                self.simulations[observables.index('tErbB2_au'), :, i] = np.log10(
                    offset_tErbB2_Cellline + scale_tErbB2_Cellline *(sol.y[V.ErbB2, :] + sol.y[V.ErbB2i, :] + 2*sol.y[V.pErbB2, :] + 2*sol.y[V.pErbB2i, :] + 
                    2*sol.y[V.pErbB2i_ph, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB12i, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB12i_ph, :] +
                    sol.y[V.pErbB32i_ph, :])
                )
                self.simulations[observables.index('tErbB3_au'), :, i] = np.log10(
                    offset_tErbB3_Cellline + scale_tErbB3_Cellline *(sol.y[V.ErbB3, :] + sol.y[V.ErbB3_HRG, :] + sol.y[V.ErbB3i, :] + 2*sol.y[V.pErbB3d, :] +
                    2*sol.y[V.pErbB3i, :] + 2*sol.y[V.pErbB3i_ph, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB13i, :] + 
                    sol.y[V.pErbB32i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetErbB3, :] + sol.y[V.pMetErbB3i, :] + sol.y[V.pMetErbB3i_ph, :])
                )
                self.simulations[observables.index('tIGF1R_au'), :, i] = np.log10(
                    offset_tIGF1R_Cellline + scale_tIGF1R_Cellline *(sol.y[V.IGF1R, :] + sol.y[V.IGF1R_IGF1, :] + sol.y[V.IGF1Ri, :] + 2*sol.y[V.pIGF1Rd, :] + 
                    2*sol.y[V.pIGF1Ri, :] + 2*sol.y[V.pIGF1Ri_ph, :])
                )

                self.simulations[observables.index('pEGFR_au'), :, i] = np.log10(
                    offset_pEGFR_Cellline + scale_pEGFR_Cellline * (2*sol.y[V.pEGFRd, :] + 2*sol.y[V.pEGFRi, :] + 2*sol.y[V.pEGFRi_ph, :] + sol.y[V.pErbB12i, :] + 
                    sol.y[V.pErbB13i, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB12i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetEGFR, :] + 
                    sol.y[V.pMetEGFRi, :] + sol.y[V.pMetEGFRi_ph, :])
                )
                self.simulations[observables.index('pErbB2_au'), :, i] = np.log10(
                    offset_pErbB2_Cellline + scale_pErbB2_Cellline * (2*sol.y[V.pErbB2, :] + 2*sol.y[V.pErbB2i, :] + 2*sol.y[V.pErbB2i_ph, :] + sol.y[V.pErbB12, :] + 
                    sol.y[V.pErbB32, :] + sol.y[V.pErbB12i, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB12i_ph, :] + sol.y[V.pErbB32i_ph, :])
                )
                self.simulations[observables.index('pErbB3_au'), :, i] = np.log10(
                    offset_pErbB3_Cellline + scale_pErbB3_Cellline * (2*sol.y[V.pErbB3d, :] + 2*sol.y[V.pErbB3i, :] + 2*sol.y[V.pErbB3i_ph, :] + sol.y[V.pErbB32, :] + 
                    sol.y[V.pErbB13, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB13i, :] + sol.y[V.pErbB32i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetErbB3, :] + 
                    sol.y[V.pMetErbB3i, :] + sol.y[V.pMetErbB3i_ph, :])
                )
                self.simulations[observables.index('pIGF1R_au'), :, i] = np.log10(
                    offset_pIGF1R_Cellline + scale_pIGF1R_Cellline * (2*sol.y[V.pIGF1Rd, :] + 2*sol.y[V.pIGF1Ri, :] + 2*sol.y[V.pIGF1Ri_ph, :])
                )
                self.simulations[observables.index('pMet_au'), :, i] = np.log10(
                    offset_pMet_Cellline + scale_pMet_Cellline * (2*sol.y[V.pMetd, :] + 2*sol.y[V.pMeti, :] + 2*sol.y[V.pMeti_ph, :] + 
                    sol.y[V.pMetEGFR, :] + sol.y[V.pMetEGFRi, :] + sol.y[V.pMetEGFRi_ph, :] + sol.y[V.pMetErbB3, :] + sol.y[V.pMetErbB3i, :] + 
                    sol.y[V.pMetErbB3i_ph, :])
                )

                self.simulations[observables.index('pMEK_au'), :, i] = np.log10(
                    offset_pMEK_Cellline + scale_pMEK_Cellline * sol.y[V.pMEK, :]
                )
                self.simulations[observables.index('pERK_au'), :, i] = np.log10(
                    offset_pERK_Cellline + scale_pERK_Cellline * sol.y[V.pERK, :]
                )

                self.simulations[observables.index('pAKT_au'), :, i] = np.log10(
                    offset_pAKT_Cellline + scale_pAKT_Cellline * sol.y[V.pAKT, :]
                )

                self.simulations[observables.index('pS6K1_au'), :, i] = np.log10(
                    offset_pS6K1_Cellline + scale_pS6K1_Cellline * sol.y[V.pS6K1, :] 
                )
                self.simulations[observables.index('pS6_au'), :, i] = np.log10(
                    offset_pS6_Cellline + scale_pS6_Cellline * sol.y[V.pS6, :]
                )

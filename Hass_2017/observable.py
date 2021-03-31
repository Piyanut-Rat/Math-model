import numpy as np

from biomass.dynamics.solver import get_steady_state, solve_ode

from name2idx import parameters as C
from name2idx import species as V
from set_model import DifferentialEquation, param_values, initial_values

observables = [
    'FACS_EGFR',
    'FACS_ErbB2',
    'FACS_ErbB3',
    'FACS_IGF1R',
    'FACS_Met',

    'tEGFR_au',
    'tErbB2_au',
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

        for observable in observables:
            self.normalization[observable] = {'timepoint': None, 'condition': []}

    t = range(250) # unit; sec.

    #x = param_values()
    #y0 = initial_values()
   

    # experimental conditions
    conditions = ['control','0.156','0.625', '2.500', '10.000']

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation = {}):
        if _perturbation:
            self.pertubation = _perturbation
        #get steady state
        x[C.Ligand] = x[C.no_ligand] # No ligand binding
        #y0[V.dose_EGF] = 0
        #y0[V.dose_HGF] = 0
        #y0[V.dose_IGF1] = 0
        #y0[V.dose_HRG] = 0
        #y0 = get_steady_state(self.diffeq, y0, tuple(x))
        y0 = initial_values()
        sol = solve_ode(self.diffeq, y0, self.t, tuple(x))
        if sol is None:
            return False

        else:
            y0 = sol.y[:, -1].tolist()

        # add ligand
        
        for i, condition in enumerate(self.conditions):
            
            if condition == 'control':
                y0[V.dose_EGF] = 0

            elif condition == '0.156':
                y0[V.dose_EGF] = 0.156*x[C.scale_Ligand]

            elif condition == '0.625':
                y0[V.dose_EGF] = 0.625*x[C.scale_Ligand]

            elif condition == '2.500':
                y0[V.dose_EGF] = 2.5*x[C.scale_Ligand]

            elif condition == '10.000':
                y0[V.dose_EGF] = 10*x[C.scale_Ligand]


            
            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                '''
                self.simulations[observables.index('FACS_EGFR'), :, i] = np.log10(sol.y[V.EGFR, :] + sol.y[V.EGFR_EGF, :] + sol.y[V.EGFR_BTC, :] + 2*sol.y[V.EGFR_BTC, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB13, :] + sol.y[V.pMetEGFR, :]
                )
                self.simulations[observables.index('FACS_ErbB2'), :, i] = np.log10(sol.y[V.ErbB2, :] + 2*sol.y[V.pErbB2, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB32, :]
                )
                self.simulations[observables.index('FACS_ErbB3'), :, i] = np.log10(sol.y[V.ErbB3, :] + sol.y[V.ErbB3_HRG, :] + 2*sol.y[V.pErbB3d, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB13, :] + sol.y[V.pMetErbB3]
                )
                self.simulations[observables.index('FACS_IGF1R'), :, i] = np.log10(sol.y[V.IGF1R, :] + sol.y[V.IGF1R_IGF1, :] + 2*sol.y[V.pIGF1Rd, :]
                )
                self.simulations[observables.index('FACS_Met'), :, i] = np.log10(sol.y[V.Met, :] + 2*sol.y[V.pMetd, :] + sol.y[V.pMetEGFR, :] + sol.y[V.pMetErbB3, :]
                )
                
                self.simulations[observables.index('tEGFR_au'), :, i] = np.log10(
                    x[C.offset_tEGFR_CelllineH322M] + x[C.scale_tEGFR_CelllineH322M]* (sol.y[V.EGFR, :] + sol.y[V.EGFR_EGF, :] + sol.y[V.EGFR_BTC, :] + sol.y[V.EGFRi, :] + 
                    2*sol.y[V.pEGFRd, :] + 2*sol.y[V.pEGFRi, :] + 2*sol.y[V.pEGFRi_ph, :] + sol.y[V.pErbB12i, :] + sol.y[V.pErbB13i, :] + 
                    sol.y[V.pErbB12, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB12i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetEGFR, :] + 
                    sol.y[V.pMetEGFRi, :] + sol.y[V.pMetEGFRi_ph, :])  
                )
                self.simulations[observables.index('tErbB2_au'), :, i] = np.log10(
                    x[C.offset_tErbB2_CelllineH322M] + x[C.scale_tErbB2_CelllineH322M] *(sol.y[V.ErbB2, :] + sol.y[V.ErbB2i, :] + 2*sol.y[V.pErbB2, :] + 2*sol.y[V.pErbB2i, :] + 
                    2*sol.y[V.pErbB2i_ph, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB12i, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB12i_ph, :] +
                    sol.y[V.pErbB32i_ph, :])
                )
                self.simulations[observables.index('tErbB3_au'), :, i] = np.log10(
                    x[C.offset_tErbB3_CelllineH322M] + x[C.scale_tErbB3_CelllineH322M] *(sol.y[V.ErbB3, :] + sol.y[V.ErbB3_HRG, :] + sol.y[V.ErbB3i, :] + 2*sol.y[V.pErbB3d, :] +
                    2*sol.y[V.pErbB3i, :] + 2*sol.y[V.pErbB3i_ph, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB13i, :] + 
                    sol.y[V.pErbB32i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetErbB3, :] + sol.y[V.pMetErbB3i, :] + sol.y[V.pMetErbB3i_ph, :])
                )
                self.simulations[observables.index('tIGF1R_au'), :, i] = np.log10(
                    x[C.offset_tIGF1R_CelllineH322M] + x[C.scale_tIGF1R_CelllineH322M] *(sol.y[V.IGF1R, :] + sol.y[V.IGF1R_IGF1, :] + sol.y[V.IGF1Ri, :] + 2*sol.y[V.pIGF1Rd, :] + 
                    2*sol.y[V.pIGF1Ri, :] + 2*sol.y[V.pIGF1Ri_ph, :])
                )
                '''
                self.simulations[observables.index('pEGFR_au'), :, i] = np.log10(x[C.offset_pEGFR_CelllineH322M] + x[C.scale_pEGFR_CelllineH322M] * (2*sol.y[V.pEGFRd, :] + 2*sol.y[V.pEGFRi, :] + 2*sol.y[V.pEGFRi_ph, :] + sol.y[V.pErbB12i, :] + 
                    sol.y[V.pErbB13i, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB12i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetEGFR, :] + sol.y[V.pMetEGFRi, :] + sol.y[V.pMetEGFRi_ph, :])
                )
                self.simulations[observables.index('pErbB2_au'), :, i] = np.log10(x[C.offset_pErbB2_CelllineH322M] + x[C.scale_pErbB2_CelllineH322M] * (2*sol.y[V.pErbB2, :] + 2*sol.y[V.pErbB2i, :] + 2*sol.y[V.pErbB2i_ph, :] + sol.y[V.pErbB12, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB12i, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB12i_ph, :] + sol.y[V.pErbB32i_ph, :])
                )
                self.simulations[observables.index('pErbB3_au'), :, i] = np.log10(x[C.offset_pErbB3_CelllineH322M] + x[C.scale_pErbB3_CelllineH322M] * (2*sol.y[V.pErbB3d, :] + 2*sol.y[V.pErbB3i, :] + 2*sol.y[V.pErbB3i_ph, :] + sol.y[V.pErbB32, :] + sol.y[V.pErbB13, :] + sol.y[V.pErbB32i, :] + sol.y[V.pErbB13i, :] + sol.y[V.pErbB32i_ph, :] + sol.y[V.pErbB13i_ph, :] + sol.y[V.pMetErbB3, :] + sol.y[V.pMetErbB3i, :] + sol.y[V.pMetErbB3i_ph, :])
                )
                '''
                self.simulations[observables.index('pIGF1R_au'), :, i] = np.log10(
                    x[C.offset_pIGF1R_CelllineH322M] + x[C.scale_pIGF1R_CelllineH322M] * (2*sol.y[V.pIGF1Rd, :] + 2*sol.y[V.pIGF1Ri, :] + 2*sol.y[V.pIGF1Ri_ph, :])
                )
                self.simulations[observables.index('pMet_au'), :, i] = np.log10(
                    x[C.offset_pMet_CelllineH322M] + x[C.scale_pMet_CelllineH322M] * (2*sol.y[V.pMetd, :] + 2*sol.y[V.pMeti, :] + 2*sol.y[V.pMeti_ph, :] + 
                    sol.y[V.pMetEGFR, :] + sol.y[V.pMetEGFRi, :] + sol.y[V.pMetEGFRi_ph, :] + sol.y[V.pMetErbB3, :] + sol.y[V.pMetErbB3i, :] + 
                    sol.y[V.pMetErbB3i_ph, :])
                )

                self.simulations[observables.index('pMEK_au'), :, i] = np.log10(x[C.offset_pMEK_CelllineH322M] + x[C.scale_pMEK_CelllineH322M] * sol.y[V.pMEK, :]
                )
                '''
                self.simulations[observables.index('pERK_au'), :, i] = np.log10(x[C.offset_pERK_CelllineH322M] + x[C.scale_pERK_CelllineH322M] * sol.y[V.pERK, :]
                )

                self.simulations[observables.index('pAKT_au'), :, i] = np.log10(x[C.offset_pAKT_CelllineH322M] + x[C.scale_pAKT_CelllineH322M] * sol.y[V.pAKT, :]
                )
                '''
                self.simulations[observables.index('pS6K1_au'), :, i] = np.log10(x[C.offset_pS6K1_CelllineH322M] + x[C.scale_pS6K1_CelllineH322M] * sol.y[V.pS6K1, :] 
                )
                '''
                self.simulations[observables.index('pS6_au'), :, i] = np.log10(x[C.offset_pS6_CelllineH322M] + (x[C.scale_pS6_CelllineH322M]* sol.y[V.pS6, :])
                )
                

class ExperimentalData(object):
    """
    Set experimental data.
    Attributes
    ----------
    experiments : list of dict
        Time series data.
    error_bars : list of dict
        Error bars to show in figures.
    """

    def __init__(self):
        self.experiments = [None] * len(observables)
        self.error_bars = [None] * len(observables)

    def set_data(self):
        pass

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return []

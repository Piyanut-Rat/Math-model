#from name2idx import C, V
from name2idx import parameters as C
from name2idx import species as V

class DifferentialEquation(object):
    def __init__(self, pertubation): #dunder-> double underscores
        self.pertubation = pertubation

    #Defined Model
    def diffeq(self, t, y, *x):
        v = {}
        # Rate equations
        ## Production
        v[1] = x[C.EGFR_prod]
        v[2] = x[C.ErbB2_prod]
        v[3] = x[C.ErbB3_prod]
        v[4] = x[C.IGF1R_prod]
        v[5] = x[C.Met_prod]

        ## Ligand binding
        v[6] = y[V.EGFR]* x[C.EGFR_lig_binding]* y[V.dose_EGF]
        v[7] = y[V.EGFR_EGF]* x[C.EGFR_lig_binding]* x[C.EGF_kD]
        v[8] = y[V.EGFR]* x[C.EGFR_BTC_binding]* y[V.dose_BTC]
        v[9] = y[V.EGFR_BTC]* x[C.EGFR_BTC_binding]* x[EGF_kD]
        v[10] = y[V.ErbB3]* x[C.ErbB3_lig_binding]* y[V.dose_HRG]
        
        v[11] = y[V.ErbB3_HRG]* x[C.ErbB3_lig_binding]* x[C.HRG_kD]









'''
# call class
if __name__ == '__main__':
    D = DifferentialEquation(3) #create new instance
    print(D.pertubation)
'''

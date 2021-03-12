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
        v[12] = y[V.IGF1R]* x[C.IGF1R_lig_binding]* y[V.dose_IGF1]
        v[13] = y[V.IGF1R_IGF1]* x[C.IGF1R_lig_binding]* x[C.IGF1_kD]
        v[14] = y[V.Met]* x[C.Met_lig_binding]* y[V.dose_HGF]
        v[15] = x[C.HGF_kD]* y[V.Met_HGF]* x[C.Met_lig_binding]

        ##Dimerization
        v[16] = (y[V.EGFR_EGF]**2)* x[C.EGFR_dimerize]
        v[17] = (y[V.EGFR_BTC]**2)* x[C.EGFR_BTC_dimerize]

        v[18] = (y[V.ErbB2]**2)* x[C.ErbB2_dimerize]

        v[19] = (y[V.ErbB3_HRG]**2)* x[C.ErbB3_dimerize]

        v[20] = (y[V.IGF1R_IGF1]**2)* x[C.IGF1R_dimerize]

        v[21] = y[V.EGFR_EGF]* x[C.EGFR_ErbB2_dimerize]* y[V.ErbB2]
        v[22] = y[V.EGFR_BTC]* x[C.EGFR_ErbB2_BTC_dimerize]* y[V.ErbB2]

        v[23] = y[V.EGFR_EGF]* x[C.EGFR_ErbB3_dimerize]* y[V.ErbB3_HRG]
        v[24] = y[V.EGFR_BTC]* x[C.EGFR_ErbB3_BTC_dimerize]* y[V.ErbB3_HRG]
        v[25] = y[V.EGFR_BTC]* x[V.EGFR_ErbB3_dimerize_noHRG]* y[V.ErbB3]

        v[26] = y[V.ErbB2]* x[C.ErbB2_ErbB3_dimerize]* y[V.ErbB3_HRG]

        v[27] = (y[V.Met_HGF]**2)* x[C.Met_dimerize]

        v[28] = y[V.ErbB3_HRG]* y[V.Met]* x[C.Met_ErbB3_dimerize]

        v[29] = y[V.ErbB3_HRG]* y[V.Met_HGF]* x[C.Met_lig_ErbB3_dimerize]

        v[30] = y[V.EGFR_EGF]* x[C.Met_EGFR_dimerize]* y[V.Met_HGF]
        v[31] = y[V.EGFR_BTC]* x[C.Met_EGFR_BTC_dimerize]* y[V.Met_HGF]

        ##Basal phosphorylations
        v[32] = (y[V.EGFR]**2)* x[C.EGFR_basal_activation]
        v[33] = (y[V.ErbB3]**2)* x[C.ErbB3_basal_activation]
        v[34] = (y[V.IGF1R]**2)* x[C.IGF1R_basal_activation]

        v[35] = y[V.EGFR]* x[C.EGFR_ErbB2_basal_act]* y[V.ErbB2]
        v[36] = y[V.EGFR]* x[C.EGFR_ErbB3_basal_act]* y[V.ErbB3]
        v[37] = y[V.ErbB2]* y[V.ErbB3]* x[C.ErbB3_ErbB2_basal_act]
        v[38] = y[V.ErbB3]* y[V.Met]* x[C.Met_ErbB3_basal_act]
        v[39] = (y[V.Met]**2)*x[C.Met_basal_act]
        v[40] = y[V.EGFR]* y[V.Met]* x[C.Met_EGFR_basal_act]

        ##Phosphorylation
        ###EGFR homodimers
        v[41] = x[C.pEGFR_internalize]* y[V.pEGFRd]
        v[42] = x[C.pEGFR_degradation]* y[V.pEGFRi]
        v[43] = y[V.RTKph]* x[C.pEGFR_phosphatase_binding]* y[X.pEGFRi]
        v[44] = x[C.pEGFRi_dephosph]* y[V.pEGFRi_ph]
        v[45] = x[C.EGFR_basal_recycle]* y[V.EGFRi]

        ###ErbB2 homodimers
        v[46] = y[V.pErbB2]* x[C.pErbB2_internalize]
        v[47] = y[V.RTKph]* y[V.pErbB2i]* x[C.pErbB2i_phosphatase]
        v[48] = x[C.pErbB2i_dephosph]* y[V.pErbB2i_ph]
        v[49] = x[C.pErbB2_degradation]* y[V.pErbB2i]
        v[50] = x[C.ErbB2_recycle]* y[V.ErbB2i]

        ###ErbB3 homodimers
        v[51] = x[C.pErbB3_internalize]* y[V.pErbB3d]
        v[52] = x[C.pErbB3_degradation]* y[V.pErbB3i]
        v[53] = y[V.RTKph]* y[V.pErbB3i]* x[C.pErbB3i_phosphatase]
        v[54] = x[C.pErbB3i_dephosph]* y[V.pErbB3i_ph]
        v[55] = x[C.ErbB3_basal_recycle]* y[V.ErbB3i]

        ###IGF1R homodimers
        v[56] = x[C.pIGF1R_internalize]* y[V.pIGF1Rd]
        v[57] = x[C.pIGF1R_degradation]* y[V.pIGF1Ri]
        v[58] = y[V.RTKph]* y[V.pIGF1Ri]* x[C.pIGF1Ri_phosphatase]
        v[59] = x[C.pIGF1Ri_dephosph]* y[V.pIGF1Ri_ph]
        v[60] = x[C.IGF1R_basal_recycle]* y[V.IGF1Ri]

        ###EGFR-ErbB2 heterodimers
        v[61] = y[V.pErbB12]* x[C.pErbB12_internalize]
        v[62] = x[C.pErbB12_degradation]* y[V.pErbB12i]
        v[63] = y[V.RTKph]* y[V.pErbB12i]* x[C.pErbB12i_phosphatase]
        v[64] = x[C.pErbB12i_dephosph]* y[V.pErbB12i_ph]

        ###ErbB2 - ErbB3 heterodimers
        v[65] = y[V.pErbB32]* x[C.pErbB32_internalize]
        v[66] = x[C.pErbB32_degradation]* y[V.pErbB32i]
        v[67] = y[V.RTKph]* y[V.pErbB32i]* x[C.pErbB32i_phosphatase]
        v[68] = x[C.pErbB32i_dephosph]* y[V.pErbB32i_ph]

        ###EGFR-ErbB3 heterodimers
        v[69] = y[V.pErbB13]* x[C.pErbB13_internalize]
        v[70] = x[C.pErbB13_degradation]* y[V.pErbB13i]
        v[71] = y[V.RTKph]* y[V.pErbB13i]* x[C.pErbB13i_phosphatase]
        v[72] = x[C.pErbB13i_dephosph]* y[V.pErbB13i_ph]

        ###MET homodimers
        v[73] = x[C.pMet_internalize]* y[V.pMetd]
        v[74] = x[C.pMet_degradation]* y[V.pMeti]
        v[75] = y[V.RTKph]* y[V.pMeti]* x[C.pMeti_phosphatase]
        v[76] = x[C.pMeti_dephosph]* y[V.pMeti_ph]
        v[77] = x[C.Met_recycle]* y[V.Meti]

        ###MET-ErbB3 heterodimers
        v[78] = y[V.pMetErbB3]* x[C.pMetErbB3_internalize]
        v[79] = x[C.pMetErbB3_degradation]* y[V.pMetErbB3i]
        v[80] = y[V.RTKph]* y[V.pMetErbB3i]* x[C.pMetErbB3i_phosphatase]
        v[81] = x[C.pMetErbB3i_dephosph]* y[V.pMetErbB3i_ph]

        ###MET-EGFR heterodimers
        v[82] = y[V.pMetEGFR]* x[C.pMetEGFR_internalize]
        v[83] = x[C.pMetEGFR_degradation]* y[V.pMetEGFRi]
        v[84] = y[V.RTKph]* y[V.pMetEGFRi]* x[C.pMetEGFRi_phosphatase]
        v[85] = x[C.pMetEGFRi_dephosph]* y[V.pMetEGFRi_ph]

        ##DOWNSTREAM
        ###ERK branch
        v[86] = 
        










'''
# call class
if __name__ == '__main__':
    D = DifferentialEquation(3) #create new instance
    print(D.pertubation)
'''

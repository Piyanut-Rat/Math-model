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
        v[86] = (y[V.MEK]* x[C.MEK_phosphorylation_pEGFR]* y[V.pEGFRd])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[87] = (y[V.MEK]* x[C.MEK_phosphorylation_pErbB12]* y[V.pErbB12])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[88] = (y[V.MEK]* x[C.MEK_phosphorylation_pErbB13]* y[V.pErbB13])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[89] = (y[V.MEK]* x[C.MEK_phosphorylation_pErbB32]* y[V.pErbB32])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[90] = (y[V.MEK]* x[C.MEK_phosphorylation_pMetd]* y[V.pMetd])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[91] = (y[V.MEK]* x[C.MEK_phosphorylation_pMetEGFR]* y[V.pMetEGFR])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[92] = (y[V.MEK]* x[C.MEK_phosphorylation_pIGF1R]* y[V.pIGF1Rd])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[93] = (y[V.MEK]* x[C.MEK_phosphorylation_pMetErbB3]* y[V.pMetErbB3])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[94] = (y[V.MEK]* x[C.MEK_internIGF1R_effect]* x[C.MEK_phosphorylation_pIGF1R]* y[V.pIGF1Ri])/(1+x[C.feedback_pAKT]* x[C.init_pAKT]+x[C.feedback_pERK]* y[V.pERK])
        v[95] = y[V.pMEK]* x[C.pMEK_dephosphorylation]

        v[96] = y[V.ERK]* x[C.ERK_phosphorylation_pMEK]* y[V.pMEK]
        v[97] = y[V.pERK]* x[C.pERK_dephosphorylation]
        
        ###AKT branch
        v[98] = (y[V.AKT]* x[C.AKT_activation_pEGFR]* y[V.pEGFRd])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[99] = (y[V.AKT]* x[C.AKT_activation_pErbB12]* y[V.pErbB12])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[100] = (y[V.AKT]* x[C.AKT_activation_pErbB13]* y[V.pErbB13])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[101] = (y[V.AKT]* x[C.AKT_activation_pErbB32]* y[V.pErbB32])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[102] = (y[V.AKT]* x[C.AKT_activation_pMetEGFR]* y[V.pMetEGFR])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[103] = (y[V.AKT]* x[C.AKT_activation_pMetd]* y[V.pMetd])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[104] = (y[V.AKT]* x[C.AKT_activation_pIGF1R]* y[V.pIGF1Rd])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[105] = (y[V.AKT]* x[C.AKT_activation_pMetErbB3]* y[V.pMetErbB3])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[106] = (y[V.AKT]* x[C.AKT_activation_pIGF1R]* x[C.AKT_internIGF1R_effect]* y[V.pIGF1Ri])/(x[C.feedback_pERK_on_AKT]* y[V.pERK] + x[C.feedback_pS6K1]* y[V.pS6K1] + 1)
        v[107] = y[V.pAKT]* x[C.pAKT_deactivation]

        ###GSK3-S6K-S6 
        v[108] = y[V.S6K1]* x[C.S6K1_phosphorylation_pAKT]* y[V.pAKT]
        v[109] = y[V.S6K1]* x[C.S6K1_phosphorylation_pERK]* y[V.pERK]
        v[110] = y[V.pS6K1]* x[C.pS6K1_dephosphorylation]

        v[111] = y[V.S6]* x[C.S6_phosphorylation_pS6K1]* y[V.pS6K1]
        v[112] = y[V.S6]* x[C.S6_phosphorylation_pERK]* y[V.pERK]
        v[113] = y[V.pS6]* x[C.pS6_dephosphorylation]

        #ODE of the law of the mass action
        dydt = [0]* V.NUM
        dydt[V.dose_EGF] = -v[6] + v[7]
        dydt[V.dose_HGF] = -v[14] + v[15]
        dydt[V.RTKph] = -v[43] + v[44] - v[47] + v[48] - v[53] + v[54] - v[58] + v[59] - v[63] + v[64] - v[67] + v[68] - v[71] + v[72] - v[75] + v[76] - v[80] + v[81] - v[84] + v[85] 
        dydt[V.dose_IGF1] = - v[12] + v[13]
        dydt[V.dose_HRG] = -v[10] + v[11]
        dydt[V.dose_BTC] = -v[8] + v[9]
        dydt[V.EGFR] = v[1] - v[6] + v[7] - v[8] + v[9] - 2*v[32] - v[35] - v[36] - v[40] + v[45]
        dydt[V.EGFR_EGF] = v[6] - v[7] - 2*v[16] - v[21] - v[23] - v[30]
        dydt[V.EGFR_BTC] = v[8] - v[9] - 2*v[17] - v[22] - v[24] - v[25] - v[31]
        dydt[V.pEGFRd] = v[16] + v[17] + v[32] - v[41]
        dydt[V.pEGFRi] = v[41] - v[42] - v[43]
        dydt[V.pEGFRi_ph] = v[43] - v[44]
        dydt[V.EGFRi] = +2*v[44] - v[45] + v[64] + v[72] + v[85]
        dydt[V.ErbB2] = v[2] - 2*v[18] - v[21] - v[22] - v[26] - v[35] - v[37] + v[50]
        dydt[V.pErbB2] = v[18] - v[46]
        dydt[V.pErbB2i] = v[46] - v[47] - v[49]
        dydt[V.ErbB2i] = +2*v[48] - v[50] + v[64] + v[68]
        dydt[V.pErbB2i_ph] = v[47] - v[48]
        dydt[V.pErbB12] = v[21] + v[22] + v[35] - v[61]
        dydt[V.pErbB12i] = v[61] - v[62] - v[63]
        dydt[V.pErbB12i_ph] = v[63] - v[64]
        dydt[V.ErbB3] = v[3] - v[10] + v[11] - v[25] - 2*v[33] - v[36] - v[37] - v[38] + v[55]
        dydt[V.ErbB3_HRG] = v[10] - v[11] - 2*v[19] - v[23] - v[24] - v[26] - v[28] - v[29]
        dydt[V.pErbB3d] = v[19] + v[33] - v[51]
        dydt[V.pErbB3i] = v[51] - v[52] - v[53]
        dydt[V.pErbB3i_ph] = v[53] - v[54]
        dydt[V.ErbB3i] = +2*v[54] - v[55] + v[68] + v[72] + v[81]
        dydt[V.pErbB13] = v[23] + v[24] + v[25] + v[36] - v[69]
        dydt[V.pErbB13i] = v[69] - v[70] - v[71]
        dydt[V.pErbB13i_ph] = v[71] - v[72]
        dydt[V.pErbB32] = v[26] + v[37] - v[65]
        dydt[V.pErbB32i] = v[65] - v[66] - v[67]
        dydt[V.pErbB32i_ph] = v[67] - v[68]
        dydt[V.IGF1R] = v[4] - v[12] + v[13] -2*v[34] + v[60]
        dydt[V.IGF1R_IGF1] = v[12] - v[13] -2*v[20]
        dydt[V.pIGF1Rd] = v[20] + v[34] - v[56]
        dydt[V.pIGF1Ri] = v[56] - v[57] - v[58]
        dydt[V.pIGF1Ri_ph] = v[58] - v[59]
        dydt[V.IGF1Ri] = +2*v[59 ] - v[60]
        dydt[V.Met] = v[5] - v[14] + v[15] - v[28] - v[38] - 2*v[39] - v[40] + v[77]
        dydt[V.Met_HGF] = v[14] - v[15] - 2*v[27] - v[29] - v[30] - v[31]
        dydt[V.pMetd] = v[27] + v[39] - v[73]
        dydt[V.pMeti] = v[73] - v[74] - v[75]
        dydt[V.pMeti_ph] = v[75] - v[76]
        dydt[V.Meti] = +2*v[76] - v[77] + v[81] + v[85]
        dydt[V.pMetErbB3] = v[28] + v[29] + v[38] - v[78]
        dydt[V.pMetErbB3i] = v[78] - v[79] - v[80]
        dydt[V.pMetErbB3i_ph] = v[80] - v[81]
        dydt[V.pMetEGFR] = v[30] + v[31] + v[40] - v[82]
        dydt[V.pMetEGFRi] = v[82] - v[83] - v[84]
        dydt[V.pMetEGFRi_ph] = v[84] - v[85]
        dydt[V.MEK] = -v[86] - v[87] - v[88] - v[89] - v[90] - v[91] - v[92] - v[93] - v[94] + v[95]
        dydt[V.pMEK] = v[86] + v[87] + v[88] + v[89] + v[90] + v[91] + v[92] + v[93] + v[94] - v[95]
        dydt[V.ERK] = -v[96] + v[97]
        dydt[V.pERK] = v[96] - v[97]
        dydt[V.AKT] = -v[98] - v[99] - v[100] - v[101] - v[102] - v[103] - v[104] - v[105] - v[106] + v[107]
        dydt[V.pAKT] = v[98] + v[99] + v[100] + v[101] + v[102] + v[103] + v[104] + v[105] + v[106] - v[107]
        dydt[V.S6K1] = -v[108] - v[109] + v[110]
        dydt[V.pS6K1] = v[108] + v[109] - v[110]
        dydt[V.S6] = -v[111] - v[112] + v[113]
        dydt[V.pS6] = v[111] + v[112] - v[113]

        return dydt

def param_values():
    x = [0]*C.NUM
    # 1-82
    x[C.AKT_activation_pEGFR] = -5.0000     #1.00*10**(−05) 
    x[C.AKT_activation_pErbB12] = -1.19400  #6.40*10**(−02) 
    x[C.AKT_activation_pErbB13] =  1.1178   #1.31*10**(+01)
    x[C.AKT_activation_pErbB32] = -0.2513   #5.61*10**(−01)
    x[C.AKT_activation_pIGF1R] = -0.1641    #6.85*10**(−01)
    x[C.AKT_activatio_pMetEGFR] = -5.0000   #1.00*10**(−05)
    x[C.AKT_activation_pMetErbB3] = -1.4321 #3.70*10**(−02)
    x[C.AKT_activation_pMetd] = -0.0479     #8.96*10**(−01)
    x[C.AKT_internIGF1R_effect] = -4.9907    #1.02*10**(−05)
    x[C.EGFR_BTC_binding] = -4.6498          #2.24*10**(−05)
    x[C.EGFR_BTC_dimerize] = 3.0000          #1.00*10**(+03)
    x[C.EGFR_ErbB2_BTC_dimerize] = 0.2072   #1.61*10**(+00)
    x[C.EGFR_ErbB2_basal_act] = -5.0000     #1.00*10**(−05)
    x[C.EGFR_ErbB2_dimerize] = -1.8041      #1.57*10**(−02)
    x[C.EGFR_ErbB3_BTC_dimerize] = -1.4539  #3.52*10**(−02)
    x[C.EGFR_ErbB3_basal_act] = -3.1243     #7.51*10**(−04)
    x[C.EGFR_ErbB3_dimerize] = -2.9275      #1.18*10**(−03)
    x[C.EGFR_ErbB3_dimerize_noHRG] = -5.0000 #1.00*10**(−05)
    x[C.EGFR_basal_activation] = -5.0000    #1.00*10**(−05)
    x[C.EGFR_basal_recycle] = +5.7152       #5.19*10**(+05)
    x[C.EGFR_dimerize] = -1.2005            #6.30*10**(−02)
    x[C.EGFR_lig_binding] = -4.7227         #1.89*10**(−05)
    x[C.EGF_kD] = +0.0000                    #1.00*10**(+00)
    x[C.ERK_phosphorylation_pMEK] = -3.5879 #2.58*10**(−04)
    x[C.ErbB2_ErbB3_dimerize] =  -1.5192    #3.03*10**(−02)
    x[C.ErbB2_dimerize] = -2.1680           #6.79*10**(−03)
    x[C.ErbB2_recycle] = -2.1824            #6.57*10**(−03)
    x[C.ErbB3_ErbB2_basal_act] = -5.0000    #1.00*10**(−05)
    x[C.ErbB3_basal_activation] = -1.5331   #2.93*10**(−02)
    x[C.ErbB3_basal_recycle] = -0.1327      #7.37*10**(−01)
    x[C.ErbB3_dimerize] = -1.3508           #4.46*10**(−02)
    x[C.ErbB3_lig_binding] = -4.2430        #5.71*10**(−05)
    x[C.HGF_kD] = -0.5229                   #3.00*10**(−01)
    x[C.HRG_kD] = -1.3010                   #5.00*10**(−02)
    x[C.IGF1R_basal_activation] = -2.9317   #1.17*10**(−03)
    x[C.IGF1R_basal_recycle] = +3.0000      #1.00*10**(+03)
    x[C.IGF1R_dimerize] = +1.2337           #1.71*10**(+01)
    x[C.IGF1R_lig_binding] = -2.8155        #1.53*10**(−03)
    x[C.IGF1_kD] = -0.5229                  #3.00*10**(−01)
    x[C.MEK_internIGF1R_effect] = -5.0000   #1.00*10**(−05)
    x[C.MEK_phosphorylation_pEGFR] = -5.0000    #1.00*10**(−05)
    x[C.MEK_phosphorylation_pErbB12] = -0.5625  #2.74*10**(−01)
    x[C.MEK_phosphorylation_pErbB13] = -5.0000  #1.00*10**(−05)
    x[C.MEK_phosphorylation_pErbB32] = -1.2178  #6.06*10**(−02)
    x[C.MEK_phosphorylation_pIGF1R] = -1.5236   #2.99*10**(−02)
    x[C.MEK_phosphorylation_pMetEGFR] = -5.0000 #1.00*10**(−05)
    x[C.MEK_phosphorylation_pMetErbB3] = -1.4164#3.83*10**(−02)
    x[C.MEK_phosphorylation_pMetd] = +0.2741    #1.88*10**(+00)
    x[C.Met_EGFR_BTC_dimerize] = -1.9539        #1.11*10**(−02)
    x[C.Met_EGFR_basal_act] = -4.7581           #1.75*10**(−05)
    x[C.Met_EGFR_dimerize] = -3.2898            #5.13*10**(−04)
    x[C.Met_ErbB3_basal_act] = +0.5177      #3.29*10**(+00)
    x[C.Met_ErbB3_dimerize] = -1.4308       #3.71*10**(−02)
    x[C.Met_basal_act] = -5.0000            #1.00*10**(−05)
    x[C.Met_dimerize] = -2.0378             #9.17*10**(−03)
    x[C.Met_lig_ErbB3_dimerize] = +2.7536   #5.67*10**(+02)
    x[C.Met_lig_binding] = -2.3452          #4.52*10**(−03)
    x[C.Met_recycle] = -0.2656              #5.42*10**(−01)
    x[C.S6K1_phosphorylation_pAKT] = -0.6024 #2.50*10**(−01)
    x[C.S6K1_phosphorylation_pERK] = -4.9718 #1.07*10**(−05)
    x[C.S6_phosphorylation_pERK] = -4.9987  #1.00*10**(−05)
    x[C.S6_phosphorylation_pS6K1] = -2.0671 #8.57*10**(−03)
    x[C.feedback_pAKT] = -4.9754            #1.06*10**(−05)
    x[C.feedback_pERK] = +3.0000            #1.00*10**(+03)
    x[C.feedback_pERK_on_AKT] = -4.9989     #1.00*10**(−05)
    x[C.feedback_pS6K1] = -4.6796           #2.09*10**(−05)
    x[C.init_AKT] = +0.4327                 #2.71*10**(+00)
    x[C.init_EGFR] = +1.2519                #1.79*10**(+01)
    x[C.init_EGFR_BTC] = +0.0000            #0.00*10**(+00)
    x[C.init_EGFR_EGF] = +0.0000            #0.00*10**(+00)
    x[C.init_ErbB2] = +0.7560               #5.70*10**(+00)
    x[C.init_ErbB3] = +0.3944               #2.48*10**(+00)
    x[C.init_ErbB3_HRG] = +0.0000           #0.00*10**(+00)
    x[C.init_IGF1R] = +0.6753               #4.73*10**(+00)
    x[C.init_IGF1R_IGF1] = +0.0000          #0.00*10**(+00)
    x[C.init_MEK] = +0.6275                 #4.24*10**(+00)
    x[C.init_Met] = +0.8978                 #7.90*10**(+00)
    x[C.init_Met_HGF] = +0.0000             #0.00*10**(+00)
    x[C.init_RTKph] = -0.2084               #6.19*10**(−01)
    x[C.init_S6] = +2.1629                  #1.46*10**(+02)
    x[C.init_pERK] = -0.4759                #3.34*10**(−01)
    x[C.init_pS6K1] = -2.9054               #1.24*10**(−03)

    #251-326
    x[C.pAKT_deactivation] = -0.5186        #3.03*10**(−01)
    x[C.pEGFR_degradation] = -5.0000        #1.00*10**(−05)
    x[C.pEGFR_internalize] = +0.7946        #6.23*10**(+00)
    x[C.pEGFR_phosphatase_binding] = +2.2663#1.85*10**(+02)
    x[C.pEGFRi_dephosph] = +1.3369          #2.17*10**(+01)
    x[C.pERK_dephosphorylation] = -0.2655   #5.43*10**(−01)
    x[C.pErbB12_degradation] = -0.7403      #1.82*10**(−01)
    x[C.pErbB12_internalize] = +0.2777      #1.90*10**(+00)
    x[C.pErbB12i_dephosph] = +3.0000        #1.00*10**(+03)
    x[C.pErbB12i_phosphatase] = -0.0730     #8.45*10**(−01)
    x[C.pErbB13_degradation] = +1.6965      #4.97*10**(+01)
    x[C.pErbB13_internalize] = +3.0000      #1.00*10**(+03)
    x[C.pErbB13i_dephosph] = +1.7665        #5.84*10**(+01)
    x[C.pErbB13i_phosphatase] = -4.7972     #1.60*10**(−05)
    x[C.pErbB2_degradation] = +2.5034       #3.19*10**(+02)
    x[C.pErbB2_internalize] = +3.0000       #1.00*10**(+03)
    x[C.pErbB2i_dephosph] = +0.9161         #8.24*10**(+00)
    x[C.pErbB2i_phosphatase] = +3.0000      #1.00*10**(+03)
    x[C.pErbB32_degradation]  = -0.2414     #5.74*10**(−01)
    x[C.pErbB32_internalize] = +0.0363      #1.09*10**(+00)
    x[C.pErbB32i_dephosph] = -1.9237        #1.19*10**(−02)
    x[C.pErbB32i_phosphatase] = -1.3117     #4.88*10**(−02)
    x[C.pErbB3_degradation] = -0.0417       #9.08*10**(−01)
    x[C.pErbB3_internalize] = +3.0000       #1.00*10**(+03)
    x[C.pErbB3i_dephosph] = +1.2126         #1.63*10**(+01)
    x[C.pErbB3i_phosphatase] = +1.9148      #8.22*10**(+01)
    x[C.pIGF1R_degradation] = -5.0000       #1.00*10**(−05)
    x[C.pIGF1R_internalize] = +3.0000       #1.00*10**(+03)
    x[C.pIGF1Ri_dephosph] = +2.5207         #3.32*10**(+02)
    x[C.pIGF1Ri_phosphatase] = +3.0000      #1.00*10**(+03)
    x[C.pMEK_dephosphorylation] = -0.4839   #3.28*10**(−01)
    x[C.pMetEGFR_degradation] = +0.1539     #1.43*10**(+00)
    x[C.pMetEGFR_internalize] = +0.1229     #1.33*10**(+00)
    x[C.pMetEGFRi_dephosph] = -0.2978       #5.04*10**(−01)
    x[C.pMetEGFRi_phosphatase] = -4.4470    #3.57*10**(−05)
    x[C.pMetErbB3_degradation] = +2.9929    #9.84*10**(+02)
    x[C.pMetErbB3_internalize] = +3.0000    #1.00*10**(+03)
    x[C.pMetErbB3i_dephosph] = +3.0000      #1.00*10**(+03)
    x[C.pMetErbB3i_phosphatase] = +3.0000   #1.00*10**(+03)
    x[C.pMet_degradation] = +0.2802         #1.91*10**(+00)
    x[C.pMet_internalize] = +3.0000         #1.00*10**(+03)
    x[C.pMeti_dephosph] = +1.1134           #1.30*10**(+01)
    x[C.pMeti_phosphatase] = +3.0000        #1.00*10**(+03)
    x[C.pS6K1_dephosphorylation] = -1.9859  #1.03*10**(−02)
    x[C.pS6_dephosphorylation] = -0.9252    #1.19*10**(−01)
    x[C.relto_A431_init_EGFR] = +0.9838     #9.63*10**(+00)
    x[C.relto_A431_init_ErbB2] = +0.0013    #1.00*10**(+00)
    x[C.relto_A431_init_ErbB3] = -0.0403    #9.11*10**(−01)
    x[C.relto_A431_init_IGF1R] = +0.4236    #2.65*10**(+00)
    x[C.relto_A431_init_Met] = -0.6500      #2.24*10**(−01)
    x[C.relto_ACHN_init_EGFR] = +0.7324     #5.40*10**(+00)
    x[C.relto_ACHN_init_ErbB2] = -0.1347    #7.33*10**(−01)
    x[C.relto_ACHN_init_ErbB3] = -0.3509    #4.46*10**(−01)
    x[C.relto_ACHN_init_IGF1R] = -0.1530    #7.03*10**(−01)
    x[C.relto_ACHN_init_Met] = +0.1113      #1.29*10**(+00)
    x[C.relto_ADRr_init_EGFR] = +0.0119     #1.03*10**(+00)
    x[C.relto_ADRr_init_ErbB2] = +0.0311    #1.07*10**(+00)
    x[C.relto_ADRr_init_ErbB3] = -0.2414    #5.74*10**(−01)
    x[C.relto_ADRr_init_IGF1R] = +0.3059    #2.02*10**(+00)
    x[C.relto_ADRr_init_Met] = -0.3171      #4.82*10**(−01)
    x[C.relto_BT20_init_EGFR] = +1.0201     #1.05*10**(+01)
    x[C.relto_BT20_init_ErbB2] = +0.1754    #1.50*10**(+00)
    x[C.relto_BT20_init_ErbB3] = +0.2661    #1.85*10**(+00)
    x[C.relto_BT20_init_IGF1R] = +0.5507    #3.55*10**(+00)
    x[C.relto_BT20_init_Met] = -0.4177      #3.82*10**(−01)
    x[C.relto_IGROV_init_EGFR] = -0.0264    #9.41*10**(−01)
    x[C.relto_IGROV_init_ErbB2] = +0.6286   #4.25*10**(+00)
    x[C.relto_IGROV_init_ErbB3] = -0.2399   #5.76*10**(−01)
    x[C.relto_IGROV_init_IGF1R] = -0.7987   #1.59*10**(−01)
    x[C.relto_IGROV_init_Met] = -0.1926     #6.42*10**(−01)
    x[C.relto_init_EGFR] = +0.2627          #1.83*10**(+00)
    x[C.relto_init_ErbB2] = -0.2169         #6.07*10**(−01)
    x[C.relto_init_ErbB3] = +0.0165         #1.04*10**(+00)
    x[C.relto_init_IGF1R] = +0.1713         #1.48*10**(+00)
    x[C.relto_init_Met] = +0.2563           #1.80*10**(+00)
    x[C.scale_Ligand] = +4.5776             #3.78*10**(+04)   

    return x


def initial_values():
    y0 = [0]*V.NUM

    return y0
''''
# call class
if __name__ == '__main__':
    D = DifferentialEquation(3) #create new instance
    print(D.pertubation)
'''

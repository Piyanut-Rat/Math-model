import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from name2idx import parameters as C
from name2idx import species as V
from set_model import initial_values, param_values

class SearchParam(object):

    # parameters with index
    idx_params = [
        C.EGFR_prod,
        C.ErbB2_prod,
        C.ErbB3_prod,
        C.IGF1R_prod,
        C.Met_prod,
        C.EGER_lig_binding,
        C.EGF_kD,
        C.EGFR_BTC_binding,
        C.ErbB3_lig_binding,
        C.IGF1R_lig_binding,
        C.Met_lig_binding,
        C.EGFR_dimerize,
        C.EGER_BTC_dimerize,
        C.ErbB2_dimerize,
        C.ErbB3_dimerize,
        C.IGF1R_dimerize,
        C.EGFR_ErbB2_dimerize,
        C.EGFR_ErbB2_BTC_dimerize,
        C.EGFR_ErbB3_dimerize,
        C.EGFR_ErbB3_BTC_dimerize,
        C.EGFR_ErbB3_dimerize_noHRG,
        C.ErbB2_ErbB3_dimerize,
        C.Met_dimerize,
        C.Met_ErbB3_dimerize,
        C.Met_lig_ErbB3_dimerize,
        C.Met_EGFR_dimerize,
        C.Met_EGFR_BTC_dimerize,
        C.EGFR_basal_activation,
        C.ErbB3_basal_activation,
        C.IGF1R_basal_activation,
        C.EGFR_ErbB2_basal_act,
        C.EGFR_ErbB3_basal_act,
        C.ErbB3_ErbB2_basal_act,
        C.Met_ErbB3_basal_act,
        C.Met_basal_act,
        C.Met_EGFR_basal_act,
        C.pEGFR_degradation,
        C.pEGFR_phosphatase_binding,
        C.pEGFRi_dephosph,
        C.EGFR_basal_recycle,
        C.pErbB2_internalize,
        C.pErbB2i_phosphatase,
        C.pErbB2i_dephosph,
        C.pErbB2_degradation,
        C.ErbB2_recycle,
        C.pErbB3_internalize,
        C.pErbB3_degradation,
        C.pErbN3i_phosphatase,
        C.pErbB3i_dephosph,
        C.ErbB3_basal_recycle,
        C.pIGF1R_internalize,
        C.pIGF1R_degradation,
        C.pIGF1Ri_phosphatase,
        C.pIGF1Ri_dephosph,
        C.IGF1R_basal_recycle,
        C.pErbB12_internalize,
        C.pErbB12_degradation,
        C.pErbB12i_phosphatase,
        C.pErbB12i_dephosph,
        C.pErbB32_internalize,
        C.pErbB32_degradation,
        C.pErbB32i_phasphatase,
        C.pErbB32i_dephosph,
        C.pErbB13_internalize,
        C.pErbB13_degradation,
        C.pErbB12i_phosphatase,
        C.pErbB13i_dephosph,
        C.pMet_internalize,
        C.pMet_degradation,
        C.pMet_phosphatase,
        C.pMet_dephosph,
        C.Met_recycle,
        C.pMetErbB3_internalize,
        C.pMetErbB3_degradation,
        C.pMetErbB3i_phasphatase,
        C.pMetErbB3i_dephosph,
        C.pMetEGFR_internalize,
        C.pMetEGFR_degradation,
        C.pMetEGFRi_phosphatase,
        C.pMetEGFRi_dephosph,
        C.MEK_phosphorylation_pErbB12,
        C.feedback_pAKT,
        C.init_pAKT,
        C.feedback_pERK,
        C.MEK_phosphorylation_pErbB13,
        C.MEK_phosphorylation_pMetd,
        C.MEK_phosphorylation_pMetEGFR,
        C.MEK_phosphorylation_pIGF1R,
        C.MEK_phosphorylation_pMetErbB3,
        C.MEK_internIGF1R_effect,
        C.MEK_phosphorylation_pIGF1R,
        C.pMEK_dephosphorylation,
        C.ERK_phosphorylation_pMEK,
        C.pERK_dephosphorylation,
        C.AKT_activation_pEGFR,
        C.feedback_pERK_on_AKT,
        C.AKT_activation_pErbB12,
        C.feedback_pS6K1,
        C.AKT_activation_pErbB13,
        C.feedback_pERK_on_AKT,
        C.AKT_activation_pErbB32,
        C.AKT_activation_pMetEGFR,
        C.AKT_activation_pMetd,
        C.AKT_activation_pIGF1R,
        C.AKT_activation_pMetErbB3,
        C.AKT_internIGF1R_effect,
        C.pAKT_deactivation,
        C.S6K1_phosphorylation_pAKT,
        C.S6K1_phosphorylation_pERK,
        C.pS6K1_dephosphorylation,
        C.S6_phosphorylation_pS6K1,
        C.S6_phosphorylation_pERK,
        C.pS6_dephosphorylation,    
    ]
    # initial values
    idx_initials = [
        # V.(specie)
    ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = initialize_search_param(
            parameters = C.NAMES,
            species = V.NAMES,
            param_values = x,
            initial_values = y0,
            estimated_params=self.idx_params,
            estimated_initials = self.idx_initials,
        )

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i]* 0.1     #lower bound
            search_rgn[1, j] = search_param[i]* 10.0    #uppwer bound

        # default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)]* 0.5     #lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)]*2.0      #upper bound

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] =[lower_bound, upper_bound]

    
        ## log10(100000) = 5; non-log of 5 is 100000 (1e+5) 

        search_rgn[:, C.AKT_activation_pEGFR] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activation_pErbB12] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activation_pErbB13] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activation_pErbB32] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activation_pIGF1R] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activatio_pMetEGFR] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activation_pMetErbB3] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_activation_pMetd] = [1e-5, 1e+3]
        search_rgn[:, C.AKT_internIGF1R_effect] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_BTC_binding] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_BTC_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_ErbB2_BTC_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_ErbB2_basal_act] = [1e-5, 1e+6]

        search_rgn[:, C.EGFR_ErbB2_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_ErbB3_BTC_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_ErbB3_basal_act] = [1e-5, 1e+6]

        search_rgn[:, C.EGFR_ErbB3_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_ErbB3_dimerize_noHRG] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_basal_activation] = [1e-5, 1e+6]
        search_rgn[:, C.EGFR_basal_recycle] = [1e-5, 1e+6]

        search_rgn[:, C.EGFR_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.EGFR_lig_binding] = [1e-5, 1e+3]
        search_rgn[:, C.EGF_kD] = [1e-5, 1e+3]
        search_rgn[:, C.ERK_phosphorylation_pMEK] = [1e-5, 1e+3]
        search_rgn[:, C.ErbB2_ErbB3_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.ErbB2_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.ErbB2_recycle] = [1e-5, 1e+3]
        search_rgn[:, C.ErbB3_ErbB2_basal_act] = [1e-5, 1e+6]
        search_rgn[:, C.ErbB3_basal_activation] = [1e-5, 1e+6]
        search_rgn[:, C.ErbB3_basal_recycle] = [1e-5, 1e+6]

        search_rgn[:, C.ErbB3_dimerizee] = [1e-5, 1e+3]
        search_rgn[:, C.ErbB3_lig_binding] = [1e-5, 1e+3]
        search_rgn[:, C.HGF_kD] = [1e-5, 1e+3]
        search_rgn[:, C.C.HRG_kD] = [1e-5, 1e+3]
        search_rgn[:, C.IGF1R_basal_activation] = [1e-5, 1e+3]
        search_rgn[:, C.IGF1R_basal_recycle] = [1e-5, 1e+3]
        search_rgn[:, C.IGF1R_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.IGF1R_lig_binding] = [1e-5, 1e+3]
        search_rgn[:, C.IGF1_kD] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_internIGF1R_effect] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pEGFR] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pErbB12] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pErbB13] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pErbB32] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pIGF1R] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pMetEGFR] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pMetErbB3] = [1e-5, 1e+3]
        search_rgn[:, C.MEK_phosphorylation_pMetd] = [1e-5, 1e+3]
        search_rgn[:, C.Met_EGFR_BTC_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.Met_EGFR_basal_act] = [1e-5, 1e+6]

        search_rgn[:, C.Met_EGFR_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.Met_ErbB3_basal_act] = [1e-5, 1e+6]

        search_rgn[:, C.Met_ErbB3_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.Met_basal_act] = [1e-5, 1e+6]

        search_rgn[:, C.Met_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.Met_lig_ErbB3_dimerize] = [1e-5, 1e+3]
        search_rgn[:, C.Met_lig_binding] = [1e-5, 1e+3]
        search_rgn[:, C.Met_recycle] = [1e-5, 1e+3]
        search_rgn[:, C.S6K1_phosphorylation_pAKT] = [1e-5, 1e+3]
        search_rgn[:, C.S6K1_phosphorylation_pERK] = [1e-5, 1e+3]
        search_rgn[:, C.S6_phosphorylation_pERK] = [1e-5, 1e+3]
        search_rgn[:, C.S6_phosphorylation_pS6K1] = [1e-5, 1e+3]
        search_rgn[:, C.feedback_pAKT] = [1e-5, 1e+3]
        search_rgn[:, C.feedback_pERK] = [1e-5, 1e+3]
        search_rgn[:, C.feedback_pERK_on_AKT] = [1e-5, 1e+3]
        search_rgn[:, C.feedback_pS6K1] = [1e-5, 1e+3]
        search_rgn[:, C.init_AKT] = [1e-5, 1e+3]
        search_rgn[:, C.init_EGFR] = [1e-5, 1e+3]
        search_rgn[:, C.init_EGFR_BTC] = [1e-5, 1e+3]
        search_rgn[:, C.init_EGFR_EGF] = [1e-5, 1e+3]
        search_rgn[:, C.init_ErbB2] = [1e-5, 1e+3]
        search_rgn[:, C.init_ErbB3] = [1e-5, 1e+3]
        search_rgn[:, C.init_ErbB3_HRG] = [1e-5, 1e+3]
        search_rgn[:, C.init_IGF1R] = [1e-5, 1e+3]
        search_rgn[:, C.init_IGF1R_IGF1] = [1e-5, 1e+3]
        search_rgn[:, C.init_MEK] = [1e-5, 1e+3]
        search_rgn[:, C.init_Met] = [1e-5, 1e+3]
        search_rgn[:, C.init_Met_HGF] = [1e-5, 1e+3]
        search_rgn[:, C.init_RTKph] = [1e-5, 1e+3]
        search_rgn[:, C.init_S6] = [1e-5, 1e+3]
        search_rgn[:, C.init_pERK] = [1e-5, 1e+3]
        search_rgn[:, C.init_pS6K1] = [1e-5, 1e+3]

        #251-326
        search_rgn[:, C.pAKT_deactivation] = [1e-5, 1e+3]
        search_rgn[:, C.pEGFR_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pEGFR_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pEGFR_phosphatase_binding] = [1e-5, 1e+3]
        search_rgn[:, C.pEGFRi_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pERK_dephosphorylation] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB12_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB12_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB12i_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB12i_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB13_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB13_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB13i_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB13i_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB2_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB2_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB2i_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB2i_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB32_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB32_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB32i_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB32i_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB3_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB3_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB3i_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pErbB3i_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pIGF1R_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pIGF1R_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pIGF1Ri_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pIGF1Ri_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pMEK_dephosphorylation] = [1e-5, 1e+3]
        search_rgn[:, C.pMetEGFR_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pMetEGFR_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pMetEGFRi_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pMetEGFRi_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pMetErbB3_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pMetErbB3_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pMetErbB3i_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pMetErbB3i_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pMet_degradation] = [1e-5, 1e+3]
        search_rgn[:, C.pMet_internalize] = [1e-5, 1e+3]
        search_rgn[:, C.pMeti_dephosph] = [1e-5, 1e+3]
        search_rgn[:, C.pMeti_phosphatase] = [1e-5, 1e+3]
        search_rgn[:, C.pS6K1_dephosphorylation] = [1e-5, 1e+3]
        search_rgn[:, C.pS6_dephosphorylation] = [1e-5, 1e+3]

        search_rgn[:, C.relto_A431_init_EGFR] = [1e-4, 1e+4]
        search_rgn[:, C.relto_A431_init_ErbB2] = [1e-4, 1e+4]
        search_rgn[:, C.relto_A431_init_ErbB3] = [1e-4, 1e+4]
        search_rgn[:, C.relto_A431_init_IGF1R] = [1e-4, 1e+4]
        search_rgn[:, C.C.relto_A431_init_Met] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ACHN_init_EGFR] = [1e-4, 1e+4]     
        search_rgn[:, C.relto_ACHN_init_ErbB2] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ACHN_init_ErbB3] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ACHN_init_IGF1R] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ACHN_init_Met] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ADRr_init_EGFR] = [1e-4, 1e+4] 
        search_rgn[:, C.relto_ADRr_init_ErbB2] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ADRr_init_ErbB3] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ADRr_init_IGF1R] = [1e-4, 1e+4]
        search_rgn[:, C.relto_ADRr_init_Met] = [1e-4, 1e+4]
        search_rgn[:, C.relto_BT20_init_EGFR] = [1e-4, 1e+4] 
        search_rgn[:, C.relto_BT20_init_ErbB2] = [1e-4, 1e+4]
        search_rgn[:, C.relto_BT20_init_ErbB3] = [1e-4, 1e+4]
        search_rgn[:, C.relto_BT20_init_IGF1R] = [1e-4, 1e+4]
        search_rgn[:, C.relto_BT20_init_Met] = [1e-4, 1e+4]
        search_rgn[:, C.relto_IGROV_init_EGFR] = [1e-4, 1e+4]
        search_rgn[:, C.relto_IGROV_init_ErbB2] = [1e-4, 1e+4]
        search_rgn[:, C.relto_IGROV_init_ErbB3] = [1e-4, 1e+4]
        search_rgn[:, C.relto_IGROV_init_IGF1R] = [1e-4, 1e+4]
        search_rgn[:, C.relto_IGROV_init_Met] = [1e-4, 1e+4]
        search_rgn[:, C.relto_init_EGFR] = [1e-4, 1e+4]
        search_rgn[:, C.relto_init_ErbB2] = [1e-4, 1e+4]
        search_rgn[:, C.relto_init_ErbB3] = [1e-4, 1e+4]
        search_rgn[:, C.relto_init_IGF1R] = [1e-4, 1e+4]
        search_rgn[:, C.relto_init_Met] = [1e-4, 1e+4]

        search_rgn[:, C.scale_Ligand] = [1e+1, 1e+6]

        search_rgn = convert_scale(
            region = search_rgn,
            parameters = C.NAMES,
            species = V.NAMES,
            estimated_params = self.idx_params,
            estimated_initials = self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10 ** (indiv_gene * (search_rgn[1, :] - search_rgn[0, :]) + search_rgn[0, :])

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (np.log10(indiv) - search_rgn[0, :]) / (search_rgn[1, :] - search_rgn[0, :])

        return indiv_gene


        
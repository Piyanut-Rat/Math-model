NAMES = [
    'EGFR_prod',
    'ErbB2_prod',
    'ErbB3_prod',
    'IGF1R_prod',
    'Met_prod',
    'EGER_lig_binding',
    'EGF_kD',
    'EGFR_BTC_binding',
    'ErbB3_lig_binding',
    'IGF1R_lig_binding',
    'Met_lig_binding',
    'EGFR_dimerize',
    'EGER_BTC_dimerize',
    'ErbB2_dimerize',
    'ErbB3_dimerize',
    'IGF1R_dimerize',
    'EGFR_ErbB2_dimerize',
    'EGFR_ErbB2_BTC_dimerize',
    'EGFR_ErbB3_dimerize',
    'EGFR_ErbB3_BTC_dimerize',
    'EGFR_ErbB3_dimerize_noHRG',
    'ErbB2_ErbB3_dimerize',
    'Met_dimerize',
    'Met_ErbB3_dimerize',
    'Met_lig_ErbB3_dimerize',
    'Met_EGFR_dimerize',
    'Met_EGFR_BTC_dimerize',
    'EGFR_basal_activation',
    'ErbB3_basal_activation',
    'IGF1R_basal_activation',
    'EGFR_ErbB2_basal_act',
    'EGFR_ErbB3_basal_act',
    'ErbB3_ErbB2_basal_act',
    'Met_ErbB3_basal_act',
    'Met_basal_act',
    'Met_EGFR_basal_act',
    'pEGFR_degradation',
    'pEGFR_phosphatase_binding',
    'pEGFRi_dephosph',
    'EGFR_basal_recycle',
    'pErbB2_internalize',
    'pErbB2i_phosphatase',
    'pErbB2i_dephosph',
    'pErbB2_degradation',
    'ErbB2_recycle',
    'pErbB3_internalize',
    'pErbB3_degradation',
    'pErbN3i_phosphatase',
    'pErbB3i_dephosph',
    'ErbB3_basal_recycle',
    ''


]
for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)
NAMES = [
    'dose_EGF',
    'dose_HGF',
    'RTKph',
    'dose_IGF1',
    'dose_HRG',
    'dose_BTC',
    'EGER',
    'EGFR_EGF',
    'EGFR_BTC',
    'pEGFRd',


    'pEGFRi',
    'pEGFRi_ph',
    'EGFRi',
    'ErbB2',
    'pErbB2',
    'pErbB2i',
    'ErbB2i',
    'pErbB2i_ph',
    'pErbB12',
    'pErbB12i',

    'pErbB12i_ph',
    'ErbB3',
    'ErbB3_HRG',
    'pErbB3d',
    'ppErbB3i',
    'pErbB3i_ph',
    'ErbB3i',
    'pErbB13',
    'pErbB13i',
    'pErbB13i_ph',

    'pErbB32',
    'pErbB32i',
    'pErbB32i_ph',
    'IGF1R',
    'IGF1R_IGF1',
    'pIGF1Rd',
    'pIGF1Ri',
    'pIGF1Ri_ph',
    'IGF1Ri',
    'Met',


    'Met_HGF',
    'pMetd',
    'pMeti',
    'pMeti_ph',
    'Meti',
    'pMetErbB3',
    'pMetErbB3i',
    'pMetErbB3i_ph',
    'pMetEGFR',
    'pMetEGFRi',


    'pMetEGRFi_ph',
    'MEK',
    'pMEK',
    'ERK',
    'pERK',
    'AKT',
    'pAKT',
    'S6K1',
    'pS6K1',
    'S6',
    'pS6',

    'EGFR',
    'EGFR_ErbB3_dimerize_noHRG',
    'pErbB3i',
    'pMetEGFRi_ph',
]
for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)
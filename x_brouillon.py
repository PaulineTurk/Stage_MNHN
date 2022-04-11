import FUNCTION.character as ch
import numpy as np
import pandas as pd


def sumLine(cond_proba, list_AA, aa_k, aa_c):
    sum_line = 0
    for aa_p in list_AA:
        sum_line += cond_proba[aa_k][aa_p][aa_c]
    return sum_line

def predictorIdentityNeighbor():
    predictor_name = "Identity Predictor Neighbor"
    list_AA = ch.characterList()
    cond_proba_id = {}
    for elem_k in list_AA:
        cond_proba_id[elem_k] = {}
        for elem_p in list_AA:
            cond_proba_id[elem_k][elem_p] = {}
            for elem_c in list_AA:
                if elem_k == elem_p:
                    cond_proba_id[elem_k][elem_p][elem_c] = 1
                else:
                    cond_proba_id[elem_k][elem_p][elem_c] = 0

    df_cond_proba_id = np.transpose(pd.DataFrame.from_dict(cond_proba_id))
    print("{}:\n".format(predictor_name), df_cond_proba_id)
    for aa_k in list_AA:
        for aa_c in list_AA:
            sum_line = sumLine(cond_proba_id, list_AA, aa_k, aa_c)   
            print("sum line (aa_k: {}, aa_c: {}): {}".format(aa_k, aa_c, sum_line)) 
    
    unit_brier_id = unitBrierNeighbor(cond_proba_id)
    return predictor_name, cond_proba_id, unit_brier_id


def unitBrierNeighbor(cond_proba):  
    list_AA = ch.characterList()
    unit_Brier = {}
    for aa_k in list_AA:
        unit_Brier[aa_k] = {}
        for aa_p in list_AA:
            unit_Brier[aa_k][aa_p] = {}
            for aa_c in list_AA:
                unit = 0
                for j in list_AA: 
                    unit += (cond_proba[aa_k][j][aa_c] - int(aa_p == j))**2    # à ne regarder que la realisation de j (tjs)
                unit_Brier[aa_k][aa_p][aa_c] = unit
    return unit_Brier

predictorIdentityNeighbor()
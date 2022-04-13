import character as ch
import numpy as np
import pandas as pd


def sumLine(cond_proba, list_AA, aa_k, aa_c):
    sum_line = 0
    for aa_p in list_AA:
        sum_line += cond_proba[aa_k][aa_p][aa_c]
    return sum_line

def predictorStationaryNeighbor(cond_proba):                      # A REPRENDRE D ICI
    predictor_name = "Stationary Predictor Neighbor"
    list_AA = ch.characterList()

    # freq of each AA knowing aa_p and the contexte     (maybe to calculate outside this function ...)
    freq_aa_p = {}
    for aa_p in list_AA:
        freq_aa_p[aa_p] = 0
        for aa_k in list_AA:
            for aa_c in list_AA: 
                freq_aa_p[aa_p] += cond_proba[aa_k][aa_p][aa_c]

    cond_proba_stationary = {}
    for elem_k in list_AA:
        cond_proba_stationary[elem_k] = {}
        for elem_p in list_AA:
            cond_proba_stationary[elem_k][elem_p] = {}
            for elem_c in list_AA:
                cond_proba_stationary[elem_k][elem_p][elem_c] = freq_aa_p[aa_p]
    df_cond_proba_stationary = np.transpose(pd.DataFrame.from_dict(cond_proba_stationary))
    print("{}:\n".format(predictor_name), df_cond_proba_stationary)
    for aa_k in list_AA:
        for aa_c in list_AA:
            sum_line = sumLine(cond_proba_stationary, list_AA, aa_k, aa_c)   
            print("sum line (aa_k: {}, aa_c: {}): {}".format(aa_k, aa_c, sum_line)) 
    
    unit_Brier_stationnaire = unitBrierNeighbor(cond_proba_stationary)
    return predictor_name, cond_proba_stationary, unit_Brier_stationnaire


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





#cond_proba = np.load("/Users/pauline/Desktop/data/NeighborRes/NeighborRes_0.05/proba_cond_(-1 , k).npy", allow_pickle='TRUE').item()
#predictorStationaryNeighbor(cond_proba)
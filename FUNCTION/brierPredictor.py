import FUNCTION.character as ch
import pandas as pd
from random import choice
#from numpy import transpose 
import numpy as np


def predictor01(liste_seq, num_accession, nom_dir_pid, Brier_count_global, count_global, pid_inf = 62):
    pid_couple = np.load(nom_dir_pid + "/" + num_accession + ".pId.npy", allow_pickle='TRUE').item()
    list_AA = ch.characterList()
    for name_1, seq_1 in liste_seq:
        for name_2 ,seq_2 in liste_seq:
            if name_1 != name_2:      # ATTANTION, ATTENTION, ATTENTION, à vérifier qu'on ne les compte pas
                if pid_couple[name_1][name_2] >= pid_inf:
                    for (aa_1, aa_2) in zip(seq_1, seq_2):
                        if aa_1 in list_AA and aa_2 in list_AA:
                            count_global += 1
                            new_listeAA = ch.characterList()
                            new_listeAA.remove(aa_2)
                            AApredit = choice(new_listeAA)
                            for j in list_AA:
                                if aa_2 == j:
                                    Brier_count_global += (0 - 1)**2
                                else:
                                    if j == AApredit:
                                        Brier_count_global += (1 - 0)**2
                                    else:
                                        Brier_count_global += (0 - 0)**2
    return Brier_count_global, count_global


def predictorPerfect(liste_seq, num_accession, nom_dir_pid, Brier_count_global, count_global, pid_inf = 62):
    pid_couple = np.load(nom_dir_pid + "/" + num_accession + ".pId.npy", allow_pickle='TRUE').item()
    list_AA = ch.characterList()
    for name_1, seq_1 in liste_seq:
        for name_2 ,seq_2 in liste_seq:
            if name_1 != name_2:      # ATTANTION, ATTENTION, ATTENTION, à vérifier qu'on ne les compte pas
                if pid_couple[name_1][name_2] >= pid_inf:
                    for (aa_1, aa_2) in zip(seq_1, seq_2):
                        if aa_1 in list_AA and aa_2 in list_AA:
                            count_global += 1
                            for j in list_AA:
                                if aa_2 == j:
                                    Brier_count_global += (1 - 1)**2
                                else:
                                    Brier_count_global += (0 - 0)**2
    return Brier_count_global, count_global



def predictorBlosum(name_matrix_cond_proba):
    predictor_name = "Blosum Predictor"
    cond_proba_Blosum = np.load(name_matrix_cond_proba ,allow_pickle='TRUE').item()
        
    #df_cond_proba_Blosum = transpose(pd.DataFrame.from_dict(cond_proba_Blosum))
    #print("{}:\n".format(predictor_name), df_cond_proba_Blosum)
    #sum_ligne = df_cond_proba_Blosum.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)

    unit_Brier_Blosum = unitBrier(cond_proba_Blosum)
    return predictor_name, cond_proba_Blosum, unit_Brier_Blosum


def predicteurEquiprobable():
    predictor_name = "Equiprobable Predictor"
    liste_AA = ch.characterList()
    nbreAA = len(liste_AA)
    cond_proba_equiproba = {}
    for elem_l in liste_AA:
        cond_proba_equiproba[elem_l] = {}
        for elem_c in liste_AA:
            cond_proba_equiproba[elem_l][elem_c] = 1/nbreAA

    #df_cond_proba_equiproba = transpose(pd.DataFrame.from_dict(cond_proba_equiproba))
    #print("{}:\n".format(predictor_name), df_cond_proba_equiproba)
    #sum_ligne = df_cond_proba_equiproba.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)

    unit_Brier_equiproba = unitBrier(cond_proba_equiproba)
    return predictor_name, cond_proba_equiproba, unit_Brier_equiproba


def predictorStationary(freq_aa):
    predictor_name = "Stationary Predictor"
    liste_AA = ch.characterList()
    cond_proba_stationary = {}
    for elem_l in liste_AA:
        cond_proba_stationary [elem_l] = {}
        for elem_c in liste_AA:
            if freq_aa != 0:
                cond_proba_stationary [elem_l][elem_c] = freq_aa[elem_c]
            else:
                cond_proba_stationary [elem_l][elem_c] = 0

    #df_cond_proba_stationary = transpose(pd.DataFrame.from_dict(cond_proba_stationary))
    #print("{}:\n".format(predictor_name), df_cond_proba_stationary)
    #sum_ligne = df_cond_proba_stationary.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)
    
    unit_Brier_stationnaire = unitBrier(cond_proba_stationary)
    return predictor_name, cond_proba_stationary, unit_Brier_stationnaire


def predictorIdentity():
    predictor_name = "Identity Predictor"
    list_AA = ch.characterList()
    cond_proba_id = {}
    for elem_l in list_AA:
        cond_proba_id[elem_l] = {}
        for elem_c in list_AA:
            if elem_l == elem_c:
                cond_proba_id[elem_l][elem_c] = 1
            else:
                cond_proba_id[elem_l][elem_c] = 0

    #df_cond_proba_id = transpose(pd.DataFrame.from_dict(cond_proba_id))
    #print("{}:\n".format(predictor_name), df_cond_proba_id)
    #sum_ligne = df_cond_proba_id.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)
    
    unit_brier_id = unitBrier(cond_proba_id)
    return predictor_name, cond_proba_id, unit_brier_id


def unitBrier(cond_proba):  
    list_AA = ch.characterList()
    unit_Brier = {}
    for aa_1 in list_AA:
        unit_Brier[aa_1] = {}
        for aa_2 in list_AA:
            unit = 0
            for j in list_AA: 
                unit += (cond_proba[aa_1][j] - int(aa_2 == j))**2 
            unit_Brier[aa_1][aa_2] = unit
    return unit_Brier






def brierMatrix(predictor_name, unit_Brier, liste_seq, accession_num, dir_pid_name, Brier_count_global, count_global, pid_inf = 62):
    pid_couple = np.load(dir_pid_name + "/" + accession_num + ".pId.npy", allow_pickle='TRUE').item()
    list_AA =  ch.characterList()
    seq_nbre = len(liste_seq)

    if predictor_name in ["Blosum Predictor", "Equiprobable Predictor",  "Stationary Predictor", "Identity Predictor"]:  
        for i in range(seq_nbre):
            name_1, seq_1 = liste_seq[i]
            for j in range(i + 1, seq_nbre):
                name_2 ,seq_2 = liste_seq[j] 
                if pid_couple[name_1][name_2] >= pid_inf:
                    for (aa_1, aa_2) in zip(seq_1, seq_2):
                        if aa_1 in list_AA and aa_2 in list_AA:
                            Brier_count_global += unit_Brier[aa_1][aa_2]
                            Brier_count_global += unit_Brier[aa_2][aa_1]
                            count_global += 2   
    return Brier_count_global, count_global
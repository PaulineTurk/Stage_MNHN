from timer import Timer
from fastaReader import readFastaMul
import numpy as np
from pathlib import Path
import pandas as pd
import os, shutil
import loadVisu

# tripletCount(path_file_fasta, triplet_count, delay_num, kp_SeqChoice, list_residu)

def tripletCount(path_file_fasta, triplet_count, delay_num, kp_SeqChoice, list_residu):    
    liSeqAliFiltre = readFastaMul(path_file_fasta)   

    if liSeqAliFiltre:
        for name_k, seq_k in liSeqAliFiltre:
            len_seq = len(seq_k)
            for name_p, seq_p in liSeqAliFiltre:
                if name_k != name_p:
                    print(f"couple seq k,p: {name_k}, {name_p}")
                    if delay_num < 0:   # before
                        index_range = - delay_num, len_seq
                    elif delay_num > 0: # after
                        index_range = 0, len_seq - delay_num
                        
                    if kp_SeqChoice == "k":
                        seq_c = seq_k
                    elif kp_SeqChoice == "p":
                        seq_c = seq_p

                    for aa_index in range(index_range[0], index_range[1]):       
                        aa_k = seq_k[aa_index] 
                        aa_p = seq_p[aa_index] 
                        if all(x in list_residu for x in [aa_k, aa_p]):     
                            index_neighbor = aa_index + delay_num
                            aa_c = seq_c[index_neighbor] 
                            if aa_c in list_residu: 
                                print(f"triplet: {aa_k},{aa_p},{aa_c}")
                                triplet_count[aa_k][aa_p][aa_c] += 1

    return triplet_count








def initialisation(list_symbol):
    triplet_count = {}
    for aa_k in list_symbol:
        triplet_count[aa_k] = {}
        for aa_p in list_symbol:
            triplet_count[aa_k][aa_p] = {}
            for aa_c in list_symbol:
                triplet_count[aa_k][aa_p][aa_c]  = 1   # to avoid issues when the triplet is not in data_train 
    return triplet_count




def conditionalProba(list_residu, triplet_count):
    # pseudo_count idea removed because we have enough data
    intra_couple_count = {}
    for aa_k in list_residu:
        intra_couple_count[aa_k] = {}
        for aa_c in list_residu: 
            intra_couple_count[aa_k][aa_c] = 0
            for aa_p in list_residu:
                intra_couple_count[aa_k][aa_c] += triplet_count[aa_k][aa_p][aa_c]
                print(f"intra_couple_count {aa_k},{aa_c}: {intra_couple_count[aa_k][aa_c]}")

    print("intra_couple_count")
    print(intra_couple_count)

    cond_proba = {}
    for aa_k in list_residu:
        cond_proba[aa_k] = {}
        for aa_p in list_residu:
            cond_proba[aa_k][aa_p] = {}
            for aa_c in list_residu: 
                if intra_couple_count[aa_k][aa_c] != 0:
                    cond_proba[aa_k][aa_p][aa_c] = (triplet_count[aa_k][aa_p][aa_c]) / (intra_couple_count[aa_k][aa_c])  
                else:
                    cond_proba[aa_k][aa_p][aa_c] = 0

    return cond_proba             



def sumLine(cond_proba, list_AA, aa_k, aa_c):
    sum_line = 0
    for aa_p in list_AA:
        sum_line += cond_proba[aa_k][aa_p][aa_c]
    return sum_line


def sumPlate(dico_triple):
    for aa_1 in dico_triple:
        sum_plateau = 0
        for aa_2 in dico_triple[aa_1]:
            for aa_3 in dico_triple[aa_1][aa_2]:
                sum_plateau += dico_triple[aa_1][aa_2][aa_3]
        print("{}, {}".format(aa_1, sum_plateau))




def simpleContextualBlosum(path_file_fasta,  delay_num, kp_SeqChoice, list_residu):

    print("Context: ({},{})".format(delay_num, kp_SeqChoice))

    triplet_count = initialisation(list_residu)
    triplet_count = tripletCount(path_file_fasta, triplet_count, delay_num, kp_SeqChoice, list_residu)

    cond_proba = conditionalProba(list_residu, triplet_count)
    print("sum per level in the cube")
    sumPlate(cond_proba)
    path_proba_cond = f"verif_proba_cond_({str(delay_num)},{kp_SeqChoice})"
    np.save(path_proba_cond, cond_proba) 
    path_proba_cond = f"{path_proba_cond}.npy"

    return path_proba_cond


if __name__ == '__main__': 
    
    path_file_fasta = "/Users/pauline/Desktop/Stage_MNHN/utils_dev/Pfam_verif_cube/PF00002.27.trim"
    #delay_num, kp_SeqChoice = -1,  "k"   # <-- validated
    #delay_num, kp_SeqChoice = -2,  "k"   # <-- validated
    delay_num, kp_SeqChoice = 3,  "p"    #  <-- validated
    list_residu = ["A", "B"]
    path_dico = simpleContextualBlosum(path_file_fasta,  delay_num, kp_SeqChoice, list_residu)
    dico_loaded = loadVisu.loadAndVisualisation(path_dico)


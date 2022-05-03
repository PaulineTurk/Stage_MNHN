from pathlib import Path
import os, shutil
import pandas as pd
from math import log2
import numpy as np
import os 
from math import sqrt
import blosum as bl
from timer import Timer
from fastaReader import readFastaMul
import seaborn as sb
import matplotlib.pyplot as plt






def countCouple_countAA(num_accession, name_folder_pid, liSeqAliFiltre, pid_inf,  count_aa_global, nbre_aa_global, 
                        count_couple_aa_global, nbre_couple_aa_global, list_residu):   
    pid_couple = np.load(f"{name_folder_pid}/{num_accession}.pId.npy", allow_pickle='TRUE').item()
    nbre_seq = len(liSeqAliFiltre)
    for i in range(nbre_seq):
        name_1, seq_1 = liSeqAliFiltre[i]
        for j in range(i + 1, nbre_seq):
            name_2 ,seq_2 = liSeqAliFiltre[j] 
            if pid_couple[name_1][name_2] >= pid_inf:
                for (aa_1, aa_2) in zip(seq_1, seq_2):
                    if aa_1 in list_residu and aa_2 in list_residu:

                        # count AA
                        count_aa_global[aa_1] += 1
                        count_aa_global[aa_2] += 1
                        nbre_aa_global += 2

                        # count couple AA
                        if aa_1 == aa_2:
                            count_couple_aa_global[aa_1][aa_2] += 2
                        else:
                            count_couple_aa_global[aa_1][aa_2] += 1
                            count_couple_aa_global[aa_2][aa_1] += 1
                        nbre_couple_aa_global += 2

    return count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global 





def multiBlosum(path_folder_BlosumResX, path_folder_fasta, path_pid_folder, list_residu, pid_inf = 62, scale_factor = 2):
    t = Timer()
    t.start()


    if os.path.isdir(path_folder_BlosumResX):
        shutil.rmtree(path_folder_BlosumResX) 
    os.mkdir(path_folder_BlosumResX)


    path_folder_fasta = Path(path_folder_fasta)
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    count_aa_global = {}  

    for aa in list_residu:
        count_aa_global[aa] = 0

    nbre_aa_global = 0

    count_couple_aa_global = {}
    for aa_1 in list_residu:
        count_couple_aa_global[aa_1] = {}
        for aa_2 in list_residu:
                count_couple_aa_global[aa_1][aa_2] = 0

    nbre_couple_aa_global = 0


    for file_name_fasta in files_in_path_folder_fasta:
            accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]

            data_train = readFastaMul(file_name_fasta)
            count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global = countCouple_countAA(accession_num, 
            path_pid_folder, data_train, pid_inf, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global, list_residu)



    freq_aa_global = {}
    for aa in list_residu:
        if nbre_aa_global != 0:
            freq_aa_global[aa] = count_aa_global[aa]/nbre_aa_global
        else: 
            freq_aa_global[aa] = 0
    path_freqAA = f"{path_folder_BlosumResX}/Blosum_freq_AA"
    np.save(path_freqAA, freq_aa_global) 
    



    freq_couple_aa_global = {}
    for aa_1 in list_residu:
        freq_couple_aa_global[aa_1] = {}
        for aa_2 in list_residu:
            if nbre_couple_aa_global != 0:
                freq_couple_aa_global[aa_1][aa_2] = count_couple_aa_global[aa_1][aa_2]/nbre_couple_aa_global
            else:
                freq_couple_aa_global[aa_1][aa_2] = 0



    # compute and save the Blosum matrix 
    matrix_blosum = {}
    for aa_1 in list_residu:
        matrix_blosum[aa_1] = {}
        for aa_2 in list_residu:
            if freq_couple_aa_global[aa_1][aa_2] != 0:
                matrix_blosum[aa_1][aa_2] = round(scale_factor*log2(freq_couple_aa_global[aa_1][aa_2]/(freq_aa_global[aa_1]*freq_aa_global[aa_2])))
            else:
                matrix_blosum[aa_1][aa_2] = 0 
    df_matrix_blosum = pd.DataFrame.from_dict(matrix_blosum)   
    #print(df_matrix_blosum)
    path_matrix = f"{path_folder_BlosumResX}/Blosum_score"
    np.save(path_matrix, matrix_blosum) 


    # compute and save the matrix of conditional probabilities
    matrix_cond_proba = {}
    for aa_1 in list_residu:
        matrix_cond_proba[aa_1] = {}
        for aa_2 in list_residu:
            if freq_couple_aa_global[aa_1][aa_2] != 0:
                matrix_cond_proba[aa_1][aa_2] = freq_couple_aa_global[aa_1][aa_2]/freq_aa_global[aa_1]
            else:
                matrix_cond_proba[aa_1][aa_2] = 0

    df_matrix_cond_proba = pd.DataFrame.from_dict(matrix_cond_proba)
    df_matrix_cond_proba = np.transpose(pd.DataFrame.from_dict(matrix_cond_proba))
    #print(df_matrix_cond_proba)
    #sum_ligne = df_matrix_cond_proba.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)
    path_proba_cond = f"{path_folder_BlosumResX}/Blosum_proba_cond"
    np.save(path_proba_cond, matrix_cond_proba)


    t.stop("Compute the substitution matrix and the conditional probability matrix")
    return matrix_blosum, matrix_cond_proba, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global 






def diffBlosum(my_blosum, list_residu, blosum_ref_num = 62):
    """Quantify the distance between my_blosum and a reference blosum"""
    blosum_ref = bl.BLOSUM(blosum_ref_num) 
    matrix_diff = {}
    average_diff = 0
    count = 0
    euclidean_d = 0
    for aa1 in list_residu:
        matrix_diff[aa1] = {}
        for aa2 in list_residu:
            matrix_diff[aa1][aa2] = int(my_blosum[aa1][aa2] - blosum_ref[aa1 + aa2])
            average_diff += matrix_diff[aa1][aa2]
            euclidean_d += (matrix_diff[aa1][aa2])**2
            count += 1

    euclidean_d = round(sqrt(euclidean_d), 2)
    average_diff = round(average_diff/count, 2)
    return matrix_diff, blosum_ref_num, euclidean_d, average_diff 






# compute and save the conditional proba matrix of a BLOSUM of reference
def probaCondReference(blosum_ref_num, residu_included, path_folder_BlosumRes_ref, scale_factor = 2):
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_BlosumRes_ref):
        shutil.rmtree(path_folder_BlosumRes_ref) 
    os.mkdir(path_folder_BlosumRes_ref)

    blosum_ref = bl.BLOSUM(blosum_ref_num)
    # reverse blosum_ref
    proba_cond_blosum_ref = {}
    for aa_1 in residu_included:
        proba_cond_blosum_ref[aa_1] = {}
        for aa_2 in residu_included:
            blosum_score = blosum_ref[aa_1 + aa_2]
            if blosum_score != 0:
                proba_cond_blosum_ref[aa_1][aa_2] = np.exp(blosum_score/scale_factor)
            else:
                proba_cond_blosum_ref[aa_1][aa_2] = 0
    
    # visualisation df not normalized    (transpose to have the known aa read at each line)
    df_proba_cond_blosum_ref = np.transpose(pd.DataFrame.from_dict(proba_cond_blosum_ref))
    #print(df_proba_cond_blosum_ref)
    
    # normalisation (sum of each line = 1)
    proba_cond_blosum_ref_normalised = {}
    for aa_1 in residu_included:
        sum_line = sum(proba_cond_blosum_ref[aa_1].values())
        #print(proba_cond_blosum_ref[aa_1])
        #print(sum_line)
        proba_cond_blosum_ref_normalised[aa_1] = {k: v/sum_line for k, v in proba_cond_blosum_ref[aa_1].items()}
    t.stop("Compute the conditional proba matrix from a Blosum of reference")

    # visualisation df normalized
    df_proba_cond_blosum_ref_normalised = np.transpose(pd.DataFrame.from_dict(proba_cond_blosum_ref_normalised)) 
    #print(df_proba_cond_blosum_ref_normalised)
    sum_ligne = df_proba_cond_blosum_ref_normalised.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)
    
    # save
    path_matrix = path_folder_BlosumRes_ref + "/Blosum_proba_cond_Ref"
    np.save(path_matrix, proba_cond_blosum_ref_normalised) 


def conditionalProbaGenerator(path_data, percentage_A, path_pid, path_BlosumRes, list_residu,  name_data_train):
    path_folder_fasta_train = f"{path_data}/PfamSplit_{str(percentage_A)}/{name_data_train}"   
    print("folder_train:", path_folder_fasta_train)

    path_folder_BlosumResX = f"{path_BlosumRes}/BlosumRes_{str(percentage_A)}_{name_data_train}"
 
    matrix_blosum, matrix_cond_proba, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global  = multiBlosum(path_folder_BlosumResX, path_folder_fasta_train , path_pid, list_residu, pid_inf = 62, scale_factor = 2) 
    heatmap_cond_proba = pd.DataFrame(matrix_cond_proba).T.fillna(0)
    heatmap = sb.heatmap(heatmap_cond_proba, annot = True, annot_kws = {"size": 3}, fmt = '.1g')
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    name_folder_Res = os.path.basename(path_folder_BlosumResX)
    name_heatmap = "Heatmap Conditional proba {}".format(name_folder_Res)   
    plt.title(name_heatmap)
    plt.close()
    path_save_fig = f"{path_BlosumRes}/{name_heatmap}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)

    # difference with Blosum
    matrix_diff, blosum_ref_num, euclidean_d, average_diff = diffBlosum(matrix_blosum, list_residu, blosum_ref_num = 62)
    df_matrix_diff = pd.DataFrame(matrix_diff).T.fillna(0)
    heatmap = sb.heatmap(df_matrix_diff, annot = True, annot_kws = {"size": 5}, fmt = 'd')
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    name_folder_Res = os.path.basename(path_folder_BlosumResX)
    heatmap_title = f"Heatmap difference {name_folder_Res} - Blosum{blosum_ref_num} ref:\n euclidean distance = {euclidean_d}, average difference = {average_diff}"
    plt.title(heatmap_title)
    plt.close()
    path_save_fig = f"{path_BlosumRes}/Heatmap difference {name_folder_Res} - Blosum{blosum_ref_num} ref.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)




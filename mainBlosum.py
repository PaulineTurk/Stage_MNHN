from pathlib import Path
import os, shutil
import pandas as pd
from math import log2
import numpy as np
import os 
from math import sqrt
import blosum as bl
from FUNCTION.timer import Timer
from FUNCTION.fastaReader import readFastaMul
import FUNCTION.character as ch



def countCouple_countAA(num_accession, name_folder_pid, liSeqAliFiltre, pid_inf,  count_aa_global, nbre_aa_global, 
                        count_couple_aa_global, nbre_couple_aa_global):   
    pid_couple = np.load(name_folder_pid + "/" + num_accession + ".pId.npy", allow_pickle='TRUE').item()
    list_AA = ch.characterList()
    nbre_seq = len(liSeqAliFiltre)
    for i in range(nbre_seq):
        name_1, seq_1 = liSeqAliFiltre[i]
        for j in range(i + 1, nbre_seq):
            name_2 ,seq_2 = liSeqAliFiltre[j] 
            if pid_couple[name_1][name_2] >= pid_inf:
                for (aa_1, aa_2) in zip(seq_1, seq_2):
                    if aa_1 in list_AA and aa_2 in list_AA:

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




def oneFasta(num_accession, name_folder_pid, file_fasta, pid_inf, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global):          
    data_train = readFastaMul(file_fasta)
    count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global = countCouple_countAA(num_accession, name_folder_pid, data_train, pid_inf, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global)
    return count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global 




def multiBlosum(path_folder_BlosumResX, path_folder_fasta, path_pid_folder, pid_inf = 62, scale_factor = 2):
    t = Timer()
    t.start()


    if os.path.isdir(path_folder_BlosumResX):
        shutil.rmtree(path_folder_BlosumResX) 
    os.mkdir(path_folder_BlosumResX)


    path_folder_fasta = Path(path_folder_fasta)
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    count_aa_global = {}  
    list_AA = ch.characterList()

    for aa in list_AA:
        count_aa_global[aa] = 0

    nbre_aa_global = 0

    count_couple_aa_global = {}
    for aa_1 in list_AA:
        count_couple_aa_global[aa_1] = {}
        for aa_2 in list_AA:
                count_couple_aa_global[aa_1][aa_2] = 0

    nbre_couple_aa_global = 0


    for file_name_fasta in files_in_path_folder_fasta:
            accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
            count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global = oneFasta(accession_num, path_pid_folder, file_name_fasta, pid_inf, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global)



    freq_aa_global = {}
    for aa in list_AA:
        freq_aa_global[aa] = count_aa_global[aa]/nbre_aa_global
    path_freqAA = path_folder_BlosumResX + "/Blosum_freq_AA"
    np.save(path_freqAA, freq_aa_global) 
    



    freq_couple_aa_global = {}
    for aa_1 in list_AA:
        freq_couple_aa_global[aa_1] = {}
        for aa_2 in list_AA:
            freq_couple_aa_global[aa_1][aa_2] = count_couple_aa_global[aa_1][aa_2]/nbre_couple_aa_global


    # compute and save the Blosum matrix 
    matrix_blosum = {}
    for aa_1 in list_AA:
        matrix_blosum[aa_1] = {}
        for aa_2 in list_AA:
            if freq_couple_aa_global[aa_1][aa_2] != 0:
                matrix_blosum[aa_1][aa_2] = round(scale_factor*log2(freq_couple_aa_global[aa_1][aa_2]/(freq_aa_global[aa_1]*freq_aa_global[aa_2])))
            else:
                matrix_blosum[aa_1][aa_2] = 0 
    df_matrix_blosum = pd.DataFrame.from_dict(matrix_blosum)   
    print(df_matrix_blosum)
    path_matrix = path_folder_BlosumResX + "/Blosum_score"
    np.save(path_matrix, matrix_blosum) 


    # compute and save the matrix of conditional probabilities
    matrix_cond_proba = {}
    for aa_1 in list_AA:
        matrix_cond_proba[aa_1] = {}
        for aa_2 in list_AA:
            if freq_couple_aa_global[aa_1][aa_2] != 0:
                matrix_cond_proba[aa_1][aa_2] = freq_couple_aa_global[aa_1][aa_2]/freq_aa_global[aa_1]
            else:
                matrix_cond_proba[aa_1][aa_2] = 0

    df_matrix_cond_proba = pd.DataFrame.from_dict(matrix_cond_proba)
    df_matrix_cond_proba = np.transpose(pd.DataFrame.from_dict(matrix_cond_proba))
    print(df_matrix_cond_proba)
    #sum_ligne = df_matrix_cond_proba.sum(axis=1)
    #print("Somme des lignes:\n", sum_ligne)
    path_proba_cond = path_folder_BlosumResX + "/Blosum_proba_cond"
    np.save(path_proba_cond, matrix_cond_proba)


    t.stop("Compute the substitution matrix and the conditional probability matrix")
    return matrix_blosum, matrix_cond_proba, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global 






def diffBlosum(my_blosum, blosum_ref_num = 62):
    """Quantify the distance between my_blosum and a reference blosum"""
    blosum_ref = bl.BLOSUM(blosum_ref_num) 
    euclidiean_d = 0
    matrix_diff = {}
    list_AA = ch.characterList()
    for aa1 in list_AA:
        matrix_diff[aa1] = {}
        for aa2 in list_AA:
            euclidiean_d += (my_blosum[aa1][aa2] - blosum_ref[aa1 + aa2])**2
            matrix_diff[aa1][aa2] = my_blosum[aa1][aa2] - blosum_ref[aa1 + aa2]
    print("Euclidean distance: {}".format(round(sqrt(euclidiean_d), 2)))
    print("My Blosum - The reference Blosum{} :".format(blosum_ref_num))
    df_matrix_diff = pd.DataFrame.from_dict(matrix_diff)
    print(df_matrix_diff)

















if __name__ == '__main__': 

    #name_matrix = "matrice_nogaps_1_seed"        # 0.19161   s    # Distance euclidienne: 48.32
    #name_matrix = "matrice_nogaps_10_seed"       # 13.0144   s    # Distance euclidienne: 34.44
    #name_matrix = "matrice_nogaps_100_seed"      # 27.53306  s    # Distance euclidienne: 26.12
    #name_matrix = "matrice_nogaps_1000_seed"     # 342.20251 s    # Distance euclidienne: 27.46
    #name_matrix = "matrice_nogaps_10000_seed"    # 3279.8419 s    # Distance euclidienne: 28.11
    #folder_fasta = "Pfam_fasta_trimAl_nogaps_header_corrected_non_redundant"



    path_pid_folder = "/Users/pauline/Desktop/data/PID_couple"
    # /Users/pauline/Desktop/data/BlosumRes    à créer en amont



    # 0.05 % Pfam
    path_folder_BlosumResX = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_0.05"  
    path_folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_0.05/PfamTrain"   # 9 seeds,  0.80627 s, euclidien distance: 77.16

    # 0.5 % Pfam
    #path_folder_BlosumResX = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_0.5"  
    #path_folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_0.5/PfamTrain"     # 91 seeds,  5.86692 s, euclidien distance: 37.34

    # 5.0 % Pfam
    #path_folder_BlosumResX = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_5.0"  
    #path_folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_5.0/PfamTrain"     # 910 seeds, 142.45058 s, euclidien distance: 36.04
    
    # 50.0 % Pfam
    #path_folder_BlosumResX = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_50.0"  
    #path_folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_50.0/PfamTrain"    # 9101 seeds, 1248.79718 s, euclidien distance: 32.28



    # manual check   #pid_inf = 62 pour le calcul de proba_cond et freq
    #path_folder_BlosumResX = "blosum_res_test"  
    #path_folder_fasta = "Pfam_test_trimmed_manuel"   

    matrix_blosum, matrix_cond_proba, count_aa_global, nbre_aa_global, count_couple_aa_global, nbre_couple_aa_global  = multiBlosum(path_folder_BlosumResX, path_folder_fasta , path_pid_folder, pid_inf = 62, scale_factor = 2) 
  
    diffBlosum(matrix_blosum, blosum_ref_num = 62)


from pathlib import Path
import os 
import numpy as np
import pandas as pd
from FUNCTION.timer import Timer
import FUNCTION.brierPredictorNeighbor as brn
from FUNCTION.fastaReader import readFastaMul


def multiBrier01(folder_fasta, dir_pid_name, pid_inf = 62):   
    t = Timer()
    t.start()
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()
    Brier_count_global = 0
    count_global = 0

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        data_test = readFastaMul(file_name_fasta)
        Brier_count_global, count_global  = brn.predictor01(data_test, accession_num, dir_pid_name, 
                                                                  Brier_count_global, count_global, pid_inf)
    Brier_Score_global = Brier_count_global/count_global
    t.stop("Brier Score with predictor 0/1 (worst case scenario)")
    return Brier_Score_global


def multiBrierPerfect(folder_fasta, dir_pid_name, pid_inf = 62):   
    t = Timer()
    t.start()
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()
    Brier_count_global = 0
    count_global = 0

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        data_test = readFastaMul(file_name_fasta)
        Brier_count_global, count_global  = brn.predictorPerfect(data_test, accession_num, dir_pid_name, 
                                                                  Brier_count_global, count_global, pid_inf)
    Brier_Score_global = Brier_count_global/count_global
    t.stop("Brier Score with perfect predictor")
    return Brier_Score_global


def multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier, delay_num, kp_SeqChoice, pid_inf = 62):  
    """
    predicteur_name is taken from the list: 
    ["Blosum Predictor", "Equiprobable Predictor",  "Stationary Predictor", "Identity Predictor"]
    """
    t = Timer()
    t.start()
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()
    
    Brier_count_global = 0
    count_global = 0

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        data_test = readFastaMul(file_name_fasta)
        Brier_count_global, count_global = brn.brierMatrixNeighbor(predictor_name, unit_Brier, data_test, accession_num, dir_pid_name, 
                                                                 Brier_count_global, count_global, delay_num, kp_SeqChoice, pid_inf)

    Brier_Score_global = Brier_count_global/count_global
    t.stop("Brier Score with {}".format(predictor_name))
    return Brier_Score_global




if __name__ == '__main__': 
    folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_0.5/PfamTrain"
    dir_pid_name = "/Users/pauline/Desktop/data/PID_couple"

    ########################## Extreme, non-matrix predictors
    # 01 predictor
    #worst_brier_score = multiBrier01(folder_fasta, dir_pid_name, pid_inf = 62)   
    #print(worst_brier_score)
    ##### nbre_block = 10    # Score: 2.0    # 88.26048 s

    # Perfect predictor
    #best_brier_score = multiBrierPerfect(folder_fasta, dir_pid_name, pid_inf = 62)  
    #print(best_brier_score)
    ##### nbre_block = 10    # Score: 0.0    # 41.22925 s 



    ########################## Matrix predictors
    # Equiprobable Predictor
    #predictor_name, cond_proba_equiproba, unit_Brier_equiproba = br.predicteurEquiprobable()
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_equiproba, pid_inf = 62) 
    #print(Brier_Score_global)
    ##### nbre_block = 10      # Score: 0.9499999995757045    # 7.97453 s


    # Identity Predictor
    #predictor_name, cond_proba_id, unit_Brier_id = br.predictorIdentity()
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_id, pid_inf = 62) 
    #print(Brier_Score_global)
    ##### nbre_block = 10      # Score: 1.2142567222866654    # 8.62109 s
    ##### nbre_block = 100     # Score: 1.1395143108363937    # 17.51399 s
    ##### nbre_block = 1000    # Score: 1.0800426063136361    # 222.40938 s


    # Stationary Predictor
    #freq_aa_global_path = "matrix/matrice_nogaps_10_seed_freq_AA.npy"

    #freq_aa_global = np.load(freq_aa_global_path, allow_pickle='TRUE').item()
    #predictor_name, cond_proba_stationary, unit_Brier_stationary = br.predictorStationary(freq_aa_global)
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_stationary, block_nbre, pid_inf = 62) 
    #print(Brier_Score_global)
    ##### nbre_block = 10      # Score: 0.9345482772799836    # 7.72772 s
    ##### nbre_block = 100     # Score: 0.935110285053936     # 17.34066 s
    ##### nbre_block = 1000    # Score: 0.9393303237972899    # 203.5844 s



    # Blosum Predictor  
    path_matrix = "/Users/pauline/Desktop/data/NeighborRes/NeighborRes_0.5/proba_cond_(1 , p).npy"
    cond_proba = np.load(path_matrix, allow_pickle='TRUE').item()
    df_cond_proba = np.transpose(pd.DataFrame.from_dict(cond_proba))  # only for visualisation
    print(df_cond_proba)                                              # only for visualisation
    delay_num = 1 
    #delay_num = 1 
    #kp_SeqChoice = "k"
    kp_SeqChoice = "p"
    predictor_name, cond_proba_blosum, unit_Brier_Blosum = brn.predictorBlosumNeighbor(path_matrix)
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_Blosum,
                                         delay_num, kp_SeqChoice) 
    print(Brier_Score_global)







    ##### nbre_block = 10      # Score: 0.8011259348825893    # 8.95182 s
    ##### nbre_block = 100     # Score: 0.7699080119866685    # 17.15007 s
    ##### nbre_block = 1000    # Score: 0.7445928173406464    # 213.03552 s

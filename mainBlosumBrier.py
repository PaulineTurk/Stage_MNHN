from FUNCTION.fastaReader import readFastaMul
from pathlib import Path
from FUNCTION.timer import Timer
import os 
import numpy as np
import FUNCTION.brierPredictor as br





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
        Brier_count_global, count_global  = br.predictor01(data_test, accession_num, dir_pid_name, 
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
        Brier_count_global, count_global  = br.predictorPerfect(data_test, accession_num, dir_pid_name, 
                                                                  Brier_count_global, count_global, pid_inf)
    Brier_Score_global = Brier_count_global/count_global
    t.stop("Brier Score with perfect predictor")
    return Brier_Score_global


def multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier, pid_inf = 62):  
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
        #print(accession_num)
        data_test = readFastaMul(file_name_fasta)
        Brier_count_global, count_global = br.brierMatrix(predictor_name, unit_Brier, data_test, accession_num, dir_pid_name, 
                                                                 Brier_count_global, count_global, pid_inf)
    Brier_Score_global = Brier_count_global/count_global
    t.stop("Brier Score with {}".format(predictor_name))
    return Brier_Score_global




if __name__ == '__main__': 

    folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_0.05/PfamTrain" 
    #folder_fasta = "Pfam_test_trimmed_manuel" 
    print("### Data TEST:", folder_fasta)
    print("")
    dir_pid_name = "/Users/pauline/Desktop/data/PID_couple"

    ########################## Extreme, non-matrix predictors
    # 01 predictor
    worst_brier_score = multiBrier01(folder_fasta, dir_pid_name, pid_inf = 62)   
    #worst_brier_score = multiBrier01(folder_fasta, dir_pid_name, pid_inf = 0)   
    print("worst_brier_score:", worst_brier_score)
    print("")

    # Perfect predictor
    best_brier_score = multiBrierPerfect(folder_fasta, dir_pid_name, pid_inf = 62)  
    #best_brier_score = multiBrierPerfect(folder_fasta, dir_pid_name, pid_inf = 0)  
    print("best_brier_score:", best_brier_score)
    print("")



    ########################## Matrix predictors
    # Equiprobable Predictor
    predictor_name, cond_proba_equiproba, unit_Brier_equiproba = br.predicteurEquiprobable()
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_equiproba, pid_inf = 62) 
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_equiproba, pid_inf = 0) 
    print("Equiprobable Predictor Brier Score:", Brier_Score_global)
    print("")



    # Identity Predictor
    predictor_name, cond_proba_id, unit_Brier_id = br.predictorIdentity()
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_id, pid_inf = 62) 
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_id, pid_inf = 0) 
    print("Identity Predictor Brier Score:", Brier_Score_global)
    print("")


    # Stationary Predictor
    path_freq_aa_global = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_0.05/Blosum_freq_AA.npy"
    #path_freq_aa_global = "blosum_res_test/Blosum_freq_AA.npy"
    print("### Data TRAIN:", path_freq_aa_global)
    freq_aa_global = np.load(path_freq_aa_global, allow_pickle='TRUE').item()
    predictor_name, cond_proba_stationary, unit_Brier_stationary = br.predictorStationary(freq_aa_global)
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_stationary, pid_inf = 62) 
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_stationary, pid_inf = 0) 
    print("Stationary Predictor Brier Score:", Brier_Score_global)
    print("")




    # Blosum Predictor  
    path_matrix_cond_proba = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_0.05/Blosum_proba_cond.npy"
    #path_matrix_cond_proba = "blosum_res_test/Blosum_proba_cond.npy"
    print("### Data TRAIN:", path_matrix_cond_proba)
    predictor_name, cond_proba_blosum, unit_Brier_Blosum = br.predictorBlosum(path_matrix_cond_proba)
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_Blosum, pid_inf = 62) 
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_Blosum, pid_inf = 0) 
    print("Blosum Predictor Brier Score:", Brier_Score_global)
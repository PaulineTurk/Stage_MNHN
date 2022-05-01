from fastaReader import readFastaMul
import brierPredictor as br
from timer import Timer
from pathlib import Path
import os 
import numpy as np






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
    if count_global != 0:
        Brier_Score_global = Brier_count_global/count_global
    else:
        Brier_Score_global = None
        print("Brier Score not computable because count_gloable = 0")
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


def multiBrierMatrix(predictor_name, path_folder_fasta, dir_pid_name, unit_Brier, list_residu, pid_inf = 62):  
    """
    predicteur_name is taken from the list: 
    ["Blosum Predictor", "Equiprobable Predictor",  "Stationary Predictor", "Identity Predictor"]
    """
    t = Timer()
    t.start()

    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()
    
    Brier_count_global = 0
    count_global = 0

    for file_name in files_in_path_folder_fasta:
        accession_num_part_1 = os.path.basename(file_name).split(".")[0]
        accession_num_part_2 = os.path.basename(file_name).split(".")[1]
        accession_num = f"{accession_num_part_1}.{accession_num_part_2}"
        data_test = readFastaMul(file_name)
        Brier_count_global, count_global = br.brierMatrix(predictor_name, unit_Brier, data_test, accession_num, dir_pid_name, 
                                                                 Brier_count_global, count_global, list_residu, pid_inf)
    if count_global != 0:
        Brier_Score_global = Brier_count_global/count_global
    else:
        Brier_Score_global = None
    t.stop("Brier Score with {}".format(predictor_name))
    return Brier_Score_global




def overfittingTest(path_data, percentage_A, path_pid, path_BlosumRes, name_data_train, name_data_test, list_residu):
    folder_fasta_train = f"{path_data}/PfamSplit_{str(percentage_A)}/{name_data_train}" 
    folder_fasta_test = f"{path_data}/PfamSplit_{str(percentage_A)}/{name_data_test}" 
    print("folder_fasta_train:", folder_fasta_train)
    print("folder_fasta_test:", folder_fasta_test)

    path_matrix_cond_proba = f"{path_BlosumRes}/BlosumRes_{str(percentage_A)}_{name_data_train}/Blosum_proba_cond.npy" 


    predictor_name, cond_proba_blosum, unit_Brier_Blosum = br.predictorBlosum(path_matrix_cond_proba, list_residu)
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta_test, path_pid, unit_Brier_Blosum, list_residu, pid_inf = 62) 
    print("Blosum Predictor Brier Score:", Brier_Score_global)
    print("")

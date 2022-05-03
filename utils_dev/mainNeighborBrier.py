from pathlib import Path
import os 
import numpy as np
import pandas as pd
from timer import Timer
import brierPredictorNeighbor as brn
from fastaReader import readFastaMul


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



def bierContextBayes(folder_fasta, dir_pid_name, pid_inf = 62):   
    t = Timer()
    t.start()
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()
    Brier_count_global = 0
    count_global = 0

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        data_test = readFastaMul(file_name_fasta)
        Brier_count_global, count_global  = brn.predictorContextBayes(data_test, accession_num, dir_pid_name, 
                                                                  Brier_count_global, count_global, pid_inf)
    Brier_Score_global = Brier_count_global/count_global
    t.stop("Brier Score with bayes context")
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


from pickle import TRUE
import dataSplit
import os, shutil
import redundancy
import Stockholm_Fasta
import PID_save as pid
import redundancy
import numpy as np
import pandas as pd
import mainBlosum
import mainBlosumBrier
import lowerToUpper
import dataCountDescription
import mainNeighbor
import ContextuelBrier


# true/false à choisir: data_treatment, new_folder (each new datasplit), data_split



#Initial configuration
# manual initialisation of a folder containing Pfam-A.seed downloaded from Pfam
path_main_folder =  "/Users/pauline/Desktop/Overfitting_test" 
#path_main_folder = "/Users/pauline/Desktop/data_Test_3"
name_file_from_Pfam = "Pfam-A.seed"

path_file_from_Pfam = path_main_folder + "/" + name_file_from_Pfam


#####################  User choices
# data treatment  (do it once !)
data_treatment = False
name_folder_stockholm = "Pfam_Stockholm"
name_folder_fasta = "Pfam_fasta"
name_folder_upper = "Pfam_upper"
name_folder_pid = "PID_couple"
name_fasta_non_redondant = "Pfam_fasta_99"
pid_sup = 99
list_AA = ch.characterList()     # defined for all the reste of the study


# description Pfam after data treatment (éventuellement le voir à chaque étape)
#descriptionPfam = True

# data_split versionning
new_folder = False
name_new_folder =  "test_1" # check that the name is not already taken
list_percentage = [50]  


# blosum generator on data_train and test (i.e split A and B)
blosumGenerator = False
name_BlosumRes = "BlosumRes"

# simple context blosum generator (i.e cubes of conditional probabilities)
cube_generator = False


##################### Brier Score (over-fitting part)
blosum_overfitting_test = False




##################### Brier Score with simple context (Naive Bayes)
brier_Score_naive_Bayes = True
path_folder_fasta_test = "/Users/pauline/Desktop/Overfitting_test/test_1/PfamSplit_50/Pfam_mini_test_context"
percentage_train = 50
path_pid_folder = "/Users/pauline/Desktop/Overfitting_test/PID_couple"
path_NeighborRes = "/Users/pauline/Desktop/Overfitting_test/test_1/NeighborRes"
list_list_len_window = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]    # [len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p]






##################### data treatment
path_folder_stockholm = path_main_folder + "/" + name_folder_stockholm
path_folder_fasta = path_main_folder + "/" + name_folder_fasta
path_folder_upper = path_main_folder + "/" + name_folder_upper
path_folder_pId = path_main_folder + "/" + name_folder_pid
path_folder_fasta_non_redondant =  path_main_folder + "/" + name_fasta_non_redondant 

if data_treatment == True:
    # separation of the Stockholm file into Stockholm files
    Stockholm_Fasta.separationStockholm(path_file_from_Pfam, path_folder_stockholm) 
    fileNumber.countFile(path_folder_stockholm)  

    # conversion from stockholm into fasta files
    Stockholm_Fasta.multiStockholmToFasta(path_folder_fasta, path_folder_stockholm)  
    dataCountDescription.dataCountDescription(path_folder_fasta)

    # conversion of all the residus lower case into upper case
    lowerToUpper.multiLowerToUpper(path_folder_fasta, path_folder_upper)
    dataCountDescription.dataCountDescription(path_folder_upper)

    # pid couple calculation
    pid.savePId(path_folder_upper, path_folder_pId)   

    # visual check
    #pid_check = np.load("/Users/pauline/Desktop/MINI_TEST/pid_couple/PF00244.23.pId.npy" ,allow_pickle='TRUE').item()  
    #df_pid_check = np.transpose(pd.DataFrame.from_dict(pid_check))
    #print(df_pid_check)  

    # redundancy issue  
    redundancy.savePIdNonRedondant(path_folder_upper, path_folder_fasta_non_redondant, pid_sup, list_AA)  # 7272.81967 s
    dataCountDescription.dataCountDescription(path_folder_fasta_non_redondant)


# description Pfam after data treatment
#if descriptionPfam == True:
#    dataCountDescription.dataCountDescription(path_folder_fasta_non_redondant)





##################### version part
path_new_folder = path_main_folder + "/" + name_new_folder
if new_folder == True:
    if os.path.isdir(path_new_folder ):
        shutil.rmtree(path_new_folder ) 
    os.mkdir(path_new_folder)
   
    name_folder_data_train = "PfamTrain"
    name_folder_data_test = "PfamTest"
    path_folder_fasta_non_redondant =  path_main_folder + "/" + name_fasta_non_redondant 
 
    for percentage_train in list_percentage:
        name_folder_data_train_test = "PfamSplit" + "_" + str(percentage_train)
        path_folder_data_train_test = path_new_folder + "/" + name_folder_data_train_test 
        dataSplit.trainTestSplit(path_folder_fasta_non_redondant, path_folder_data_train_test, name_folder_data_train, name_folder_data_test, percentage_train) 









##################### blosum generator
path_BlosumRes = path_new_folder + "/" + name_BlosumRes 

if blosumGenerator == True:
    if os.path.isdir(path_BlosumRes):
        shutil.rmtree(path_BlosumRes) 
    os.mkdir(path_BlosumRes)
    for percentage_train in list_percentage:
        mainBlosum.conditionalProbaGenerator(path_new_folder, percentage_train, path_folder_pId, path_BlosumRes, train_test_reverse = False)    
        mainBlosum.conditionalProbaGenerator(path_new_folder, percentage_train, path_folder_pId, path_BlosumRes, train_test_reverse = True)




##################### cube generator (i.e blosum cond proba with one contextual aa +/- near)
if cube_generator == True:
    list_delay_number = [-1, 1]
    name_NeighborRes = "NeighborRes"
    path_folder_pid = "/Users/pauline/Desktop/data/PID_couple"
    percentage_train = 50
    path_folder_fasta_train = path_new_folder + "/PfamSplit_" + str(percentage_train) +  "/PfamTrain"

    path_NeighborRes = path_new_folder + "/" + name_NeighborRes


    for delay_num in list_delay_number:
        for kp_SeqChoice in ["k", "p"]:
            print("{}, {}".format(delay_num, kp_SeqChoice))
            path_proba_cond =  mainNeighbor.simpleContextualBlosum(path_folder_fasta_train, percentage_train, path_folder_pid, path_NeighborRes, delay_num, kp_SeqChoice, pid_inf = 62, scale_factor = 2)
  

##################### Brier Score (over-fitting part)
if blosum_overfitting_test == True:
    for percentage_train in list_percentage:
        for test_is_train in [True, False]:
            for train_test_reverse in [True, False]:
                print("test_is_train:", test_is_train)
                print("train_test_reverse:", train_test_reverse)
                mainBlosumBrier.overfittingTest(path_new_folder, percentage_train, test_is_train, train_test_reverse, path_folder_pId, path_BlosumRes)





##################### Brier Score with simple context (Naive Bayes)
if brier_Score_naive_Bayes == True:
    for list_len_window in list_list_len_window:
        list_bloc = []
        for position, len_window in enumerate(list_len_window):
            if len_window != 0:
                list_neighborResSelection_name, list_neighborResSelection = ContextuelBrier.neighborResSelection(position, len_window, percentage_train, path_NeighborRes)         
                list_bloc.append(list_neighborResSelection)
            else:
                list_bloc.append([])
        ContextuelBrier.multriContextBayes(path_folder_fasta_test, path_pid_folder, path_NeighborRes, list_len_window, list_bloc, pid_inf = 62)
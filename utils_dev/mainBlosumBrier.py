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
    if count_global != 0:
        Brier_Score_global = Brier_count_global/count_global
    else:
        Brier_Score_global = None
    t.stop("Brier Score with {}".format(predictor_name))
    return Brier_Score_global



# overfitting ?

def overfittingTest(path_data, percentage_train, test_is_train, train_test_reverse, path_pid, path_BlosumRes):

    #print("percentage_train:", percentage_train)
    # intial data_train/test
    folder_fasta_train = path_data + "/PfamSplit_" + str(percentage_train) + "/PfamTrain"   # ATTENTION changer l'appel, le généraliser plus
    folder_fasta_test = path_data + "/PfamSplit_" +  str(percentage_train) + "/PfamTest" 
    #print("folder_fasta_train:", folder_fasta_train)
    #print("folder_fasta_test:", folder_fasta_test)

    # possible combination of data_train/test
    if train_test_reverse == False:
        folder_train = folder_fasta_train
        path_matrix_cond_proba = path_BlosumRes + "/BlosumRes_" + str(percentage_train) + "_A" +  "/Blosum_proba_cond.npy"   # à calculer les matrices nécessaires +  adapter les paths
        if test_is_train == True:
            folder_test = folder_train        
        else:
            folder_test = folder_fasta_test
    else:
        folder_train = folder_fasta_test
        path_matrix_cond_proba = path_BlosumRes + "/BlosumRes_" + str(percentage_train) + "_B" + "/Blosum_proba_cond.npy"   # à calculer les matrices nécessaires +  adapter les paths

        if test_is_train == True:
            folder_test = folder_train
        else:
            folder_test = folder_fasta_train

    print("folder_train:", folder_train)
    print("folder_test:", folder_test)

    predictor_name, cond_proba_blosum, unit_Brier_Blosum = br.predictorBlosum(path_matrix_cond_proba)
    Brier_Score_global = multiBrierMatrix(predictor_name, folder_test, path_pid, unit_Brier_Blosum, pid_inf = 62) 
    print("Blosum Predictor Brier Score:", Brier_Score_global)
    print("")



if __name__ == '__main__': 

    #folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_0.05/PfamTrain" 
    #folder_fasta = "Pfam_test_trimmed_manuel" 
    #print("### Data TEST:", folder_fasta)
    #print("")
    #dir_pid_name = "/Users/pauline/Desktop/data/PID_couple"

    ########################## Extreme, non-matrix predictors
    # 01 predictor
    #worst_brier_score = multiBrier01(folder_fasta, dir_pid_name, pid_inf = 62)     
    #print("worst_brier_score:", worst_brier_score)
    #print("")

    # Perfect predictor
    #best_brier_score = multiBrierPerfect(folder_fasta, dir_pid_name, pid_inf = 62)   
    #print("best_brier_score:", best_brier_score)
    #print("")



    ########################## Matrix predictors
    # Equiprobable Predictor
    #predictor_name, cond_proba_equiproba, unit_Brier_equiproba = br.predicteurEquiprobable()
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_equiproba, pid_inf = 62) 
    #print("Equiprobable Predictor Brier Score:", Brier_Score_global)
    #print("")



    # Identity Predictor
    #predictor_name, cond_proba_id, unit_Brier_id = br.predictorIdentity()
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_id, pid_inf = 62) 
    #print("Identity Predictor Brier Score:", Brier_Score_global)
    #print("")


    # Stationary Predictor
    #path_freq_aa_global = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_0.05/Blosum_freq_AA.npy"
    #path_freq_aa_global = "blosum_res_test/Blosum_freq_AA.npy"
    #print("### Data TRAIN:", path_freq_aa_global)
    #freq_aa_global = np.load(path_freq_aa_global, allow_pickle='TRUE').item()
    #predictor_name, cond_proba_stationary, unit_Brier_stationary = br.predictorStationary(freq_aa_global)
    #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_stationary, pid_inf = 62) 
    #print("Stationary Predictor Brier Score:", Brier_Score_global)
    #print("")





















    # Blosum Predictor 
    #dir_pid_name = "/Users/pauline/Desktop/data/PID_couple"
    #list_percentage = [0.05, 0.5, 5]   # 50 not already done
    #for percentage in list_percentage:
        #print("")
        #print(percentage)
        #folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_" + str(percentage) + "/PfamTrain" 
        #path_matrix_cond_proba = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_" + str(percentage) + "/Blosum_proba_cond.npy"
        #predictor_name, cond_proba_blosum, unit_Brier_Blosum = br.predictorBlosum(path_matrix_cond_proba)
        #Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_Blosum, pid_inf = 62) 
        #print("Blosum Predictor Brier Score:", Brier_Score_global)

    # résultat sur données non trimmée et avec data_test = data_train:
    # 0.05: 1.31896 s, 0.7911168602683646
    # 0.5: 100.05061 s, 0.8696871995814536
    # 5: 810.4333 s, 0.8059479425896071






    path_data = "/Users/pauline/Desktop/data"
    percentage_train = 0.05
    path_pid = "/Users/pauline/Desktop/data/PID_couple"
    path_BlosumRes = "/Users/pauline/Desktop/data/BlosumResTEST"
    version = 0
    for test_is_train in [True, False]:
        for train_test_reverse in [True, False]:
            print("test_is_train:", test_is_train)
            print("train_test_reverse:", train_test_reverse)
            overfittingTest(path_data, percentage_train, test_is_train, train_test_reverse, path_pid, path_BlosumRes)



    # Blosum Predictor avec BLOSUM62 de référence
    #dir_pid_name = "/Users/pauline/Desktop/data/PID_couple"
    #list_percentage = [0.05, 0.5, 5]   # 50 not already done
    #for percentage in list_percentage:
    #    print("")
    #    print(percentage)
    #    folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_" + str(percentage) + "/PfamTrain" 
    #    path_matrix_cond_proba = "/Users/pauline/Desktop/data/BlosumRes/BlosumRes_ref62/Blosum_proba_cond_Ref.npy"
    #    predictor_name, cond_proba_blosum, unit_Brier_Blosum = br.predictorBlosum(path_matrix_cond_proba)
    #    Brier_Score_global = multiBrierMatrix(predictor_name, folder_fasta, dir_pid_name, unit_Brier_Blosum, pid_inf = 62) 
    #    print("Blosum Predictor Brier Score:", Brier_Score_global)

    # résultat sur données non trimmée et avec data_test = data_train utilisé pour calculer Blosum meme si 
    # on utilise cette fois-ci (BLOSUM62) Rq. ce n'étais pas la peine de recalculer les unit de Brier pour
    # chaque data_Test comme on part tjs de BOSUM62 quel que soit le data_test testé:
    # 0.05: 1.67675 s, 0.8459335064120398
    # 0.5: 116.56326 s, 0.9929872469534892
    # 5: 949.35 s, 0.8666057902074663’
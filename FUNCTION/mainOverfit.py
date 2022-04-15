import mainBlosumBrier
import brierPredictor




def overfittingTest(path_data, percentage_train, test_is_train, train_test_reverse, path_pid, path_matrix_cond_proba):

    print(percentage_train)
    # intial data_train/test
    folder_fasta_train = path_data + "/PfamSplit_" + str(percentage_train) + "/PfamTrain" 
    folder_fasta_test = path_data + "/PfamSplit_" +  str(percentage_train) + "/PfamTest" 
    print("folder_fasta_train:", folder_fasta_train)
    print("folder_fasta_test:", folder_fasta_test)

    # possible combination of data_train/test
    if train_test_reverse == False:
        folder_train = folder_fasta_train
        if test_is_train == True:
            folder_test = folder_train
        else:
            folder_test = folder_fasta_test
    else:
        folder_train = folder_fasta_test
        if test_is_train == True:
            folder_test = folder_test
        else:
            folder_test = folder_fasta_train

    print("folder_train:", folder_train)
    print("folder_test:", folder_test)

    path_matrix_cond_proba = path_BlosumRes + "/BlosumRes_" + str(percentage_train) + "/Blosum_proba_cond.npy"   # à calculer les matrices nécessaires +  adapter les paths
    predictor_name, cond_proba_blosum, unit_Brier_Blosum = brierPredictor.predictorBlosum(path_matrix_cond_proba)
    Brier_Score_global = mainBlosumBrier.multiBrierMatrix(predictor_name, folder_test, path_pid, unit_Brier_Blosum, pid_inf = 62) 
    print("Blosum Predictor Brier Score:", Brier_Score_global)
    print("")


if __name__ == '__main__': 
    path_data = "/Users/pauline/Desktop/data"
    percentage_train = 0.05
    path_pid = "/Users/pauline/Desktop/data/PID_couple"
    path_BlosumRes = "/Users/pauline/Desktop/data/BlosumRes"
    for test_is_train in [True, False]:
        for train_test_reverse in [True, False]:
            print("test_is_train:", test_is_train)
            print("train_test_reverse:", train_test_reverse)
            overfittingTest(path_data, percentage_train, test_is_train, train_test_reverse, path_pid, path_BlosumRes)

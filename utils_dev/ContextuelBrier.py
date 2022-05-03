import fastaReader
import numpy as np
import timer
import os
from path import Path
from pathlib import Path
import matplotlib.pyplot as plt



def contextCatcher(seq_k, seq_p, len_seq,  position_k, len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p):
    # seq_k et seq_p après selection selon PID et diff name
    # len_seq = len(seq_k) = len(seq_p)
    # position_k: position de aa_k dans seq_k

    if (position_k - len_window_left_k) >= 0 and (position_k + len_window_right_k +1 <= len_seq) and (position_k - len_window_left_p >= 0) and (position_k + len_window_right_p +1 <= len_seq):   
        if len_window_left_k > 0:
            window_left_k = seq_k[position_k - len_window_left_k: position_k][::-1] 
            #[::-1]    # reverse the order so the residus are ordered by increasing distance from the couple of interest
        else:
            window_left_k = []

        if len_window_right_k > 0:
            window_right_k = seq_k[position_k +1: position_k + len_window_right_k +1]
        else:
            window_right_k = []

        if len_window_left_p > 0:
            window_left_p = seq_p[position_k - len_window_left_p: position_k][::-1] 
        else:
            window_left_p = []

        if len_window_right_p > 0:
            window_right_p = seq_p[position_k +1: position_k + len_window_right_p +1]
        else:
            window_right_p = []

    else:
        window_left_k = []
        window_right_k = [] 
        window_left_p = [] 
        window_right_p = []

    return window_left_k, window_right_k, window_left_p, window_right_p


# probability vector computed for each aligned couple of valid aa
def neighborResSelection(position, len_window, percentage_train, path_NeighborRes):
    list_neighborResSelection = []
    list_neighborResSelection_name = []

    list_delay_num = [k for k in range(1, len_window + 1)]
    kp_SeqChoice = "p"
    if position in [0, 2]:
        list_delay_num = [-k for k in list_delay_num]
    if position in [0, 1]:
        kp_SeqChoice = "k"
    for delay in list_delay_num:
        path_proba_cond = f"{path_NeighborRes}/proba_cond_({str(delay)},{kp_SeqChoice})_percentage_train_{str(percentage_train)}.npy"
        if os.path.isfile(path_proba_cond):
            list_neighborResSelection_name.append(path_proba_cond)
            cond_proba = np.load(path_proba_cond , allow_pickle='TRUE').item()
            list_neighborResSelection.append(cond_proba)


    # TODO: stoper le score de brier sinon pour ce contexte !!!!!!!!!!!!!!!!!

    return list_neighborResSelection_name, list_neighborResSelection




def probaVector(seq_k, seq_p, index, list_bloc, list_residu,  list_len_window):
    len_window_left_k = list_len_window[0] 
    len_window_right_k = list_len_window[1]  
    len_window_left_p = list_len_window[2]   
    len_window_right_p = list_len_window[3]

    vector_proba = {}

    len_seq = len(seq_k)
    if seq_k[index] in list_residu and seq_p[index] in list_residu:
        window_left_k, window_right_k, window_left_p, window_right_p = contextCatcher(seq_k, seq_p, len_seq, index, len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p)
        
        list_window = [window_left_k, window_right_k, window_left_p, window_right_p]

        for position, window in enumerate(list_window):
            for index_window, contextual_residu in enumerate(window):
                if contextual_residu in list_residu:
                    for aa_x in list_residu:
                        vector_proba[aa_x] = 1
            
                    for aa_x in list_residu:   # to compute the vector of proba for the prediction                    
                        bloc = list_bloc[position][index_window] 
                        vector_proba[aa_x] *= bloc[seq_k[index]][aa_x][contextual_residu]
        
        normalisation_term = sum(vector_proba.values())
        #print("before normalisation check sum = ? : ", sum(vector_proba.values()))
        #print("before normalisation:", vector_proba)
        vector_proba = {k: v / normalisation_term for k, v in vector_proba.items()}
        #print("after normalisation check sum = ? : ", sum(vector_proba.values()))

    #print("after normalisation:", vector_proba)
    return vector_proba



def predictorContextBayes(data_test, accession_num, path_pid_folder, Brier_count_global, count_global, count_seed, 
                          list_bloc, list_len_window, pid_inf, list_residu): 
    pid_couple = np.load(f"{path_pid_folder}/{accession_num}.pId.npy", allow_pickle='TRUE').item()  

    count_seed += 1
    if data_test:
        for name_k, seq_k in data_test:
            len_seq = len(seq_k)
            for name_p, seq_p in data_test:
                if name_k != name_p:
                    if pid_couple[name_k][name_p] >= pid_inf:
                        for index in range(len_seq):
                            vector_proba = probaVector(seq_k, seq_p, index, list_bloc, list_residu,  list_len_window)
                            if vector_proba != {}:
                                count_global += 1
                                for aa_x in list_residu:
                                    Brier_count_global += (vector_proba[aa_x] - int(seq_p[index] == aa_x))**2
    return Brier_count_global, count_global, count_seed


                        


    

def multriContextBayes(path_folder_fasta, path_pid_folder, path_NeighborRes, list_len_window, list_bloc, list_residu, list_list_len_window, pid_inf = 62):   
    t = timer.Timer()
    t.start()
 
    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()
    Brier_count_global = 0
    count_global = 0
    count_seed = 0

    list_count_seed = []
    list_Brier_score = []

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        data_test = fastaReader.readFastaMul(file_name_fasta)
        Brier_count_global, count_global, count_seed  = predictorContextBayes(data_test, accession_num, path_pid_folder, 
                                                                  Brier_count_global, count_global, count_seed, list_bloc, list_len_window, pid_inf, list_residu)
           
        if count_global != 0:      
            list_Brier_score.append(Brier_count_global/count_global)   
            list_count_seed.append(count_seed)              

    Brier_Score_global = Brier_count_global/count_global

    plt.scatter(list_count_seed, list_Brier_score, label = list_len_window, alpha = 0.75, s = 5)
    plt.xlabel('Seed Number')
    plt.ylabel('Brier Score')
    count_global_readable = '{:,.2f}'.format(count_global)
    title = f"Average Brier Score on seeds from {os.path.basename(path_folder_fasta)}"
    plt.title(f"{title} \nover {count_global_readable} valid triplets")
    path_image = path_NeighborRes
    plt.legend()
    plt.savefig(f"{path_image}/{title}_{list_list_len_window}.png")
    
    t.stop("Brier Score with contextual Bayes")
    print(Brier_Score_global)
    return Brier_Score_global
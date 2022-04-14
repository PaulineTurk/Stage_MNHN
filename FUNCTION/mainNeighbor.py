from fileNumber import countFile
from timer import Timer
import character as ch
from fastaReader import readFastaMul
import numpy as np
from pathlib import Path
import pandas as pd
import os, shutil



def tripletCount(accession_num, path_folder_pid, file_fasta, pid_inf, triplet_count, delay_num, kp_SeqChoice):    
    # k: known
    # p: predict
    # c: context
    # kp_SeqChoice: choice the reference sequence to look at its neighbor 
    liSeqAliFiltre = readFastaMul(file_fasta)
    pid_couple = np.load(path_folder_pid + "/" + accession_num+ ".pId.npy", allow_pickle='TRUE').item()
    list_AA = ch.characterList()


    if liSeqAliFiltre:
        for name_k, seq_k in liSeqAliFiltre:
            len_seq = len(seq_k)
            for name_p, seq_p in liSeqAliFiltre:
                if name_k != name_p:
                    if pid_couple[name_k][name_p] >= pid_inf:
                        if delay_num < 0:   # before
                            index_range = - delay_num, len_seq
                        elif delay_num > 0: # after
                            index_range = delay_num, len_seq - 1
                        
                        if kp_SeqChoice == "k":
                            seq_c = seq_k
                        elif kp_SeqChoice == "p":
                            seq_c = seq_p

                        for aa_index in range(index_range[0], index_range[1]):       
                            aa_k = seq_k[aa_index] 
                            aa_p = seq_p[aa_index] 
                            if all(x in list_AA for x in [aa_k, aa_p]):
                                index_neighbor = aa_index + delay_num
                                aa_c = seq_c[index_neighbor] 
                                if aa_c in list_AA: 
                                    triplet_count[aa_k][aa_p][aa_c] += 1
    else:
        print(accession_num)
    return triplet_count




def initialisation(list_symbol):
    triplet_count = {}
    for aa_k in list_symbol:
        triplet_count[aa_k] = {}
        for aa_p in list_symbol:
                triplet_count[aa_k][aa_p] = {}
                for aa_c in list_symbol:
                    triplet_count[aa_k][aa_p][aa_c]  = 1   # to avoid issues when the triplet is not in data_train  (Rq. len + 1 pas encore changé)
    return triplet_count




def conditionalProba(list_AA, triplet_count, pseudo_count = 0):
    """pseudo_count = 0, d'Alembert rule not needed because we are working with enough data"""
    intra_couple_count = {}
    intra_couple_number = 0
    len_list = len(list_AA)
    for aa_k in list_AA:
        intra_couple_count[aa_k] = {}
        for aa_c in list_AA: 
            intra_couple_count[aa_k][aa_c] = 0
            for aa_p in list_AA:
                intra_couple_count[aa_k][aa_c] += triplet_count[aa_k][aa_p][aa_c]
                intra_couple_number += triplet_count[aa_k][aa_p][aa_c]
    print("number of internal couples:", '{:,.2f}'.format(intra_couple_number))

    cond_proba = {}
    for aa_k in list_AA:
        cond_proba[aa_k] = {}
        for aa_p in list_AA:
            cond_proba[aa_k][aa_p] = {}
            for aa_c in list_AA: 
                if intra_couple_count[aa_k][aa_c] != 0:
                    cond_proba[aa_k][aa_p][aa_c] = (triplet_count[aa_k][aa_p][aa_c] + pseudo_count) / (intra_couple_count[aa_k][aa_c] + len_list * pseudo_count)  
                else:
                    cond_proba[aa_k][aa_p][aa_c] = 0

    #for aa_k in list_AA:                   # the problem could only happen with very small databases (not the case here) ----- not sure ...
    #    for aa_c in list_AA:
    #        sum_line =sumLine(cond_proba, list_AA, aa_k, aa_c)
    #        if sum_line == 0:
    #            for aa_p in list_AA:
    #                estimation_proba = 0
    #                for aa_c in list_AA:
    #                    estimation_proba += cond_proba[aa_k][aa_p][aa_c]
    #            cond_proba[aa_k][aa_p][aa_c] = estimation_proba

    return cond_proba                 



def sumLine(cond_proba, list_AA, aa_k, aa_c):
    sum_line = 0
    for aa_p in list_AA:
        sum_line += cond_proba[aa_k][aa_p][aa_c]
    return sum_line


def sumPlate(dico_triple):
    for aa_1 in dico_triple:
        sum_plateau = 0
        for aa_2 in dico_triple[aa_1]:
            for aa_3 in dico_triple[aa_1][aa_2]:
                sum_plateau += dico_triple[aa_1][aa_2][aa_3]
        print("{}, {}".format(aa_1, sum_plateau))




def simpleContextualBlosum(folder_fasta, path_folder_pid, path_NeighborResX, delay_num, kp_SeqChoice, 
                           pid_inf = 62, scale_factor = 2):
    t = Timer()
    t.start()

    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    #nbre_file = countFile(path_folder_fasta)
    #countfile = 0

    list_AA =  ch.characterList()
    triplet_count = initialisation(list_AA)

    for file_name_fasta in files_in_path_folder_fasta:
            #countfile += 1
            accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
            triplet_count = tripletCount(accession_num, path_folder_pid, file_name_fasta, pid_inf, triplet_count, delay_num, kp_SeqChoice)
            #print(100*countfile/nbre_file)
    t.stop("Compute the conditional probability matrix")

    print("({},{}) - conditional proba".format(delay_num, kp_SeqChoice))
    cond_proba = conditionalProba(list_AA, triplet_count)
    sumPlate(cond_proba)
    path_proba_cond = path_NeighborResX + "/proba_cond_(" + str(delay_num) + " , " + kp_SeqChoice + ")"
    np.save(path_proba_cond, cond_proba) 
    path_proba_cond = path_proba_cond + ".npy"

    return path_proba_cond
    




if __name__ == '__main__': 

    # /Users/pauline/Desktop/data/NeighborRes    folder to create
    path_folder_pid = "/Users/pauline/Desktop/data/PID_couple"
    list_percentage = [0.05, 0.5]
    for percentage in list_percentage:
        folder_fasta = "/Users/pauline/Desktop/data/PfamSplit_" + str(percentage) +"/PfamTrain"
        path_NeighborResX = "/Users/pauline/Desktop/data/NeighborRes/NeighborRes_" + str(percentage)  # folder to create
        for delay_num in [-1, 1]:
            for kp_SeqChoice in ["k", "p"]:
                print(percentage)
                print("{}, {}".format(delay_num, kp_SeqChoice))
                path_proba_cond = simpleContextualBlosum(folder_fasta, path_folder_pid, path_NeighborResX, delay_num, kp_SeqChoice, pid_inf = 62, scale_factor = 2)
  
                # visual check
                cond_proba = np.load(path_proba_cond, allow_pickle='TRUE').item()
                df_cond_proba = np.transpose(pd.DataFrame.from_dict(cond_proba))
                print(df_cond_proba)  

 
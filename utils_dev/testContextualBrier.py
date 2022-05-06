from ast import If
import os, random
import fastaReader
from pathlib import Path
import timer
import numpy as np



# regrouper les couples possibles....+ context enlarge

# for the convertion !!!!!
#https://www.projectpro.io/recipes/convert-dictionary-matrix-or-narray-in-python


def seedSelection(path_folder_fasta):

    random_seed = random.choice(os.listdir((Path(path_folder_fasta))))
    path_random_seed = f"{path_folder_fasta}/{random_seed}"
    len_random_seed = len(fastaReader.readFastaMul(path_random_seed))

    while len_random_seed<=1: # il aurait fallu les enlever dès le début car on sait qu'ils ne seront pas informatifs (à enlever à terme),
        path_random_seed, len_random_seed = seedSelection(path_folder_fasta)

    return path_random_seed, len_random_seed 



# ajouter un doc qui pour chaque pid_couple contient son pid max excluant la diagonale pour savoir si le seed peut etre utilisé

def seedSeqSelection(path_folder_fasta, path_pid_folder):  
    path_random_seed, len_random_seed = seedSelection(path_folder_fasta)

    # load pid file of the seed 
    accession_num = os.path.basename(path_random_seed).split(".")[0] + '.' + os.path.basename(path_random_seed).split(".")[1]
    pid_couple = np.load(f"{path_pid_folder}/{accession_num}.pId.npy", allow_pickle='TRUE').item()
    # test on pid_max



    values_pid_couple = pid_couple.values()
    print(values_pid_couple)



    return path_random_seed, len_random_seed  



def seqCoupleSelection(path_random_seed, len_random_seed, path_pid_folder, pid_inf = 62):

    seed = fastaReader.readFastaMul(path_random_seed)
    (name_1, seq_1), (name_2, seq_2)= random.sample(seed, 2)

    # load pid file of the seed 
    accession_num = os.path.basename(path_random_seed).split(".")[0] + '.' + os.path.basename(path_random_seed).split(".")[1]
    pid_couple = np.load(f"{path_pid_folder}/{accession_num}.pId.npy", allow_pickle='TRUE').item()

    pid_selected_couple = pid_couple[name_1][name_2]

    while pid_selected_couple < pid_inf:
        (name_1, seq_1), (name_2, seq_2)= random.sample(seed, 2)
        pid_selected_couple = pid_couple[name_1][name_2]
        #print(pid_selected_couple)
    return (name_1, seq_1), (name_2, seq_2), pid_selected_couple




if __name__ == '__main__': 
    path_folder_fasta = "/Users/pauline/Desktop/Overfitting_test/Pfam_fasta_99"
    path_pid_folder = "/Users/pauline/Desktop/Overfitting_test/PID_couple"

    for i in range(10):
        path_random_seed, len_random_seed  = seedSeqSelection(path_folder_fasta, path_pid_folder)
        print(f"seed selected, nbre of seq: {path_random_seed}, {len_random_seed}")
        #(name_1, seq_1), (name_2, seq_2), pid_selected_couple = seqCoupleSelection(path_random_seed, len_random_seed, path_pid_folder)
        #print("pid couple:", pid_selected_couple)
        #print("seq_1:", (name_1, seq_1))
        #print("seq_2:", (name_2, seq_2))




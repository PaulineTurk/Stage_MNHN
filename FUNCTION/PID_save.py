import os, shutil
from pathlib import Path
import FUNCTION.PID as PID
import numpy as np
import pandas as pd
from FUNCTION.timer import Timer
from FUNCTION.fastaReader import readFastaMul


def pIdCouple(path_fasta_file, path_file_pId):
    liste_seq = readFastaMul(path_fasta_file)
    pid_couple = {}
    for name_1, seq_1 in liste_seq:
        pid_couple[name_1] = {}
        for name_2, seq_2 in liste_seq:
            pid_couple[name_1][name_2] = PID.pId(seq_1, seq_2)
    np.save(path_file_pId, pid_couple) 



def savePId(path_folder_fasta, path_folder_pId):
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_pId):
        shutil.rmtree(path_folder_pId) 
    os.mkdir(path_folder_pId)
    
    path_folder_fasta = Path(path_folder_fasta)
    path_fasta_files = path_folder_fasta.iterdir()

    for path_fasta_file in path_fasta_files:
        num_accession = os.path.basename(path_fasta_file).split(".")[0] + '.' + os.path.basename(path_fasta_file).split(".")[1]
        path_file_pId = path_folder_pId + "/" + num_accession + '.pId'
        pIdCouple(path_fasta_file, path_file_pId)
    t.stop("Compute and save the pId files")




if __name__ == '__main__': 
    path_folder_fasta = "/Users/pauline/Desktop/data/Pfam_fasta"
    path_folder_pId = "/Users/pauline/Desktop/data/PID_couple"
    savePId(path_folder_fasta, path_folder_pId)     #    14629.03248 s    
    # 
    #pid_check = np.load("/Users/pauline/Desktop/data/Pfam_test_PID/PF00002.27.pId.npy" ,allow_pickle='TRUE').item()  
    #df_pid_check = np.transpose(pd.DataFrame.from_dict(pid_check))
    #print(df_pid_check)  
 




    # MINI TEST
    #path_folder_fasta = "Pfam_test"
    #path_folder_pId = "Pfam_test_PID"
    #savePId(path_folder_fasta, path_folder_pId)        
    # 
    #pid_check = np.load("Pfam_test_PID/PF00002.27.pId.npy" ,allow_pickle='TRUE').item()  
    #df_pid_check = np.transpose(pd.DataFrame.from_dict(pid_check))
    #print(df_pid_check)
        
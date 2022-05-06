import numpy as np
import os, shutil
from pathlib import Path
import fastaReader
from timer import Timer

# give same proba to each valid couple to be selected (1 file par valid couple ??? comme ca revient à faire des seed à 2seq,
# reste à en choisir un à chaque fois) et dans le nom du fichier mettre l'accession number



def writeValidCouple(path_seed, path_pid_folder, path_couple_folder, pid_inf = 62):
    """
    path_seed: taken from Pfam_B (fraction of Pfam from which the test exemples are taken)
    """

    # load pid file of the seed + add accession_num to the name
    accession_num = os.path.basename(path_seed).split(".")[0] + '.' + os.path.basename(path_seed).split(".")[1]
    pid_couple = np.load(f"{path_pid_folder}/{accession_num}.pId.npy", allow_pickle='TRUE').item()
    
    seed = fastaReader.readFastaMul(path_seed)
    nbre_Seq = len(seed)

    couple_count = 0
    for i in range(nbre_Seq):
        name_1 = seed[i][0]
        for j in range(i+1, nbre_Seq):
            name_2 = seed[j][0]

            if pid_couple[name_1][name_2] >= pid_inf:
                #print(f"{name_1}, {name_2}, {pid_couple[name_1][name_2]}")  #pid couple check
                couple_count += 1
                path_file = f"{path_couple_folder}/{accession_num}_{couple_count}"

                input_handle = open(path_seed)
                output_handle = open(path_file, "w")

                for l in input_handle:
                    if l[0] == ">":  
                        if l == f">{name_1}\n" or l == f">{name_2}\n":
                            write = True
                        else:
                            write = False
                    if write == True:
                        output_handle.write(l)

                output_handle.close()
                input_handle.close()









def writeAllValidCouple(path_folder_fasta, path_couple_folder, path_pid_folder, pid_inf = 62):
    t = Timer()
    t.start()

    if not os.path.exists(path_couple_folder):
        os.mkdir(path_couple_folder)

    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()

    for path_seed in files_in_path_folder_fasta:
            writeValidCouple(path_seed, path_pid_folder, path_couple_folder, pid_inf)

    t.stop("Compute the conditional probability matrix")









def countValidCouple(path_seed, path_pid_folder, list_accession_num, list_couple_nbre, pid_inf):
  
    # load pid file of the seed + add accession_num to the name
    accession_num = os.path.basename(path_seed).split(".")[0] + '.' + os.path.basename(path_seed).split(".")[1]
    list_accession_num.append(accession_num)
    print(accession_num)
    pid_couple = np.load(f"{path_pid_folder}/{accession_num}.pId.npy", allow_pickle='TRUE').item()
    nbre_Seq = len(path_seed)
    #print(nbre_Seq)

    seed = fastaReader.readFastaMul(path_seed)

    couple_count = 0
    for i in range(nbre_Seq):
        name_1 = seed[i][0]
        #print("i", i)
        for j in range(i+1, nbre_Seq):
            #print("j", j)
              
            name_2 = seed[j][0]
            if pid_couple[name_1][name_2] >= pid_inf:
                couple_count += 1
    list_couple_nbre.append(couple_count)
    return list_accession_num, list_couple_nbre


def countAllValidCouple(path_folder_fasta, path_pid_folder, pid_inf = 62):
    t = Timer()
    t.start()

    list_accession_num = []
    list_couple_nbre = []

    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()

    for path_seed in files_in_path_folder_fasta:
        print(path_seed)
        list_accession_num, list_couple_nbre = countValidCouple(path_seed, path_pid_folder, list_accession_num, list_couple_nbre, pid_inf)
    return list_accession_num, list_couple_nbre






if __name__ == "__main__":
    path_folder_fasta = "/Users/pauline/Desktop/test_dev_v2/test_1/PfamSplit_50/Pfam_B"  # test on mini Pfam_B
    #path_seed = "/Users/pauline/Desktop/test_dev_v2/test_1/PfamSplit_50/Pfam_B/PF00002.27.B" # test on 1 seed from mini Pfam_B
    path_pid_folder = "/Users/pauline/Desktop/test_dev_v2/PID_couple"
    path_couple_folder = "/Users/pauline/Desktop/test_dev_v2/test_1/testCoupleSelection"


    writeAllValidCouple(path_folder_fasta, path_couple_folder, path_pid_folder, pid_inf = 62)





    #prob out of range ...

    #list_accession_num, list_couple_nbre = countAllValidCouple(path_folder_fasta, path_pid_folder)
    #print(list_couple_nbre)
    #list_accession_num, list_couple_nbre  = [], []
    #path_seed = "/Users/pauline/Desktop/Stage_MNHN/PF07966.15.B"
    #pid_inf = 62

    #countValidCouple(path_seed, path_pid_folder, list_accession_num, list_couple_nbre, pid_inf)






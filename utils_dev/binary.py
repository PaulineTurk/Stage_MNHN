import struct
import sys
import os
import numpy as np
import fastaReader
from timer import Timer
from pathlib import Path






def writeValidCoupleBinary(path_seed, path_pid_folder, path_new_file, pid_inf):
     # load pid file of the seed + add accession_num to the name
     accession_num = os.path.basename(path_seed).split(".")[0] + '.' + os.path.basename(path_seed).split(".")[1]
     pid_couple = np.load(f"{path_pid_folder}/{accession_num}.pId.npy", allow_pickle='TRUE').item()

     seed = fastaReader.readFastaMul(path_seed)
     nbre_Seq = len(seed)
     count = 0

     for i in range(nbre_Seq):
          name_1, seq_1 = seed[i]
          for j in range(i+1, nbre_Seq):
               name_2, seq_2 = seed[j]
              
               if pid_couple[name_1][name_2] >= pid_inf:

                    count += 1

                    # écriture
                    with open(path_new_file, "ab") as fb:   # path file binary
                         #print(f"{name_1}\n,{name_2}\n,{pid_couple[name_1][name_2]}\n")
                         # a: append (create the file if it doesn't exist, else append to what is already written)
                         # b: binary

                         # ecriture du nombre de caractères encodés
                         fb.write(struct.pack("i", len(name_1)))  # i: integer
                         fb.write(struct.pack("i", len(name_2)))

                         fb.write(struct.pack("i", len(seq_1)))
                         fb.write(struct.pack("i", len(seq_2)))
  
                         # il faut convertir les caractères en bytes
                         name_1_ascii = name_1.encode("ascii")
                         name_2_ascii = name_2.encode("ascii")

                         seq_1_ascii = seq_1.encode("ascii")
                         seq_2_ascii = seq_2.encode("ascii")

                         # ecriture en indiquant le format "s" précédé du nombre de caractères encodés
                         monFormat=str(len(name_1))+"s"
                         fb.write(struct.pack(monFormat, name_1_ascii))
                         monFormat=str(len(name_2))+"s"
                         fb.write(struct.pack(monFormat, name_2_ascii))

                         monFormat=str(len(seq_1))+"s"
                         fb.write(struct.pack(monFormat, seq_1_ascii))
                         monFormat=str(len(seq_2))+"s"
                         fb.write(struct.pack(monFormat, seq_2_ascii))
     print(f"{accession_num}: nbre de seq, nbre de couples valides {nbre_Seq}, {count}")
     
                    

def writeAllValidCoupleBinary(path_folder_seed, path_pid_folder, path_new_file, pid_inf = 62):
     t = Timer()
     t.start()

     files_in_path_folder_seed = Path(path_folder_seed).iterdir()

     for path_seed in files_in_path_folder_seed:
          writeValidCoupleBinary(path_seed, path_pid_folder, path_new_file, pid_inf)

     t.stop("Compute binary valid couples over Pfam_test")
     


def readValidCouplebinary(path_seed_binary):
     # path_seed_binary: equivalent to previous path_new_file
     pass


if __name__ == '__main__': 
    path_folder_seed = "/Users/pauline/Desktop/Stage_MNHN/test_write_binary"  # test on mini Pfam_B
    #path_seed = "/Users/pauline/Desktop/test_dev_v2/test_1/PfamSplit_50/Pfam_B/PF00002.27.B" # test on 1 seed from mini Pfam_B
    path_pid_folder = "/Users/pauline/Desktop/test_dev_v2/PID_couple"
    #path_new_file = "/Users/pauline/Desktop/test_dev_v2/test_1/testCoupleSelection"
    #path_new_file = "/Users/pauline/Desktop/test_dev_v2/testCoupleSelection.bin" # appear when searched from a terminal
    path_new_file = "testCoupleSelection.bin" 


    writeAllValidCoupleBinary(path_folder_seed, path_pid_folder, path_new_file, pid_inf = 62)

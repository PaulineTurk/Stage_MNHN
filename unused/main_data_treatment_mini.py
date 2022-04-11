import FUNCTION.Stockholm_Fasta as SF
import FUNCTION.fileNumber as cfn
import FUNCTION.PID_save as pid
import FUNCTION.redundancy
import FUNCTION.trimming
import FUNCTION.character as ch
import numpy as np
import pandas as pd

# /Users/pauline/Desktop/MINI_TEST à creer en amont et contenant pFAM_A.Seed

# separation of the Stockholm file into Stockholm files
file_name = "/Users/pauline/Desktop/MINI_TEST/pFAM_A.Seed"
folder_name = "/Users/pauline/Desktop/MINI_TEST/pFAM_Stockholm"
SF.separationStockholm(file_name, folder_name)  #  0.00653 s
cfn.countFile(folder_name)  # 4 files


# conversion from stockholm into fasta files
folder_fasta = "/Users/pauline/Desktop/MINI_TEST/pFAM_fasta"
folder_stockholm = "/Users/pauline/Desktop/MINI_TEST/pFAM_Stockholm"
SF.multiStockholmToFasta(folder_fasta, folder_stockholm)  # 0.07315 s
cfn.countFile(folder_fasta)   # 4 files


# pid couple calculation
path_folder_fasta = "/Users/pauline/Desktop/MINI_TEST/pFAM_fasta"
path_folder_pId = "/Users/pauline/Desktop/MINI_TEST/pid_couple"
pid.savePId(path_folder_fasta, path_folder_pId)     #   0.86083 s   

# checking the visualisation
#pid_check = np.load("/Users/pauline/Desktop/MINI_TEST/pid_couple/PF00244.23.pId.npy" ,allow_pickle='TRUE').item()  
#df_pid_check = np.transpose(pd.DataFrame.from_dict(pid_check))
#print(df_pid_check)  



# dealing with the redundancy issue
folder_fasta = "/Users/pauline/Desktop/MINI_TEST/pFAM_fasta"
folder_fasta_non_redondant =  "/Users/pauline/Desktop/MINI_TEST/pFAM_fasta_99"
pid_sup = 99
forbidden_symbol = ["B", "Z", "X", "-"]
FUNCTION.redundancy.savePIdNonRedondant(folder_fasta, folder_fasta_non_redondant, pid_sup, forbidden_symbol)  # 0.62418 s


# data trimming 
path_folder_fasta =  "/Users/pauline/Desktop/MINI_TEST/pFAM_fasta_99"   # the non-redundant one
path_folder_fasta_trimmed = "/Users/pauline/Desktop/MINI_TEST/pFAM_fasta_99_trimmed"
list_AA = ch.characterList()
exclusion_fraction = 0.5
trim_fraction_tolerated = 0.5
FUNCTION.trimming.trimSave(path_folder_fasta, path_folder_fasta_trimmed, list_AA, exclusion_fraction, trim_fraction_tolerated) # 207.76572 s
    

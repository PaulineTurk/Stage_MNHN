import FUNCTION.Stockholm_Fasta as SF
import FUNCTION.fileNumber as cfn
import FUNCTION.PID_save as pid
import FUNCTION.redundancy
import FUNCTION.trimming
import FUNCTION.character as ch
import numpy as np
import pandas as pd



# separation of the Stockholm file into Stockholm files
file_name = "/Users/pauline/Desktop/data/Pfam-A.seed"
folder_name = "/Users/pauline/Desktop/data/Pfam_Stockholm"
SF.separationStockholm(file_name, folder_name) 
cfn.countFile(folder_name)  


# conversion from stockholm into fasta files
folder_fasta = "/Users/pauline/Desktop/data/Pfam_fasta"
folder_stockholm = "/Users/pauline/Desktop/data/Pfam_Stockholm"
SF.multiStockholmToFasta(folder_fasta, folder_stockholm) 
cfn.countFile(folder_fasta)  


# pid couple calculation
path_folder_fasta =  "/Users/pauline/Desktop/data/Pfam_fasta"
path_folder_pId =  "/Users/pauline/Desktop/data/PID_couple"
pid.savePId(path_folder_fasta, path_folder_pId)   

# checking the visualisation
#pid_check = np.load("/Users/pauline/Desktop/MINI_TEST/pid_couple/PF00244.23.pId.npy" ,allow_pickle='TRUE').item()  
#df_pid_check = np.transpose(pd.DataFrame.from_dict(pid_check))
#print(df_pid_check)  



# dealing with the redundancy issue
folder_fasta =  "/Users/pauline/Desktop/data/Pfam_fasta"
folder_fasta_non_redondant =   "/Users/pauline/Desktop/data/Pfam_fasta_99"
pid_sup = 99
forbidden_symbol = ["B", "Z", "X", "-"]   # ATTENTION: il me semble qu'il y a des symboles etranges dans Pfam ex. O et U ...
FUNCTION.redundancy.savePIdNonRedondant(folder_fasta, folder_fasta_non_redondant, pid_sup, forbidden_symbol)  


# data trimming 
path_folder_fasta =  "/Users/pauline/Desktop/data/Pfam_fasta_99"   # the non-redundant one
path_folder_fasta_trimmed = "/Users/pauline/Desktop/data/Pfam_fasta_99_trimmed"
list_AA = ch.characterList()
exclusion_fraction = 0.5
trim_fraction_tolerated = 0.5
FUNCTION.trimming.trimSave(path_folder_fasta, path_folder_fasta_trimmed, list_AA, exclusion_fraction, trim_fraction_tolerated) # 207.76572 s
    

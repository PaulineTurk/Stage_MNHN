import os, shutil
from pickle import FALSE
from pathlib import Path
import FUNCTION.PID as percentageID 
from FUNCTION.timer import Timer
from FUNCTION.fastaReader import readFastaMul



def nonRedundant(file_fasta, pid_sup, file_seq_non_redundant, forbidden_symbol):
    liste_seq = readFastaMul(file_fasta)
    cluster = percentageID.clusterAntiRedundancy(liste_seq, pid_sup, file_seq_non_redundant, forbidden_symbol)
    seq_non_redundant = percentageID.representativeNonRedundant(cluster)

    input_handle = open(file_fasta)
    output_handle = open(file_seq_non_redundant, "w")

    flag_write = False
    for l in input_handle:
        if l[0] == ">":   
            if l[1:-1].split(" ")[0] in seq_non_redundant:    # keep the name only 
                flag_write = True
            else:
                flag_write = False
        if flag_write == True:
            output_handle.write(l)
    output_handle.close()
    input_handle.close()



def savePIdNonRedondant(folder_fasta, folder_fasta_non_redondant, pid_sup, forbidden_symbol = ["B", "Z", "X", "-"]):
    t = Timer()
    t.start()
    path_folder_fasta_non_redondant = folder_fasta_non_redondant + '/'
    if os.path.isdir(path_folder_fasta_non_redondant):
        shutil.rmtree(path_folder_fasta_non_redondant) 
    os.mkdir(path_folder_fasta_non_redondant)
    
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        file_name_fasta_non_redondant = accession_num + '.fasta.non.redundant'
        nonRedundant(file_name_fasta, pid_sup, path_folder_fasta_non_redondant + file_name_fasta_non_redondant, forbidden_symbol)
    t.stop("Compute and save non-redundant files")




if __name__ == '__main__': 
    #savePIdNonRedondant("Pfam_fasta_trimAl_gappyout_header_corrected", "Pfam_fasta_trimAl_gappyout_header_corrected_non_redundant", 99) #3090.13537 s
    #savePIdNonRedondant("Pfam_fasta_trimAl_nogaps_header_corrected", "Pfam_fasta_trimAl_nogaps_header_corrected_non_redundant", 99)   # 1568.50542 s  PF07703.17.fasta.non.redundant seul à seq vide --> fichier vide
    

    ### mini test
    savePIdNonRedondant("Pfam_test", "Pfam_test_99", 99)
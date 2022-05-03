import os, shutil
from pickle import FALSE
from pathlib import Path
from timer import Timer



def LowerToUpper(file_fasta, file_seq_upper):
    input_handle = open(file_fasta)
    output_handle = open(file_seq_upper, "w")

    for l in input_handle:
        if l[0] != ">":   
            l = l.upper()
        output_handle.write(l)
    output_handle.close()
    input_handle.close()



def multiLowerToUpper(path_folder_fasta, path_folder_fasta_upper):
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_fasta_upper):
        shutil.rmtree(path_folder_fasta_upper) 
    os.mkdir(path_folder_fasta_upper)
    
    path_folder_fasta = Path(path_folder_fasta)
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num_part_1 = os.path.basename(file_name_fasta).split(".")[0]
        accession_num_part_2 = os.path.basename(file_name_fasta).split(".")[1]
        accession_num = f"{accession_num_part_1}.{accession_num_part_2}"
        path_file_name_fasta_upper = f"{path_folder_fasta_upper }/{accession_num}.fasta.upper"
        LowerToUpper(file_name_fasta, path_file_name_fasta_upper)
    t.stop("Correction upper files")



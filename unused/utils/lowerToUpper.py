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



def multiLowerToUpper(folder_fasta, folder_fasta_upper):
    t = Timer()
    t.start()
    path_folder_fasta_upper = folder_fasta_upper + '/'
    if os.path.isdir(path_folder_fasta_upper ):
        shutil.rmtree(path_folder_fasta_upper ) 
    os.mkdir(path_folder_fasta_upper)
    
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    for file_name_fasta in files_in_path_folder_fasta:
        accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        file_name_fasta_upper = accession_num + '.fasta.upper'
        LowerToUpper(file_name_fasta, path_folder_fasta_upper  + file_name_fasta_upper)
    t.stop("Correction upper files")



if __name__ == '__main__': 
    multiLowerToUpper("/Users/pauline/Desktop/test_upper/Pfam_fasta","/Users/pauline/Desktop/test_upper/Pfam_fasta_upper" )  # 20s sur Pfam

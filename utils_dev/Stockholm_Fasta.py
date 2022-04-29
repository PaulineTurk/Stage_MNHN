from Bio import AlignIO
import os, shutil
from pathlib import Path
from timer import Timer




def separationStockholm(path_file_name, path_folder_to_save):
    """
    Separate a multiStockholm file into monoStockholm files.

    Args:
        path_file_name: the multiStockholm file
        path_folder_to_save: folder where the monoStockholm files generated are saved 
    """
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_to_save):
        shutil.rmtree(path_folder_to_save) 
    os.mkdir(path_folder_to_save)

    # collect the accession number of each alignment
    list_accession_num = []

    with open(path_file_name) as file_name:
        for l in file_name:
            if l[0:7] == "#=GF AC":
                init_accession_num = l.index('PF')
                accession_num = l[init_accession_num: -1]
                list_accession_num.append(accession_num)

    nbre_file = len(list_accession_num)

    # generate the monoStockholm files that are named after the accession number's alignment ID
    file_out_nbre = 1
    with open(path_file_name) as file_name:
        path_file_out = f"{path_folder_to_save}/{list_accession_num[file_out_nbre]}.stockholm"
        with open(path_file_out) as file_out:
            for l in file_name:
                while l[0:2] != "//" and file_out_nbre <= nbre_file - 1: # avoid generating an empty file at the end:
                    file_out.write(l)
                else:
                    file_out_nbre += 1
                    path_file_out = f"{path_folder_to_save}/{list_accession_num[file_out_nbre]}.stockholm"
        print("i m lost")
                    

    
    t.stop("Separation of the multiStockholm file into monoStockholm files")







def stockholmToFasta(file_name_fasta, file_name_stockholm) :
    """Convert a Stockholm file into Fasta file"""

    input_handle = open(file_name_stockholm)
    output_handle = open(file_name_fasta, "w")
    alignments = AlignIO.parse(input_handle,  "stockholm")
    for alignment in alignments:
        AlignIO.write([alignment], output_handle, "fasta")
    output_handle.close()
    input_handle.close()






def multiStockholmToFasta(folder_fasta, folder_stockholm):
    """Convert Stockholm files into Fasta files"""
    t = Timer()
    t.start()

    path_folder_fasta = folder_fasta + '/'
    if os.path.isdir(path_folder_fasta):
        shutil.rmtree(path_folder_fasta) 
    os.mkdir(path_folder_fasta)
    
    path_folder_stockholm = Path(folder_stockholm + '/')
    files_in_path_folder_stockholm = path_folder_stockholm.iterdir()

    for file_name_stockholm in files_in_path_folder_stockholm:
        accession_num = file_name_stockholm.name.split(".")[0] + '.' + file_name_stockholm.name.split(".")[1]
        file_name_fasta = accession_num + '.fasta'
        stockholmToFasta(path_folder_fasta + file_name_fasta, file_name_stockholm)
    t.stop("Conversion of Stockholm files into Fasta files")







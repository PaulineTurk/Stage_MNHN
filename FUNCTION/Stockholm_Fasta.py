from Bio import AlignIO
import os, shutil
from pathlib import Path
from timer import Timer
import fileNumber as cfn




def separationStockholm(file_name, folder_name):
    """
    Separate a multiStockholm file into monoStockholm files.

    Args:
        file_name: name of the multiStockholm file downloaded from Pfam
        folder_name: folder name where the monoStockholm files generated are saved 
    """
    t = Timer()
    t.start()

    input_handle = open(file_name)

    path_folder = folder_name + '/'
    if os.path.isdir(path_folder):
        shutil.rmtree(path_folder) 
    os.mkdir(path_folder)

    # collect the accession number of each alignment
    list_accession_num = []
    for l in input_handle:
        if l[0:7] == "#=GF AC":
            init_accession_num = l.index('PF')
            accession_num = l[init_accession_num: -1]
            list_accession_num.append(accession_num)
    input_handle.close()
    nbre_file = len(list_accession_num)

    # generate the monoStockholm files named after the accession number's alignment
    input_handle = open(file_name)
    file_out_nbre = 0
    path_file_out = path_folder + list_accession_num[file_out_nbre] + ".stockholm"
    output_handle = open(path_file_out, "w")

    for l in input_handle:
        output_handle.write(l)
        if l[0:2] == "//" and file_out_nbre <= nbre_file - 2: # avoid generating an empty file at the end
            output_handle.close()
            file_out_nbre += 1
            path_file_out = path_folder + list_accession_num[file_out_nbre] + ".stockholm"
            output_handle = open(path_file_out, "w")

    output_handle.close()
    input_handle.close()
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






if __name__ == '__main__': 
    # Pfam 35.0 (November 2021, 19632 entries)  
    # http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
    # Pfam-A.seed.gz,  15-Nov-2021 12:23,    153197377   (file downloaded in march 2022)

    separationStockholm("Pfam-A.seed", "Pfam_Stockholm")  #  12.94936 s
    cfn.countFile("Pfam_Stockholm")  # 19632 files

    multiStockholmToFasta("Pfam_fasta", "Pfam_Stockholm")  # 62.50208 s
    cfn.countFile("Pfam_fasta")   # 19632 files




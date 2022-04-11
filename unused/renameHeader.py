import os, shutil
from pathlib import Path
from FUNCTION.timer import Timer


def removePB(header):
    """Remove from the header of an alignment the information 
       on its length automatically added by trimAl"""
    liste_header = header.split(" ")
    new_header = str.join(" ", liste_header[0:-2]) + "\n"
    return new_header




def correctHeader(file_trimAl_original, file_trimAl_header_corrected):
    """Rewrite a file by removing from the header of each alignment the 
       information on its length automatically added by trimAl"""
    input_handle = open(file_trimAl_original)
    output_handle = open(file_trimAl_header_corrected, "w")
    for l in input_handle:
        if l[0] == ">":  
            l = removePB(l)
            output_handle.write(l)
        else:
            output_handle.write(l)
    output_handle.close()
    input_handle.close()



def saveTrimHeaderCorrected(folder_trimAl_original, folder_trimAl_header_corrected):
    """ For each file in folder_trimAl_original, rewrite it by removing from the header 
        of each alignment the information on its length automatically added by trimAl"""
    t = Timer()
    t.start()
    path_folder_trimAl_header_corrected = folder_trimAl_header_corrected+ '/'
    if os.path.isdir(path_folder_trimAl_header_corrected):
        shutil.rmtree(path_folder_trimAl_header_corrected) 
    os.mkdir(path_folder_trimAl_header_corrected)
    
    path_folder_trimAl_original = Path(folder_trimAl_original+ '/')
    files_in_path_trimAl_original = path_folder_trimAl_original.iterdir()

    for file_trimAl_original in files_in_path_trimAl_original:
        num_accession = os.path.basename(file_trimAl_original).split(".")[0] + '.' + os.path.basename(file_trimAl_original).split(".")[1]
        file_trimAl_header_corrected = num_accession + '.fasta.trim'
        correctHeader(file_trimAl_original, path_folder_trimAl_header_corrected + file_trimAl_header_corrected)
    t.stop("Headers correction after trimAl rename process")




if __name__ == '__main__': 
    #saveTrimHeaderCorrected("Pfam_fasta_trimAl_gappyout", "Pfam_fasta_trimAl_gappyout_header_corrected")   # 20.75643 s 
    saveTrimHeaderCorrected("Pfam_fasta_trimAl_nogaps", "Pfam_fasta_trimAl_nogaps_header_corrected")  # 17.33124 s
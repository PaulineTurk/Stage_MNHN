from fastaReader import readFastaMul
from timer import Timer
import character as ch
import os, shutil
from pathlib import Path


def trimDecision(path_fasta_file, exclusion_fraction, trim_fraction_tolerated, list_AA):
    """
    Management of gaps and ambiguous amino-acids:
    - exclusion_fraction: if the fraction of residus at a position in the sequence is > exclusion_fraction 
                            --> remove this position from each sequence of the seed
    - trim_fraction_tolarated: if the fraction of positions removed is > trim_fraction_tolerated
                            --> remove the seed from de dataset
    """
    list_seq = readFastaMul(path_fasta_file)
    nbre_seq = len(list_seq)
    len_seq = len(list_seq[0][1])
    count_excluded = [0]*len_seq

    for _, seq in list_seq:
        for index_aa in range(0, len_seq):
            if seq[index_aa] not in list_AA:
                count_excluded[index_aa] += 1
    fraction_excluded = [count/nbre_seq for count in count_excluded]

    to_exclude = []
    for fraction in fraction_excluded:
        if fraction > exclusion_fraction:
            to_exclude.append(1)  
        else:
            to_exclude.append(0)

    seed_trimed_kept = sum(to_exclude)/len_seq < trim_fraction_tolerated

    return seed_trimed_kept, to_exclude






def trimApplication(path_fasta_file, path_fasta_file_trimmed, seed_trimed_kept, to_exclude, count_seed_excluded):
    
    if seed_trimed_kept == False:
        count_seed_excluded += 1

    else:
        input_handle = open(path_fasta_file)
        output_handle = open(path_fasta_file_trimmed, "w")
        list_seq = readFastaMul(path_fasta_file)
        len_seq = len(list_seq[0][1])

        for l in input_handle:
            if l[0] == ">":  
                output_handle.write(l)
                aa_counter = 0
            else:
                for elem in l[0:-1]: 
                    if to_exclude[aa_counter] == 0:
                        output_handle.write(elem)
                    aa_counter += 1
                if aa_counter == len_seq:      # if the sequence is finished, go back to the line
                    output_handle.write(l[-1])
        output_handle.close()
        input_handle.close()
    return count_seed_excluded




def trimSave(path_folder_fasta, path_folder_fasta_trimmed, list_AA, exclusion_fraction, trim_fraction_tolerated):
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_fasta_trimmed):
        shutil.rmtree(path_folder_fasta_trimmed) 
    os.mkdir(path_folder_fasta_trimmed)

    path_folder_fasta = Path(path_folder_fasta)
    path_fasta_files = path_folder_fasta.iterdir()

    count_seed_excluded = 0

    for path_fasta_file in path_fasta_files:
        seed_trimed_kept, to_exclude = trimDecision(path_fasta_file, exclusion_fraction, trim_fraction_tolerated, list_AA)
        num_accession = os.path.basename(path_fasta_file).split(".")[0] + '.' + os.path.basename(path_fasta_file).split(".")[1]
        name_fasta_file_trimmed = num_accession + '.trim'
        path_fasta_file_trimmed = path_folder_fasta_trimmed + "/"+ name_fasta_file_trimmed
        count_seed_excluded = trimApplication(path_fasta_file, path_fasta_file_trimmed, seed_trimed_kept, to_exclude, count_seed_excluded)
    print("count_seed_excluded:", count_seed_excluded)

    t.stop("Compute and save the trimmed files")



if __name__ == '__main__': 
    #path_folder_fasta = "/Users/pauline/Desktop/data/Pfam_fasta"
    #path_folder_fasta_trimmed = "/Users/pauline/Desktop/data/Pfam_fasta_trimmed"
    #list_AA = ch.characterList()
    #exclusion_fraction = 0.5
    #trim_fraction_tolerated = 0.5
    #trimSave(path_folder_fasta, path_folder_fasta_trimmed, list_AA, exclusion_fraction, trim_fraction_tolerated) # 207.76572 s



    # mini test   
    path_folder_fasta = "Pfam_test"
    path_folder_fasta_trimmed = "Pfam_test_trimmed"
    list_AA = ch.characterList()
    exclusion_fraction = 0.5
    trim_fraction_tolerated = 0.5
    trimSave(path_folder_fasta, path_folder_fasta_trimmed, list_AA, exclusion_fraction, trim_fraction_tolerated) # 207.76572 s
    

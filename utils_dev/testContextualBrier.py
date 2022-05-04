import os, random
import fastaReader
from pathlib import Path



def seedSelection(path_folder_fasta):
    random_seed = random.choice(os.listdir((Path(path_folder_fasta))))
    path_random_seed = f"{path_folder_fasta}/{random_seed}"
    len_random_seed = len(fastaReader.readFastaMul(path_random_seed))
    while len_random_seed <=1:  # il aurait fallu les enlever dès le début car on sait qu'ils ne seront pas informatifs (à enlever à terme)
        print(f"{path_random_seed}, {len_random_seed}")
        random_seed = random.choice(os.listdir((Path(path_folder_fasta))))
        path_random_seed = f"{path_folder_fasta}/{random_seed}"
        print(f"{path_random_seed}, {len_random_seed}")
        print(" ")
    return path_random_seed
    





if __name__ == '__main__': 
    path_folder_fasta = "/Users/pauline/Desktop/Overfitting_test/Pfam_fasta_99"
    for i in range(10):
        selected_seed = seedSelection(path_folder_fasta)
        print(selected_seed)


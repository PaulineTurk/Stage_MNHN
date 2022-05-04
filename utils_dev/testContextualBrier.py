import os, random
import fastaReader




def seedSelection(path_folder_fasta):
    path_random_seed = random.choice(os.listdir((path_folder_fasta))
    while len(fastaReader.readFastaMul(path_random_seed))<=1:  # il aurait fallu les enlever dès le début car on sait qu'ils ne seront pas informatifs

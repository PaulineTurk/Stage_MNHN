import os, shutil
from pathlib import Path
import numpy as np
import pandas as pd
from timer import Timer
from fastaReader import readFastaMul





def pId(seq_1, seq_2):  
    """Return the percentage of identity between the raw sequences seq_1 and seq_2"""
    pid = 0
    len_seq = len(seq_1)
    for indice_aa in range(len_seq):
        if seq_1[indice_aa] == seq_2[indice_aa]:
            pid += 1
    return 100*pid/len_seq




def pIdCoupleSeq(liste_seq):
    """Return the percentage of identity between each couple of raw sequences in liste_seq"""
    pid_couple = {}
    for seq_1 in liste_seq:
        pid_couple[seq_1] = {}
        for seq_2 in liste_seq:
            pid_couple[seq_1][seq_2] = pId(seq_1, seq_2)
    return pid_couple




def lenSeqReal(Seq, included_residue):
    """Return the lenght of Seq considering only residues from the included_residue list"""
    len_seq_real = 0
    for aa in Seq:
        if aa in included_residue:
            len_seq_real += 1
    return len_seq_real
 


def clusterAntiRedundancy(liste_seq, pid_sup, file_seq_non_redondant, included_residue):    
    """Return a partition of liste_seq of sequences with a percentage of identity greater or equal than pid_sup """
    cluster = {}
    if liste_seq:   # if the list is not empty
        name_0, seq_0 = liste_seq[0] 
        len_seq_real_0 = lenSeqReal(seq_0, included_residue)
        cluster[0] = [(name_0, seq_0, len_seq_real_0)]

        for name_1, seq_1 in liste_seq:
            len_seq_real_1 = lenSeqReal(seq_1, included_residue)
            group = 0
            indice = 0

            while group <= len(cluster) - 1 and indice <= len(cluster[group]) - 1:
                seq_2 = cluster[group][indice][1]
                pourcentage_id = pId(seq_1, seq_2) 
                if pourcentage_id < pid_sup:
                    group += 1
                    indice = 0
                else:
                    if indice == len(cluster[group]) - 1:
                        cluster[group].append((name_1, seq_1, len_seq_real_1))
                        indice += 2 # avoid infinite loop
                    else:
                        indice += 1
            if group == len(cluster):
                cluster[group] = [(name_1, seq_1, len_seq_real_1)]
    else:
        print(file_seq_non_redondant)
    return cluster



def representativeNonRedundant(cluster):
    """Select the first sequence with the longest length in the cluster as the cluster representative"""
    seq_non_redundant = []
    for group in cluster:
        current_group = cluster[group]
        representative = current_group[0]   
        for elem in current_group:
            if elem[2] > representative[2]: # 2 stands for the lenght of a sequence
                                            # wihtout the forbidden symbols
                representative = elem
        seq_non_redundant.append(representative[0])   # 0 stands for the name of the sequence
    return seq_non_redundant






def pIdCouple(path_fasta_file, path_file_pId):
    liste_seq = readFastaMul(path_fasta_file)
    pid_couple = {}
    for name_1, seq_1 in liste_seq:
        pid_couple[name_1] = {}
        for name_2, seq_2 in liste_seq:
            pid_couple[name_1][name_2] = pId(seq_1, seq_2)
    np.save(path_file_pId, pid_couple) 




def savePId(path_folder_fasta, path_folder_pId):
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_pId):
        shutil.rmtree(path_folder_pId) 
    os.mkdir(path_folder_pId)
    
    path_folder_fasta = Path(path_folder_fasta)
    path_fasta_files = path_folder_fasta.iterdir()

    for path_fasta_file in path_fasta_files:
        num_accession = os.path.basename(path_fasta_file).split(".")[0] + '.' + os.path.basename(path_fasta_file).split(".")[1]
        path_file_pId = f"{path_folder_pId}/{num_accession}.pId"
        pIdCouple(path_fasta_file, path_file_pId)
        
    t.stop("Compute and save the pId files")
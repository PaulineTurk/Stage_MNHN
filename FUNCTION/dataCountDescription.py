from pathlib import Path
import pandas as pd
from timer import Timer
import character as ch
from fastaReader import readFastaMul



def dataCountDescription(folder_fasta):
    t = Timer()
    t.start()
    path_folder_fasta = Path(folder_fasta + '/')
    files_in_path_folder_fasta = path_folder_fasta.iterdir()
    list_AA =  ch.characterList()

    count_couple_context = {}
    for aa_1 in list_AA:
        count_couple_context[aa_1] = {}
        for aa_2 in list_AA:
            count_couple_context[aa_1][aa_2] = 0

    nbre_seed = 0
    nbre_seq = 0
    total_position = 0
    total_residu = 0
    residu_count_distribution = {}   # consider all residus (dico construction along the way)



    for file_name_fasta in files_in_path_folder_fasta:
        #accession_num = os.path.basename(file_name_fasta).split(".")[0] + '.' + os.path.basename(file_name_fasta).split(".")[1]
        data_Pfam = readFastaMul(file_name_fasta)
        nbre_seed += 1
        len_seq = len(data_Pfam[0][1])
        total_position += len_seq
        for name, seq in data_Pfam:
            nbre_seq += 1
            total_residu += len_seq 
            for aa_index in range(len_seq - 1):
                if seq[aa_index] in list_AA and seq[aa_index + 1] in list_AA:
                    count_couple_context[seq[aa_index]][seq[aa_index + 1]] += 1 
            for aa in seq:
                if aa in residu_count_distribution:
                    residu_count_distribution[aa] += 1
                else:
                    residu_count_distribution[aa] = 1
    t.stop("Count couple context")


    print("nbre_seed:", '{:,.2f}'.format(nbre_seed))
    print("nbre_seq:", '{:,.2f}'.format(nbre_seq))
    print("nbre_position:", '{:,.2f}'.format(total_position))
    print("total_residu:", '{:,.2f}'.format(total_residu))

    # mean len seq
    mean_len_seq = round(total_residu/nbre_seq, 2)
    print("mean_len_seq:", '{:,.2f}'.format(mean_len_seq))
    # mean nbre seq /seed
    mean_nbre_seq = round(nbre_seq/nbre_seed, 2)
    print("mean_nbre_seq:", '{:,.2f}'.format(mean_nbre_seq))

    # aa percentage distribution
    residu_percentage_distribution = {k: round(100*v / total_residu, 2) for k, v in residu_count_distribution.items()}
    print("residu_count_distribution:", residu_count_distribution)
    print("residu_percentage_distribution:", residu_percentage_distribution)

    minCount, maxCount = minMaxCount(count_couple_context)
    print("minCount couple aa:", '{:,.2f}'.format(minCount))
    print("maxCount couple aa:", '{:,.2f}'.format(maxCount))
    df_matrixCount = pd.DataFrame.from_dict(count_couple_context)  

    return df_matrixCount

def minMaxCount(double_dico):
    minCount = double_dico["A"]["A"]
    maxCount = double_dico["A"]["A"]

    for aa_1 in double_dico:
        for aa_2 in double_dico:
            if double_dico[aa_1][aa_2] < minCount:
                minCount = double_dico[aa_1][aa_2]
            if double_dico[aa_1][aa_2] > maxCount:
                maxCount = double_dico[aa_1][aa_2]
    return minCount, maxCount




if __name__ == '__main__':  
    #folder_pfam_name = "/Users/pauline/Desktop/data/Pfam_fasta"     
    #folder_pfam_name = "/Users/pauline/Desktop/data/Pfam_fasta_99"   
    #folder_pfam_name = "/Users/pauline/Desktop/data/Pfam_fasta_99_trimmed" 
    folder_pfam_name = "/Users/pauline/Desktop/data/PfamSplit_0.05/PfamTrain"
    #folder_pfam_name = "/Users/pauline/Desktop/data/PfamSplit_0.5/PfamTrain"
    #folder_pfam_name = "/Users/pauline/Desktop/data/PfamSplit_5.0/PfamTrain"
    #folder_pfam_name = "/Users/pauline/Desktop/data/PfamSplit_50.0/PfamTrain"
    print(folder_pfam_name)
    df_matrixCount = dataCountDescription(folder_pfam_name)
    #print(df_matrixCount)

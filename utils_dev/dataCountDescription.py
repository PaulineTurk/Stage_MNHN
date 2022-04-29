from pathlib import Path
import pandas as pd
from timer import Timer
from fastaReader import readFastaMul
import matplotlib.pyplot as plt
import os


def countFile(folder):
    path, dirs, files = next(os.walk(folder))
    file_count = len(files)
    print(file_count)
    return file_count

def dataCountDescription(path_folder_to_describe):
    t = Timer()
    t.start()
    path_folder_fasta = Path(path_folder_to_describe + '/')
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
    #count_file = 0


    for file_name_fasta in files_in_path_folder_fasta:
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
        #count_file += 1
        #percentage_file = 100*count_file/19632
        #print(percentage_file)


    print("nbre_seed:", '{:,.2f}'.format(nbre_seed))
    print("nbre_seq:", '{:,.2f}'.format(nbre_seq))
    print("total_residu:", '{:,.2f}'.format(total_residu))
    print("nbre_position:", '{:,.2f}'.format(total_position))


    # mean len seq
    mean_len_seq = round(total_residu/nbre_seq, 2)
    print("mean_len_seq:", '{:,.2f}'.format(mean_len_seq))
    # mean nbre seq /seed
    mean_nbre_seq = round(nbre_seq/nbre_seed, 2)
    print("mean_nbre_seq:", '{:,.2f}'.format(mean_nbre_seq))

    # aa percentage distribution
    #df_residu_count_distribution =  pd.DataFrame.from_dict(residu_count_distribution, orient='index')
    #print("residu_count_distribution:", df_residu_count_distribution)

    plt.bar(list(residu_count_distribution.keys()), residu_count_distribution.values(), color='g')
    plt.xlabel('Residus')
    plt.ylabel('Number')
    path_image = os.path.dirname(path_folder_to_describe)
    title = 'Residu number in ' + os.path.basename(path_folder_to_describe)
    plt.title(title)
    plt.savefig(path_image + "/" + title)
    plt.close()

    residu_percentage_distribution = {k: round(100*v / total_residu, 2) for k, v in residu_count_distribution.items()}
    #df_residu_percentage_distribution =  pd.DataFrame.from_dict(residu_percentage_distribution, orient='index')
    #print("residu_percentage_distribution:", df_residu_percentage_distribution)


    plt.bar(list(residu_percentage_distribution.keys()), residu_percentage_distribution.values(), color='g')
    plt.xlabel('Residus')
    plt.ylabel('Percentage')
    path_image = os.path.dirname(path_folder_to_describe)
    title = 'Residu percentage in ' + os.path.basename(path_folder_to_describe)
    plt.title(title)
    plt.savefig(path_image + "/" + title)
    plt.close()

    #minCount, maxCount = minMaxCount(count_couple_context)                # pas sure à garder !
    #print("minCount couple aa:", '{:,.2f}'.format(minCount))
    #print("maxCount couple aa:", '{:,.2f}'.format(maxCount))
    #df_matrixCount = pd.DataFrame.from_dict(count_couple_context)  
    
    #return df_matrixCount
    t.stop("Time for data description")




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
    #list_percentage = [0.05, 0.5, 5, 50]
    #for percentage in list_percentage:
    #    folder_pfam_name = "/Users/pauline/Desktop/data/PfamSplit_" + str(percentage) + "/PfamTrain"
    #    print(folder_pfam_name)
    #    df_matrixCount = dataCountDescription(folder_pfam_name)
    #    #print(df_matrixCount)



    folder_pfam_name = "/Users/pauline/Desktop/test_upper/Pfam_fasta_upper"   # 286.56985 s
    #folder_pfam_name = "/Users/pauline/Desktop/data/Pfam_fasta_99"  
    #folder_pfam_name = "/Users/pauline/Desktop/data/Pfam_fasta_99_trimmed"  
    print(folder_pfam_name)
    dataCountDescription(folder_pfam_name)


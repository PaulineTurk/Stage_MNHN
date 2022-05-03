import mainNeighbor
import sys



path_folder_pId = sys.argv[1]
path_folder_data_split = sys.argv[2]
path_new_folder = sys.argv[3]



list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
#list_delay_number = [k for k in range(-10, 11) if k!=0]
list_delay_number = [-1, 1, -2, 2]
name_NeighborRes = "NeighborRes"
name_data_train = "Pfam_Train"
path_folder_fasta_train = f"{path_folder_data_split}/{name_data_train}" 
path_NeighborRes = f"{path_new_folder}/{name_NeighborRes}"
percentage_train = 50





for delay_num in list_delay_number:
    for kp_SeqChoice in ["k", "p"]:   
        path_proba_cond =  mainNeighbor.simpleContextualBlosum(path_folder_fasta_train, percentage_train, path_folder_pId, path_NeighborRes, delay_num, kp_SeqChoice, list_residu, pid_inf = 62)
        pass
print(path_folder_pId)
print("script_1")
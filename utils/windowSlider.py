import fastaReader
import numpy as np
import character as ch


def contextCatcher(seq_k, seq_p, len_seq,  position_k, len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p):
    # seq_k et seq_p après selection selon PID et diff name
    # len_seq = len(seq_k) = len(seq_p)
    # position_k: position de aa_k dans seq_k

    #print("position_k", position_k)
    #print("position_k - len_window_left_k", position_k - len_window_left_k)
    #print("position_k - len_window_left_k < 0", position_k - len_window_left_k < 0)
    #print("int(position_k - len_window_left_k < 0)", int((position_k - len_window_left_k) < 0))


    if int((position_k - len_window_left_k) >= 0) * int(position_k + len_window_right_k +1 <= len_seq) * int(position_k - len_window_left_p >= 0) * int(position_k + len_window_right_p +1 <= len_seq):
        # on suppose qu'on ne considère un aa et son voisinage que si tous ces voisinages sont définis (selon la taille des fenetres choisies en avance)
        # Rq. aa_k + context with context depending on the windows (tout le context sera écrit, mais pour le calcul de sa proba, ne tenir compte des résidus inclus)
        print(position_k)
        if len_window_left_k > 0:
            window_left_k = seq_k[position_k - len_window_left_k: position_k]
        else:
            window_left_k = []

        if len_window_right_k > 0:
            window_right_k = seq_k[position_k +1: position_k + len_window_right_k +1]
        else:
            window_right_k = []

        if len_window_left_p > 0:
            window_left_p = seq_p[position_k - len_window_left_p: position_k]
        else:
            window_left_k = []

        if len_window_right_p > 0:
            window_right_p = seq_p[position_k +1: position_k + len_window_right_p +1]
        else:
            window_right_p = []

        print("window_left_k", window_left_k)
        print("window_right_k", window_right_k)
        print("window_left_p", window_left_p)
        print("window_right_p", window_right_p)
    else:
        print("invalid position:", position_k)
    print("")







def windowSlider(accession_num, path_folder_pid, file_fasta, pid_inf, windowContext_visited, delay_num, kp_SeqChoice):    
    # k: known
    # p: predict
    # c: context
    # kp_SeqChoice: choice the reference sequence to look at its neighbor 
    liSeqAliFiltre = fastaReader.readFastaMul(file_fasta)
    pid_couple = np.load(path_folder_pid + "/" + accession_num+ ".pId.npy", allow_pickle='TRUE').item()
    list_AA = ch.characterList()


    if liSeqAliFiltre:
        for name_k, seq_k in liSeqAliFiltre:
            len_seq = len(seq_k)
            for name_p, seq_p in liSeqAliFiltre:
                if name_k != name_p:
                    if pid_couple[name_k][name_p] >= pid_inf:
                        for aa_index in range(len_seq):       
                            aa_k = seq_k[aa_index] 
                            aa_p = seq_p[aa_index] 
                            if all(x in list_AA for x in [aa_k, aa_p]):
                                index_neighbor = aa_index + delay_num
                                aa_c = seq_c[index_neighbor] 
                                if aa_c in list_AA: 
                                    triplet_count[aa_k][aa_p][aa_c] += 1
    else:
        print(accession_num)
    return triplet_count



if __name__ == '__main__': 
    seq_p = "ARNDCQEGHN"
    seq_k = "AR-LCQEGHI"
    len_seq = len(seq_k)
    len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p = 3, 2, 5, 1
    for index in range(len_seq):
        contextCatcher(seq_k, seq_p, len_seq, index, len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p)
    
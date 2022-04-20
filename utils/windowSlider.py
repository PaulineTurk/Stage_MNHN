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





if __name__ == '__main__': 
    seq_p = "ARNDCQEGHN"
    seq_k = "AR-LCQEGHI"
    len_seq = len(seq_k)
    len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p = 3, 2, 5, 1
    for index in range(len_seq):
        contextCatcher(seq_k, seq_p, len_seq, index, len_window_left_k, len_window_right_k, len_window_left_p, len_window_right_p)
    
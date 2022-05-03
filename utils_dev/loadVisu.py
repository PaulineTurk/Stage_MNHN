import numpy as np
import pandas as pd
import mainNeighbor

def loadAndVisualisation(path_dico):
    print(path_dico)
    dico_loaded = np.load(path_dico, allow_pickle='TRUE').item()
    df_dico_loaded = np.transpose(pd.DataFrame.from_dict(dico_loaded))  # not sure about transpose
    print(df_dico_loaded)
    return dico_loaded


import numpy as np
import pandas as pd
cond_proba = np.load("/Users/pauline/Desktop/data_Test_3/test_1/NeighborRes/proba_cond_(3,p) _percentage_train_90.npy", allow_pickle='TRUE').item()
df_cond_proba = np.transpose(pd.DataFrame.from_dict(cond_proba))
print(df_cond_proba)  

import numpy as np
import pandas as pd
from neuroHarmonize import harmonizationLearn
import pickle
import openpyxl

# load and split data by sex
data = pd.read_csv('data/deriv/HC_data.csv')
is_male = (data['sex'] == 'MALE').to_list()
is_female = np.logical_not(is_male)

data_M = data[is_male]
data_F = data[is_female]

## Prepare covariates
covars_M = data_M.iloc[:, 2:4] # select site and age
covars_M.columns = ['SITE', 'AGE']
covars_F = data_F.iloc[:, 2:4] # select site and age
covars_F.columns = ['SITE', 'AGE']

# Prepare data for harmonization Step 1
icv_M = np.array(data_M.iloc[:, 5])
icv_F = np.array(data_F.iloc[:, 5])

# Run harmonization Step 1
mod_icv_M, icv_M_adj = harmonizationLearn(np.stack((icv_M, icv_M)).T, covars_M, smooth_terms=['AGE'], eb=False) # stack ICV to fit dimensions required by harmonizationLearn
mod_icv_F, icv_F_adj = harmonizationLearn(np.stack((icv_F, icv_F)).T, covars_F, smooth_terms=['AGE'], eb=False) 

# prep data for step 2
covars_M['ICV_adj'] = icv_M_adj[:, 0]
covars_F['ICV_adj'] = icv_F_adj[:, 0]

rois_M = np.array(data_M.iloc[:, -145:])
rois_F = np.array(data_F.iloc[:, -145:])

# run  harmonization Step 2
mod_rois_M, rois_M_adj = harmonizationLearn(rois_M, covars_M, smooth_terms=['AGE'], eb=True) 
mod_rois_F, rois_F_adj = harmonizationLearn(rois_F, covars_F, smooth_terms=['AGE'], eb=True)

# join data into DataFrame
data_M_adj = pd.concat([data_M.iloc[:, :5],
                        pd.DataFrame(icv_M_adj[:, 0], 
                                     index=data_M.index, 
                                     columns=[data_M.columns[5]]
                                     ),
                        pd.DataFrame(rois_M_adj, 
                                     index=data_M.index,
                                     columns=data_M.columns[-145:]
                                     )
                        ], axis=1)

data_F_adj = pd.concat([data_F.iloc[:, :5],
                        pd.DataFrame(icv_F_adj[:, 0], 
                                     index=data_F.index, 
                                     columns=[data_F.columns[5]]
                                     ),
                        pd.DataFrame(rois_F_adj, 
                                     index=data_F.index,
                                     columns=data_F.columns[-145:]
                                     )
                        ], axis=1)
                        
data_adj = pd.concat([data_M_adj, data_F_adj], axis=0).sort_index()

# save adjusted DataFrame as .csv
data_adj.to_csv('data/deriv/HC_data_adj-ComBat-GAM.csv')

# save models for later use
with open('results/models/ComBat-GAM_icv_M', mode = 'wb') as file:
    pickle.dump(mod_icv_M, file)
with open('results/models/ComBat-GAM_icv_F', mode = 'wb') as file:
    pickle.dump(mod_icv_F, file)
with open('results/models/ComBat-GAM_rois_M', mode = 'wb') as file:
    pickle.dump(mod_rois_M, file)    
with open('results/models/ComBat-GAM_rois_F', mode = 'wb') as file:
    pickle.dump(mod_rois_F, file)

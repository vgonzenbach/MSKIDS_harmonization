import os
import pickle
import pandas as pd
import numpy as np
from neuroHarmonize import harmonizationApply

# Specify directories
DATA_DIR = os.path.abspath(os.path.join('', '..', 'data/deriv'))
MODELS_DIR = os.path.abspath(os.path.join('', '..', 'results/models'))

# Load MS data and split by sex
data = pd.read_csv(os.path.join(DATA_DIR, 'MS_data.csv'))
is_male = (data['sex'] == 'MALE').to_list()
is_female = np.logical_not(is_male)

data_M = data[is_male]
data_F = data[is_female]

## Prepare covariates

# Load PNC data to match number of sites
data_PNC = pd.read_csv(os.path.join(DATA_DIR, 'HC_data.csv'))
a_PNC_male = data_PNC[data_PNC['site'] == 'PNC'].iloc[0, :] # a Male
a_PNC_female = data_PNC[data_PNC['site'] == 'PNC'].iloc[1, :] # a Female

data_M = data_M.append(a_PNC_male)
data_F = data_F.append(a_PNC_female)

covars_M = data_M.iloc[:, 2:4] # select site and age
covars_M.columns = ['SITE', 'AGE']
covars_F = data_F.iloc[:, 2:4] # select site and age
covars_F.columns = ['SITE', 'AGE']

## Prepare data for harmonization Step 1
icv_M = np.array(data_M.iloc[:, 5])
icv_F = np.array(data_F.iloc[:, 5])

# Load ICV model for males and females
with open(os.path.join(MODELS_DIR, 'ComBat-GAM_icv_M'), 'rb') as f:
    model_icv_M = pickle.load(f)
with open(os.path.join(MODELS_DIR, 'ComBat-GAM_icv_F'), 'rb') as f:
    model_icv_F = pickle.load(f)

# Apply harmonization model Step 1
icv_M_adj = harmonizationApply(icv_M, covars_M, model_icv_M)
icv_F_adj = harmonizationApply(icv_F, covars_F, model_icv_F)

## Prepare data for harmonization Step 2
covars_M['ICV_adj'] = icv_M_adj[:, 0]
covars_F['ICV_adj'] = icv_F_adj[:, 0]

rois_M = np.array(data_M.iloc[:, -145:])
rois_F = np.array(data_F.iloc[:, -145:])

# Load ROI model for males and females
with open(os.path.join(MODELS_DIR, 'ComBat-GAM_rois_M'), 'rb') as f:
    model_rois_M = pickle.load(f)
with open(os.path.join(MODELS_DIR, 'ComBat-GAM_rois_F'), 'rb') as f:
    model_rois_F = pickle.load(f)

# run  harmonization Step 2
rois_M_adj = harmonizationApply(rois_M, covars_M, model_rois_M) 
rois_F_adj = harmonizationApply(rois_F, covars_F, model_rois_F)

# join data into DataFrame
data_M_adj = pd.concat([data_M.iloc[:, :5], # Demographics
                        pd.DataFrame(icv_M_adj[:, 0],  # ICV
                                     index=data_M.index, 
                                     columns=[data_M.columns[5]]
                                     ),
                        pd.DataFrame(rois_M_adj,  # ROIS
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
data_adj.to_csv(os.path.join(DATA_DIR, 'MS_data_adj-ComBat-GAM.csv'), index=False)
import os
import pickle
import pandas as pd
import numpy as np
from neuroHarmonize import harmonizationApply

# Specify directories
PROJECT_ROOT = "/Users/vgonzenb/PennSIVE/MSKIDS/"
DATA_DIR = os.path.join("data", "deriv")
MODELS_DIR = os.path.join("results", "models")

def apply_harmonization_by_sex(data_file = "MS_data.csv", mod = "GAM"):
     # Load MS data and split by sex
    data = pd.read_csv(os.path.join(PROJECT_ROOT + DATA_DIR, data_file))

    data_M = data[data.sex == "MALE"]
    data_F = data[data.sex == "FEMALE"]

    ## Prepare covariates

    # Load PNC data to match number of sites
    data_HC = pd.read_csv(os.path.join(PROJECT_ROOT + DATA_DIR, 'HC_data.csv'))

    a_PNC_male = data_HC[(data_HC.site == "PNC") & (data_HC.sex == "MALE")].iloc[0, :]
    a_PNC_female = data_HC[(data_HC.site == "PNC") & (data_HC.sex == "FEMALE")].iloc[0, :]

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
    if mod == "GAM":
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-GAM_ICV_from-HC_data_MALE_eb-True.pickle'), 'rb') as f:
            model_icv_M = pickle.load(f)
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-GAM_ICV_from-HC_data_FEMALE_eb-True.pickle'), 'rb') as f:
            model_icv_F = pickle.load(f)
    elif mod == "Linear":
        covars_M["AGE_SQUARED"] = (covars_M.AGE - covars_M.AGE.mean()) ** 2 # center before squaring (second moment)
        covars_F["AGE_SQUARED"] = (covars_F.AGE - covars_F.AGE.mean()) ** 2

        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-Linear_ICV_from-HC_data_MALE_eb-True.pickle'), 'rb') as f:
            model_icv_M = pickle.load(f)
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-Linear_ICV_from-HC_data_FEMALE_eb-True.pickle'), 'rb') as f:
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
    if mod == "GAM":
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-GAM_ROIs_from-HC_data_MALE_eb-True.pickle'), 'rb') as f:
            model_rois_M = pickle.load(f)
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-GAM_ROIs_from-HC_data_FEMALE_eb-True.pickle'), 'rb') as f:
            model_rois_F = pickle.load(f)
    if mod == "Linear":
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-Linear_ROIs_from-HC_data_MALE_eb-True.pickle'), 'rb') as f:
            model_rois_M = pickle.load(f)
        with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-Linear_ROIs_from-HC_data_FEMALE_eb-True.pickle'), 'rb') as f:
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

    # concatenate Data Frame for Males and Females
    data_adj = pd.concat([data_M_adj, data_F_adj], axis=0).sort_index()

    # drop PNC rows from data to save
    data_adj = data_adj.drop(index = [data_M_adj.index[-1], data_F_adj.index[-1]])
    
    # save adjusted DataFrame as .csv
    data_adj.to_csv(os.path.join(PROJECT_ROOT + DATA_DIR, f'{data_file[:-4]}_adj-ComBat-{mod}_from-HC_split-MF_eb-True.csv'), index=False)

apply_harmonization_by_sex(mod="Linear")

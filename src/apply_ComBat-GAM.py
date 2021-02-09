import os
import pickle
import pandas as pd
import numpy as np
from neuroHarmonize import harmonizationApply

# Specify directories
PROJECT_ROOT = "/Users/vgonzenb/PennSIVE/MSKIDS/"
DATA_DIR = os.path.join("data", "deriv")
MODELS_DIR = os.path.join("results", "models")

def apply_harmonization(data_file = "MS_data.csv"):
    """
    Applies harmonization from Healthy Controls to MS data. Will become a generalized function for future projects.
    """
     # Load MS data and split by sex
    data = pd.read_csv(os.path.join(PROJECT_ROOT + DATA_DIR, data_file))

    ## Prepare covariates

    # Load PNC data to match number of sites
    data_HC = pd.read_csv(os.path.join(PROJECT_ROOT + DATA_DIR, 'HC_data.csv'))
    
    a_PNC_male = data_HC[(data_HC.site == "PNC") & (data_HC.sex == "MALE")].iloc[0, :]
    a_PNC_female = data_HC[(data_HC.site == "PNC") & (data_HC.sex == "FEMALE")].iloc[0, :]

    data = data.append(a_PNC_male)
    data = data.append(a_PNC_female)

    covars = data.iloc[:, 2:4] # select site, age and sex
    covars.columns = ['SITE', 'AGE']
    covars = covars.join(pd.get_dummies(data.sex))
    
    ## Prepare data for harmonization Step 1
    icv = np.array(data.iloc[:, 5])

    # Load ICV model for males and females
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-GAM_ICV_from-HC_data_eb-True.pickle'), 'rb') as f:
        model_icv = pickle.load(f)

    # Apply harmonization model Step 1
    icv_adj = harmonizationApply(icv, covars, model_icv)

    ## Prepare data for harmonization Step 2
    covars['ICV_adj'] = icv_adj[:, 0]

    ROIs = np.array(data.iloc[:, -145:])

    # Load ROI model for males and females
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, 'ComBat-GAM_ROIs_from-HC_data_eb-True.pickle'), 'rb') as f:
        model_ROIs = pickle.load(f)

    # run  harmonization Step 2
    ROIs_adj = harmonizationApply(ROIs, covars, model_ROIs) 

    # join data into DataFrame
    data_adj = pd.concat([data.iloc[:, :5], # Demographics
                          pd.DataFrame(icv_adj[:, 0],  # ICV
                                       index=data.index, 
                                       columns=[data.columns[5]]
                                       ),
                          pd.DataFrame(ROIs_adj,  # ROIS
                                       index=data.index,
                                       columns=data.columns[-145:]
                                       )
                            ], axis=1)

    # drop PNC rows from data to save
    data_adj = data_adj.drop(index = data_adj.index[-2:])
    
    # save adjusted DataFrame as .csv
    data_adj.to_csv(os.path.join(PROJECT_ROOT + DATA_DIR, 'MS_data_adj-ComBat-GAM_from-HC_split-None_eb-True.csv'), index=False)

apply_harmonization()
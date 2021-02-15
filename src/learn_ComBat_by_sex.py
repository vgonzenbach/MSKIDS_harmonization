import os
import numpy as np
import pandas as pd
from neuroHarmonize import harmonizationLearn
import pickle

PROJECT_ROOT = "/Users/vgonzenb/PennSIVE/MSKIDS/"
MODELS_DIR = os.path.join("results", "models")
DATA_DIR = os.path.join("data", "deriv")


def harmonize_data_by_sex(data_file = 'HC_data.csv', mod="GAM", eb=True):
    # Load and split data by sex
    data = pd.read_csv(os.path.join(PROJECT_ROOT + DATA_DIR, data_file))
    data_M = data[data.sex == "MALE"]
    data_F = data[data.sex == "FEMALE"]

    # Prepare covariates
    covars_M = data_M.iloc[:, 2:4] # select site and age
    covars_M.columns = ['SITE', 'AGE']
    covars_F = data_F.iloc[:, 2:4] # select site and age
    covars_F.columns = ['SITE', 'AGE']

    # Prepare data for harmonization Step 1
    icv_M = np.array(data_M.iloc[:, 5])
    icv_F = np.array(data_F.iloc[:, 5])

    # Run harmonization Step 1
    if mod == "GAM":
        mod_icv_M, icv_M_adj = harmonizationLearn(np.stack((icv_M, icv_M)).T, covars_M, smooth_terms=['AGE'], eb=False) # stack ICV to fit dimensions required by harmonizationLearn
        mod_icv_F, icv_F_adj = harmonizationLearn(np.stack((icv_F, icv_F)).T, covars_F, smooth_terms=['AGE'], eb=False) 
    elif mod == "Linear":
        covars_M["AGE_SQUARED"] = (covars_M.AGE - covars_M.AGE.mean()) ** 2 # center before squaring (second moment)
        covars_F["AGE_SQUARED"] = (covars_F.AGE - covars_F.AGE.mean()) ** 2

        mod_icv_M, icv_M_adj = harmonizationLearn(np.stack((icv_M, icv_M)).T, covars_M, eb=False) # stack ICV to fit dimensions required by harmonizationLearn
        mod_icv_F, icv_F_adj = harmonizationLearn(np.stack((icv_F, icv_F)).T, covars_F, eb=False) 
    else:
        print("Please specify a type of model")
        return None

    # prep data for step 2
    covars_M['ICV_adj'] = icv_M_adj[:, 0]
    covars_F['ICV_adj'] = icv_F_adj[:, 0]

    rois_M = np.array(data_M.iloc[:, -145:])
    rois_F = np.array(data_F.iloc[:, -145:])

    # run  harmonization Step 2
    if mod == "GAM":
        mod_rois_M, rois_M_adj = harmonizationLearn(rois_M, covars_M, smooth_terms=['AGE'], eb=eb) 
        mod_rois_F, rois_F_adj = harmonizationLearn(rois_F, covars_F, smooth_terms=['AGE'], eb=eb)
    elif mod == "Linear":
        mod_rois_M, rois_M_adj = harmonizationLearn(rois_M, covars_M, eb=eb) 
        mod_rois_F, rois_F_adj = harmonizationLearn(rois_F, covars_F, eb=eb)

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
    data_adj.to_csv(os.path.join(PROJECT_ROOT + DATA_DIR, f'{data_file[:-4]}_adj_ComBat-{mod}_split-MF_eb-{eb}.csv'), index=False)

    # save models for later use
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, f'ComBat-{mod}_ICV_from-{data_file[:-4]}_MALE_eb-{eb}.pickle'), mode = 'wb') as file:
        pickle.dump(mod_icv_M, file)
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, f'ComBat-{mod}_ICV_from-{data_file[:-4]}_FEMALE_eb-{eb}.pickle'), mode = 'wb') as file:
        pickle.dump(mod_icv_F, file)
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, f'ComBat-{mod}_ROIs_from-{data_file[:-4]}_MALE_eb-{eb}.pickle'), mode = 'wb') as file:
        pickle.dump(mod_rois_M, file)    
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, f'ComBat-{mod}_ROIs_from-{data_file[:-4]}_FEMALE_eb-{eb}.pickle'), mode = 'wb') as file:
        pickle.dump(mod_rois_F, file)

harmonize_data_by_sex(mod="Linear")
import os
import numpy as np
import pandas as pd
from neuroHarmonize import harmonizationLearn
import pickle

PROJECT_ROOT = "/Users/vgonzenb/PennSIVE/MSKIDS/"
MODELS_DIR = os.path.join("results", "models")
DATA_DIR = os.path.join("data", "deriv")


def harmonize_data(data_file = 'HC_data.csv', mod="GAM", eb=True):

    data = pd.read_csv(os.path.join(PROJECT_ROOT + DATA_DIR, data_file))

    # Prepare covariates
    covars = data.iloc[:, 2:4] # select site, age and sex
    covars.columns = ['SITE', 'AGE']
    covars = covars.join(pd.get_dummies(data.sex))

    # Prepare data for harmonization Step 1
    icv = np.array(data.iloc[:, 5])

    # Run harmonization Step 1
    if mod == "GAM":
        mod_ICV, ICV_adj = harmonizationLearn(np.stack((icv, icv)).T, covars, smooth_terms=['AGE'], eb=False) # stack ICV to fit dimensions required by harmonizationLearn
    elif mod == "Linear":
        mod_ICV, ICV_adj = harmonizationLearn(np.stack((icv, icv)).T, covars, eb=False) # stack ICV to fit dimensions required by harmonizationLearn
    else:
        print("Please specify a type of model")
        return None
    # prep data for step 2
    covars['ICV_adj'] = ICV_adj[:, 0]

    ROIs = np.array(data.iloc[:, -145:])

    # run  harmonization Step 2
    if mod == "GAM":
        mod_ROIs, ROIs_adj = harmonizationLearn(ROIs, covars, smooth_terms=['AGE'], eb=eb) 
    elif mod == "Linear":
        mod_ROIs, ROIs_adj = harmonizationLearn(ROIs, covars, eb=eb) 
    
    # join data into DataFrame
    data_adj = pd.concat([data.iloc[:, :5],
                          pd.DataFrame(ICV_adj[:, 0], 
                                       index=data.index, 
                                       columns=[data.columns[5]]
                                       ),

                          pd.DataFrame(ROIs_adj, 
                                       index=data.index,
                                       columns=data.columns[-145:]
                                       )
                          ], axis=1)

    # save adjusted DataFrame as .csv
    data_adj.to_csv(os.path.join(PROJECT_ROOT + DATA_DIR, f'{data_file[:-4]}_adj_ComBat-{mod}_split-None_eb-{eb}.csv'), index=False)

    # save models for later use
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, f'ComBat-{mod}_ICV_from-{data_file[:-4]}_eb-{eb}.pickle'), mode = 'wb') as file:
        pickle.dump(mod_ICV, file)
    with open(os.path.join(PROJECT_ROOT + MODELS_DIR, f'ComBat-{mod}_ROIs_from-{data_file[:-4]}_eb-{eb}.pickle'), mode = 'wb') as file:
        pickle.dump(mod_ROIs, file)    

harmonize_data()
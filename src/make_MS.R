library(readxl)
library(dplyr)
source('src/match_obs.R')

## =============================================CLEAN MSKIDS DATA========================================================
MSK_demogr_df = readxl::read_xlsx("data/base/Copy of MSKIDS MRI Inventory Penn_10272020.xlsx") # MSKIDS demographics
MSK_demogr_df = MSK_demogr_df[-(nrow(MSK_demogr_df):(nrow(MSK_demogr_df)-3)),] # exclude 'legend' rows

MSK_roivol_df = read.csv("data/base/MSKIDS_3.5D_2020_DLICV_N4N4_MUSE_Features_DerivedVolumes.csv") # MSKIDS volumetric data

# Subset demographics healthy controls 
MS = MSK_demogr_df %>% 
  filter(diagnosis_category != "HEALTHY CONTROL" & `MNI QC Result` == "ACCEPTED") %>% 
  select(`Study ID`, Scanner, `Site of Acquisition`, `Age at Scan`, sex, diagnosis_category, Visit)

MS %>% colnames

# Exclude unused sites
MS = MS[MS$`Site of Acquisition` %in% c("CHP", "HSC"),]

# Match demographics table with volumes
res = match_obs(ref_id = paste0(MS$`Study ID`, "-", tolower(MS$Visit)), target_id = MSK_roivol_df$ID)

MS_data = cbind(MS[!is.na(res$row_matches),], # select rows where there was a match and
                MSK_roivol_df[res$row_matches[!is.na(res$row_matches)], -1]) # bind to specific row

# Select columns for ComBat adjustment
MS_data = MS_data[, c(1:5, 8, 122:266)]

colnames(MS_data)[1:5] = c("ID", "scanner", "site", "age", "sex")

# Check consistency between site and scanner
table(MS_data$scanner, MS_data$site)

# Delete n=18 HSC-SIEMENSTIMTRIO
MS_data = MS_data[!(MS_data$site == "HSC" & MS_data$scanner == "SIEMENSTIMTRIO"),]

write.csv(MS_data, 'data/deriv/MS_data.csv', row.names = FALSE)

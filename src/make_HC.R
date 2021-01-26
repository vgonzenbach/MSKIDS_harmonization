#' 
#' Applying ComBat to ROI volumes in MSKIDS data
#' Virgilio Gonzenbach  
#' PennSIVE
#' 
library(readxl)
library(dplyr)
source('src/match_obs.R')

## =============================================CLEAN MSKIDS DATA========================================================
MSK_demogr_df = readxl::read_xlsx("data/base/Copy of MSKIDS MRI Inventory Penn_10272020.xlsx") # MSKIDS demographics
MSK_demogr_df = MSK_demogr_df[-(nrow(MSK_demogr_df):(nrow(MSK_demogr_df)-3)),] # exclude 'legend' rows

MSK_roivol_df = read.csv("data/base/MSKIDS_3.5D_2020_DLICV_N4N4_MUSE_Features_DerivedVolumes.csv") # MSKIDS volumetric data

# Subset demographics healthy controls 
MSK_HC = MSK_demogr_df %>% 
  filter(diagnosis_category == "HEALTHY CONTROL" & `MNI QC Result` == "ACCEPTED") %>% 
  select(`Study ID`, Scanner, `Site of Acquisition`, `Age at Scan`, sex)

MSK_MS = MSK_demogr_df %>% 
  filter(diagnosis_category != "HEALTHY CONTROL" & `MNI QC Result` == "ACCEPTED") %>% 
  select(`Study ID`, Scanner, `Site of Acquisition`, `Age at Scan`, sex, diagnosis_category, Visit)

saveRDS(MSK_MS, "MS_data.rds")

# Map IDs from ref data.frame to target data.frame----

res = match_obs(ref_id = MSK_HC$`Study ID`, target_id = MSK_roivol_df$ID)
MSK_HC = cbind(MSK_HC[!is.na(res$row_matches),], MSK_roivol_df[res$row_matches[!is.na(res$row_matches)], -1])
   
## =============================================CLEAN PNC DATA========================================================

PNC_demogr_df = read.csv("data/base/PNC/demographics_from_20160207_dataRelease_update20161114.csv", stringsAsFactors = FALSE)
PNC_roivol_df = read.csv("data/base/PNC/GO-BBL_muse_dramms+ants_C1.2_Features_DerivedVolumes.csv", stringsAsFactors = FALSE) # PNC volumetric data

# Load filters
ltnExclude = read.csv("data/base/PNC/healthexclude_ltn.csv")[, "ltnExclude"]
t1Exclude = read.csv("data/base/PNC/n1601_t1QaData_v2.csv")[, "t1Exclude"]

# Subject-level exclusion and select columns
PNC = PNC_demogr_df %>% filter(!ltnExclude & !t1Exclude) %>% select(bblid, scanid, ageAtGo1Scan, sex)

ref_id = paste(PNC$bblid, PNC$scanid, sep = "-")
res = match_obs(ref_id = ref_id, target_id = PNC_roivol_df$ID)

PNC = cbind(PNC, PNC_roivol_df[res$row_matches, -1])

# Handle duplicate for 114738-4738
res$duplicates

# Manually select first match as correct, and map this row
index = match(res$duplicates$dup_targets[[1]][1], PNC_roivol_df$ID)


PNC[which(res$row_matches %>% is.na()), 5:ncol(PNC)] = PNC_roivol_df[index, -1] # Note: Use new.env for cleaner code


# create scanner and site column for PNC data 
# have only one id for PNC data
colnames(MSK_HC)[1:10]
colnames(PNC)[1:10]


row_id = data.frame(ID = paste0("PNC-", PNC$bblid), 
                    scanner = "SIEMENSTIMTRIO", 
                    site = "PNC")

PNC = cbind(row_id, PNC[,-(1:2)])

# change age column name
colnames(PNC)[4] = "age"
colnames(MSK_HC)[1:4] = c("ID","scanner","site", "age")


data = rbind(MSK_HC, PNC)
data$sex = recode(data$sex, `1` = "MALE", `2` = "FEMALE")

# Exclude n=7 batch
data = data[!(data$site == "HSC" & data$scanner == "SIEMENTIMSTRIO"),]

# Exclude extra ROIs
data = data[,c(1:6, 120:ncol(data))]

write.csv(data, "data/deriv/HC_data.csv", row.names = FALSE)








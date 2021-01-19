## Harmonize ICV and ROIs by sex
library(neuroCombat)
library(dplyr)

df = readRDS('data/deriv/HC_data.rds')
M_filter = df$sex == "MALE"
F_filter = df$sex == "FEMALE"

colnames(df)[6] = "ICV"
df$site = as.factor(df$site)

df_adj = df
df_adj[,6:ncol(df_adj)] = matrix(nrow = 1286, ncol = 146)

# Combat Harmonization Step 1
mod_M = model.matrix(~ age, data = df[M_filter,])
ICV_M = neuroCombat::neuroCombat(dat = t(df$ICV[M_filter]),
                                 batch = df$site[M_filter], 
                                 mod = mod_M,
                                 eb = FALSE,
                                 parametric = TRUE,
                                 ref.batch = "PNC")

mod_F = model.matrix(~ age, data = df[F_filter,])
ICV_F = neuroCombat::neuroCombat(dat = t(df$ICV[F_filter]), 
                                 batch = df$site[F_filter], 
                                 mod = mod_F,
                                 eb = FALSE,
                                 parametric=TRUE,
                                 ref.batch = "PNC")

df_adj$ICV[M_filter] = t(ICV_M$dat.combat)
df_adj$ICV[F_filter] = t(ICV_F$dat.combat)


## ComBat Harmonization Step 2
mod_M = model.matrix(~ age + ICV, data = df[M_filter,])
resCombat_M = neuroCombat::neuroCombat(dat = t(df[M_filter, 7:ncol(df)]), 
                                       batch = df$site[M_filter], 
                                       mod = mod_M,
                                       eb = TRUE,
                                       parametric= TRUE,
                                       ref.batch = "PNC")

mod_F = model.matrix(~ age + ICV, data = df[F_filter,])
resCombat_F = neuroCombat::neuroCombat(dat = t(df[F_filter, 7:ncol(df)]), 
                                       batch = df$site[F_filter], 
                                       mod = mod_F,
                                       eb = TRUE,
                                       parametric= TRUE,
                                       ref.batch = "PNC")

df_adj[M_filter, 7:ncol(df_adj)] = t(resCombat_M$dat.combat)
df_adj[F_filter, 7:ncol(df_adj)] = t(resCombat_F$dat.combat)

# Save adjusted data
write.csv(df_adj, 'data/deriv/HC_data_adj-ComBat-Linear.csv')

# Save models
saveRDS(ICV_M, 'results/models/ComBat-Linear_ICV_M.rds')
saveRDS(ICV_F, 'results/models/ComBat-Linear_ICV_F.rds')
saveRDS(resCombat_M, 'results/models/ComBat-Linear_ROIs_M.rds')
saveRDS(resCombat_F, 'results/models/ComBat-Linear_ROIs_F.rds')



library(factoextra)
library(dplyr)
library(FactoMineR)

df = readRDS("data/deriv/HC_data.rds")
dict = readxl::read_xlsx("MUSE_ROI_Dict.xlsx")

df = df[df$site == "PNC",]

# PCA
indexes = paste0("X", 1:300)
indexes = indexes[indexes %in% colnames(df)] # keep indexes that match a column name

vol = data.frame(df[,indexes])


resPCA = FactoMineR::PCA(vol, scale.unit = FALSE, graph = FALSE)

## 1 and 2 ----
fviz_pca_var(resPCA)

fviz_pca_ind(resPCA, geom = "point", col.ind = df$sex, alpha.ind = 0.5, palette = "npg")

˜## 3 and 4 ----
fviz_pca_var(resPCA, axes = c(3,4))
fviz_pca_ind(resPCA, axes = c(3,4), geom = "point", col.ind = df$sex, alpha.ind = 0.5, palette = "npg")

# Sex differences in components
dat = data.frame(pc1 = resPCA$ind$coord[,1], 
                 pc2 = resPCA$ind$coord[,2], 
                 pc3 = resPCA$ind$coord[,3],
                 pc4 = resPCA$ind$coord[,4],
                 sex = df$sex, 
                 age = df$age)


lm(scale(pc1) ~ sex, data = dat) %>% summary()
lm(scale(pc2) ~ sex, data = dat) %>% summary()
lm(scale(pc3) ~ sex, data = dat) %>% summary()
lm(scale(pc3) ~ sex, data = dat) %>% summary()

# Account for age
lm(scale(pc1) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc2) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc3) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc4) ~ sex + scale(age), data = dat) %>% summary()

# Test interactions
lm(scale(pc1) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc2) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc3) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc4) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()

# Normalize with ICV -----
norm_vol = vol / df$X702

resPCA = FactoMineR::PCA(norm_vol, 
                         scale.unit = FALSE, 
                         graph = FALSE, 
                         ncp = 8)

## 1 and 2 ----
plot.PCA(resPCA, choix = "var", cex = 1, autoLab = "yes")
fviz_pca_ind(resPCA, geom = "point", 
             col.ind = df$sex, 
             alpha.ind = 0.5, 
             palette = "npg")

## 3 and 4 ----
plot.PCA(resPCA, axes = c(3,4), choix = "var", cex = 1, autoLab = "yes")
fviz_pca_ind(resPCA, axes = c(3,4), 
             geom = "point", 
             col.ind = df$sex, 
             alpha.ind = 0.5, 
             palette = "npg")

## 5 and 6 ----
plot.PCA(resPCA, axes = c(5,6), choix = "var", cex = 1, autoLab = "yes")
fviz_pca_ind(resPCA, axes = c(5,6), 
             geom = "point", 
             col.ind = df$sex, 
             alpha.ind = 0.5, 
             palette = "npg")

## 7 and 8
plot.PCA(resPCA, axes = c(7,8), choix = "var", cex = 1, autoLab = "yes")
fviz_pca_ind(resPCA, axes = c(7,8), 
             geom = "point", 
             col.ind = df$sex, 
             alpha.ind = 0.5, 
             palette = "npg")

# Sex differences in components
dat = data.frame(pc1 = resPCA$ind$coord[,1], 
                 pc2 = resPCA$ind$coord[,2], 
                 pc3 = resPCA$ind$coord[,3],
                 pc4 = resPCA$ind$coord[,4],
                 pc5 = resPCA$ind$coord[,5],
                 pc6 = resPCA$ind$coord[,6],
                 pc7 = resPCA$ind$coord[,7],
                 pc8 = resPCA$ind$coord[,8],
                 sex = df$sex, 
                 age = df$age)


lm(scale(pc1) ~ sex, data = dat) %>% summary()
lm(scale(pc2) ~ sex, data = dat) %>% summary()
lm(scale(pc3) ~ sex, data = dat) %>% summary()
lm(scale(pc4) ~ sex, data = dat) %>% summary()
lm(scale(pc5) ~ sex, data = dat) %>% summary()
lm(scale(pc6) ~ sex, data = dat) %>% summary()
lm(scale(pc7) ~ sex, data = dat) %>% summary()
lm(scale(pc8) ~ sex, data = dat) %>% summary()

# Account for age
lm(scale(pc1) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc2) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc3) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc4) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc5) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc6) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc7) ~ sex + scale(age), data = dat) %>% summary()
lm(scale(pc8) ~ sex + scale(age), data = dat) %>% summary()

# Test interactions
lm(scale(pc1) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc2) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc3) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc4) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc5) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc6) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc7) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()
lm(scale(pc8) ~ sex + scale(age) + sex*scale(age), data = dat) %>% summary()



# Add bootstrapped means
# Plot PC values on brain
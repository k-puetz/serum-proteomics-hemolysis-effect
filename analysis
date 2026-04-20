#loading libraries
library(readxl)
library(dplyr)
library(scales)

#loading serum proteomics dataset
serum_proteomics <- read_excel('data.xlsx', sheet = 'preprocessed_data')

#removing HeLa QC samples
serum_proteomics <- serum_proteomics %>% 
  select(-contains("HeLa"))

#loading list of proteins known to be connected to hemolysis
hemolysis <- read_excel('hemolysis_proteins.xlsx')

#find hemolysis proteins in serum_proteomics
hemo_prots <- serum_proteomics[
  sapply(serum_proteomics$First.Protein.Description, function(x) {
    any(grepl(paste(hemolysis$hemolysis, collapse = "|"), x, ignore.case = TRUE))
  }),
]

#build a matrix of the hemolysis proteins
mat <- as.matrix(hemo_prots[, 7:ncol(hemo_prots)]) 

#pca of only the hemolysis proteins
pca_hemo <- prcomp(t(mat), scale. = TRUE)
plot(pca_hemo$x[,1], pca_hemo$x[,2])

#generating a hemo_score by using the average PC1 of the hemolysis proteins 
hemo_score <- pca_hemo$x[,1]

#generating a matrix of all serum proteomics
all_mat <- as.matrix(serum_proteomics[, 7:ncol(serum_proteomics)])
pca_all <- prcomp(t(all_mat), scale. = TRUE)

#coloring the total serum proteomics PCA by hemo_score (gradient from blue -> red)
cols <- col_numeric(
  palette = c("blue", "red"),
  domain = hemo_score
)

plot(pca_all$x[,1], pca_all$x[,2],
     col = cols(hemo_score),
     pch = 20)

#Correlation: is there a connection between the hemo_score an total variance?
cor(hemo_score, pca_all$x[,1])
cor(hemo_score, pca_all$x[,2])

#Regression: Does hemo_score explain the variance in the serum proteomic measurements?
summary(lm(pca_all$x[,1] ~ hemo_score))

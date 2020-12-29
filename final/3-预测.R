# 预测
library(caret)
library(randomForest)
library(tidyverse)
library(RANN)
library(DMwR)
library(pROC)
library(mixOmics)
## install BiocManager if not installed if (!requireNamespace("BiocManager", quietly = TRUE))     install.packages("BiocManager") 
## install mixOmics 
BiocManager::install('mixOmics')


# 导入标准化后的 T1 & T3 数据
pred.norm <- read.csv(file = "final/Metab_data_normalized_t_T1&T3.csv",stringsAsFactors = F)
pred.norm <- pred.norm %>%
  select(-sample) %>%
  mutate(group=factor(group,labels = c("T1","T3")))

# Tunning
set.seed(998)
inTraining <- createDataPartition(pred.norm$group, p = .75, list = FALSE)
training <- pred.norm[ inTraining,]
testing  <- pred.norm[-inTraining,]


# 1. upregulated protein

# 2. downregulated protein

# 3. all significant protein

# 4. top predicable protein


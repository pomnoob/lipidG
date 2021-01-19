# Try GBM for proteomic data
library(rsample)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # a java-based platform
library(pdp)          # model visualization
library(ggplot2)      # model visualization
library(lime)         # model visualization

# Model Tunning

pomic <- read.csv(file = "final/proteome data with z score.csv",
                  stringsAsFactors = F)

# Exclude data from Beijing 
pomic.exBJ <- pomic %>%
  filter(city!="北京")

pomic.exBJ.GBM <- pomic.exBJ %>%
  select(haz,2:473) %>%
  filter(!is.na(haz))

# Include data from Beijing
pomic.GBM <- pomic %>%
  select(waz,2:473) %>%
  filter(!is.na(waz))

# for reproducibility
set.seed(123)

# train GBM model
gbm.fit <- gbm(
  formula = haz ~ .,
  distribution = "gaussian",
  data = pomic.exBJ.GBM,
  n.trees = 10000,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.05,
  cv.folds = 10
)  

best <- which.min(gbm.fit$cv.error)
sqrt(gbm.fit$cv.error[best])
gbm.perf(gbm.fit, method = "cv")

# for reproducibility
set.seed(123)

# train GBM model
gbm.fitBJ <- gbm(
  formula = waz ~ .,
  distribution = "gaussian",
  data = pomic.GBM,
  n.trees = 10000,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.05,
  cv.folds = 10
)  

best <- which.min(gbm.fitBJ$cv.error)
sqrt(gbm.fitBJ$cv.error[best])
gbm.perf(gbm.fitBJ, method = "cv")


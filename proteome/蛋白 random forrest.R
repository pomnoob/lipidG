library(caret)
library(randomForest)
library(tidyverse)
library(RANN)
library(DMwR)
library(pROC)
################################################################################

# 先做 waz
load(file = "data/proteomic-infant weight for age.Rdata")
wazRF <- waz %>%
  select(3:139,twaz)%>%
  mutate(twaz=factor(twaz,labels = c("T1","T3")))

#missing value imputation
set.seed(1212)
wazRF.i <- rfImpute(twaz~.,data=wazRF,iter = 5)

##train data
set.seed(3456)
trainIndex <- createDataPartition(wazRF.i$twaz, p = .67, 
                                  list = FALSE, 
                                  times = 1)
training <- wazRF.i[ trainIndex,]
testing  <- wazRF.i[-trainIndex,]

# cross-validation method
fitControl <- trainControl(
  #Monte-Carlo CV
  method = "LGOCV",
  p=2/3,
  savePredictions = TRUE,
  classProbs = TRUE)

# Random Forrest
set.seed(825)
rfFit1 <- train(twaz ~ ., data = training, 
                method = "rf", 
                trControl = fitControl,
                summaryFunction = twoClassSummary)
rfFit1

rf.pred <- predict(rfFit1, testing, type="prob")


result.roc <- roc(testing$twaz, rf.pred$T3)
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")

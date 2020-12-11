# 使用差异显著的蛋白质做一下RF
library(caret)
library(randomForest)
library(tidyverse)
library(RANN)
library(DMwR)
library(pROC)

# 1. 使用未做校正差异显著的蛋白值做RF
## 筛选数据
### 选择差异显著的蛋白质
load(file = "data/不做校正差异显著的蛋白质 WAZ T1和T3两组之间t检验.Rdata")
load(file = "data/proteomic-infant weight for age.Rdata")
pid.sig <- protWAZ.p$pid
waz.sig <- waz %>%
  select(twaz,all_of(pid.sig))%>%
  mutate(twaz=factor(twaz,labels = c("T1","T3")))

#missing value imputation
set.seed(1212)
wazRF.i <- rfImpute(twaz~.,data=waz.sig,iter = 5)
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

# 2. 使用校正后差异仍然显著的蛋白做RF
## 筛选数据
### 选择差异显著的蛋白质
load(file = "data/FDR校正后差异显著的蛋白质 WAZ T1和T3两组之间t检验.Rdata")
load(file = "data/proteomic-infant weight for age.Rdata")
pid.adj <- protWAZ.adj$pid
waz.adj <- waz %>%
  select(twaz,all_of(pid.adj))%>%
  mutate(twaz=factor(twaz,labels = c("T1","T3")))

#missing value imputation
set.seed(1212)
wazRF.i <- rfImpute(twaz~.,data=waz.adj,iter = 5)
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
# 3. 对比1条件下RF筛选出的蛋白质与FDR校正后差异显著的蛋白质是否一致

# 4. 预测模型
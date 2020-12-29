# 预测
library(caret)
library(randomForest)
library(tidyverse)
library(RANN)
library(DMwR)
library(pROC)

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
# CV
fitControl <- trainControl(## 10-fold CV
  method = "LGOCV",
  number = 20)
#training
set.seed(825)
plsFit1 <- train(group ~ ., data = pred.norm, 
                 method = "pls", 
                 trControl = fitControl,
                 summaryFunction = twoClassSummary)
plsFit1

pls.pred <- predict(plsFit1, testing, type="prob")
result.roc <- roc(testing$group, pls.pred$T1)
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
# 1. upregulated protein

# 2. downregulated protein

# 3. all significant protein

# 4. top predicable protein


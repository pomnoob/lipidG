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
proNorm <- read.csv(file = "Download/data_normalized_id.csv",stringsAsFactors = F)
# 选择火山图筛选出的pid
tPID <- read.csv(file="Download/volcano.csv",stringsAsFactors = F)
tPID2 <- tPID$pid

t.test.pro <- proNorm %>%
  select(id,Label,all_of(tPID2))%>%
  mutate(Label=factor(Label,labels = c("Normal","Low")))

t.test.pro <- column_to_rownames(t.test.pro,var = "id")

# 选择VIP大于2的pid
plsPID <- read.csv(file="Download/plsda_vip.csv",stringsAsFactors = F)
plsPID <- plsPID %>%
  filter(Comp..1 > 2)
plsPID2 <- plsPID$pid

pls.pro <- proNorm %>%
  select(Label,all_of(plsPID2))%>%
  mutate(Label=factor(Label,labels = c("Normal","Low")))

## 火山图数据做预测
## 数据分成训练集和预测集
# cross-validation method
fitControl <- trainControl(
  #Monte-Carlo CV
  method = "LGOCV",
  p=2/3,
  savePredictions = TRUE,
  classProbs = TRUE)

# 火山图id预测
pls.pred.tall <- data.frame(id=character(),Label=integer(),Low=double())
for (i in 1:500) {
  #put this within the loop to get different training and testing data each time
  trainIndex.t <- createDataPartition(t.test.pro$Label, p = .67, 
                                    list = FALSE, 
                                    times = 1)
  training.t <- t.test.pro[ trainIndex.t,]
  testing.t  <- t.test.pro[-trainIndex.t,]
  ## case control 均衡
  smote_train.t <- SMOTE(Label ~ ., data  = training.t)                         
                           
  plsFit.t <- train(Label ~ ., data = training.t, 
                   method = "pls", 
                   trControl = fitControl,
                   summaryFunction = twoClassSummary)
  
  #to get prediction for ROC
  #predict the membership using testing data
  pls.pred.t <- predict(plsFit.t, testing.t, type="prob")
  pls.test.t <- testing.t %>%
    select(Label)%>%
    rownames_to_column(var="id")
  #get probability
  pls.pred.t <- pls.pred.t %>%
    select(Low) %>%
    rownames_to_column(var="id")
  pls.pt <- left_join(pls.test.t,pls.pred.t,by="id")
  pls.pred.tall <- rbind(pls.pred.tall,pls.pt)
}

#plot ROC
library(pROC)
par(pty="s")
result.roc <- roc(pls.pred.tall$Label, 
                  pls.pred.tall$Low,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, 
                  print.auc=TRUE,
                  legacy.axes=TRUE,
                  percent = TRUE,
                  xlab="False Positive Percentage",
                  ylab="True Positive Percentage",
                  col="black")

par(pty="m")

roc.t <- roc(pls.pred.tall$Label, 
                  pls.pred.tall$Low)

groc.t <- ggroc(roc.t,
              legacy.axes = TRUE,size=1)+
  theme_minimal()+
  
  theme(
    legend.position = c(0.5,0.5),
    panel.grid = element_blank(),
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(color = "black",size=16),
    axis.line = element_line(color="black"),
    axis.text = element_text(color="black",size=14),
    legend.text = element_text(color="black",size=16),
    legend.title = element_text(color="black",size=20)
  )+
  xlab("False Positive")+ylab("True Positive")+
  
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="grey", linetype="dashed",size=1.5)+
  geom_text(aes(x=0.6,y=0.4,label="AUC=0.875 (CI: 0.868-0.883)"),
            size=5,color="black")
groc.t

ggsave(filename = "figures/ROC of t test.svg",
       plot = groc.t,
       width=12,
       height = 10,
       dpi = 75)

# PLS id 预测
pls.pred.plsall <- data.frame(id=character(),Label=integer(),Low=double())
for (i in 1:500) {
  #put this within the loop to get different training and testing data each time
  trainIndex.pls <- createDataPartition(pls.pro$Label, p = .67, 
                                    list = FALSE, 
                                    times = 1)
  training.pls <- t.test.pro[ trainIndex.pls,]
  testing.pls  <- t.test.pro[-trainIndex.pls,]
  ## case control 均衡
  smote_train.pls <- SMOTE(Label ~ ., data  = training.pls)                         
  
  plsFit.pls <- train(Label ~ ., data = training.pls, 
                    method = "pls", 
                    trControl = fitControl,
                    summaryFunction = twoClassSummary)
  
  #to get prediction for ROC
  #predict the membership using testing data
  pls.pred.pls <- predict(plsFit.pls, testing.pls, type="prob")
  pls.test.pls <- testing.pls %>%
    select(Label)%>%
    rownames_to_column(var="id")
  #get probability
  pls.pred.pls <- pls.pred.pls %>%
    select(Low) %>%
    rownames_to_column(var="id")
  pls.p.pls <- left_join(pls.test.pls,pls.pred.pls,by="id")
  pls.pred.plsall <- rbind(pls.pred.plsall,pls.p.pls)
}

#plot ROC
library(pROC)
par(pty="s")
result.roc.pls <- roc(pls.pred.plsall$Label, 
                  pls.pred.plsall$Low,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.5, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, 
                  print.auc=TRUE,
                  legacy.axes=TRUE,
                  percent = TRUE,
                  xlab="False Positive Percentage",
                  ylab="True Positive Percentage",
                  col="black")
ci.se(result.roc.pls)
par(pty="m")

roc.pls <- roc(pls.pred.plsall$Label, 
             pls.pred.plsall$Low)

groc.pls <- ggroc(roc.pls,
                legacy.axes = TRUE,size=1)+
  theme_minimal()+
  
  theme(
    legend.position = c(0.5,0.5),
    panel.grid = element_blank(),
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(color = "black",size=16),
    axis.line = element_line(color="black"),
    axis.text = element_text(color="black",size=14),
    legend.text = element_text(color="black",size=16),
    legend.title = element_text(color="black",size=20)
  )+
  xlab("False Positive")+ylab("True Positive")+
  
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="grey", linetype="dashed",size=1.5)+
  geom_text(aes(x=0.6,y=0.4,label="AUC=0.873 (CI: 0.865-0.880)"),
            size=5,color="black")
groc.pls

ggsave(filename = "figures/ROC of PLS.svg",
       plot = groc.pls,
       width=12,
       height = 10,
       dpi = 75)
save(pls.pred.plsall,file = "proteome/PLS预测 样本为VIP大于2.Rdata")
save(pls.pred.tall,file = "proteome/PLS预测 样本为火山图结果.Rdata")
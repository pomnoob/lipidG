library(tidyverse)
library(mixOmics)
library(caret)
library(pROC)
# 导入标准化后的 T1 & T3 数据
pred.norm <- read.csv(file = "final/Metab_data_normalized_t_T1&T3.csv",stringsAsFactors = F)
pred.norm <- pred.norm %>%
  dplyr::select(-sample) %>%
  mutate(group=factor(group,labels = c("T1","T3")))

# All proteins
# 分割数据
trainMeats <- createDataPartition(pred.norm$group, p = 2/3,list = F)
# 准备数据用于pls-da
X <- as.matrix(pred.norm[-1])
Y <- pred.norm$group  
# 训练集和测试集
Xtrain <- X[ trainMeats, ]
Ytrain <- Y[ trainMeats ]
Xtest <- X[ -trainMeats, ]
Ytest <- Y[ -trainMeats ]

## PLS-DA function
plsda.train <- plsda(Xtrain, Ytrain, ncomp = 6) # where ncomp is the number of components wanted

# 预测
plsda.pred <- predict(plsda.train, Xtest, dist = "max.dist")

prediction <- plsda.pred$class$max.dist[,4] 

confusion.mat <- get.confusion_matrix(truth = Ytest, predicted = prediction)

#######################################
# all data

plsda.res <- plsda(X, Y, ncomp = 10)

set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 10, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 100) 
# perf.plsda.srbct$error.rate  # error rates
plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

auc.plsda = mixOmics::auroc(plsda.res, roc.comp = 2)

plotIndiv(plsda.res , comp = 1:2,
          group = pred.norm$group, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA on SRBCT')

# with background
background = background.predict(plsda.res, comp.predicted=2, dist = "max.dist") 
#optional: xlim = c(-40,40), ylim = c(-30,30))

plotIndiv(plsda.res, comp = 1:2,
          group = pred.norm$group, ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

p.vip <- vip(plsda.res)
########################################################################################
# spls-da
set.seed(2543) # for reproducibility
# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10, seq(20, 400, 10))

tune.splsda.maBER <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                    progressBar = TRUE, dist = 'mahalanobis.dist', measure = "BER",
                                    test.keepX = list.keepX, nrepeat = 10) # ncomp=1

tune.splsda.maAUC <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                 progressBar = TRUE, dist = 'mahalanobis.dist', measure = "AUC",
                                 test.keepX = list.keepX, nrepeat = 10) # ncomp=3

tune.splsda.maoverall <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                 progressBar = TRUE, dist = 'mahalanobis.dist', measure = "overall",
                                 test.keepX = list.keepX, nrepeat = 10) # ncomp=4 # BEST so far

tune.splsda.maxBER <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                     progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                     test.keepX = list.keepX, nrepeat = 10) # ncomp=3

tune.splsda.maxAUC <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "AUC",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=2
set.seed(2543)
tune.splsda.maxoverall <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "overall",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=1

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.maxoverall$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX <- tune.splsda.maxoverall$choice.keepX[1:ncomp]
plot(tune.splsda.maxoverall, col = color.jet(5))

splsda.srbct <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX)
plotIndiv(splsda.srbct, comp = c(1,2),
          group = Y, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')
auc.splsda = auroc(splsda.srbct, roc.comp = 2)

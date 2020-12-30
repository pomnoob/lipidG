library(tidyverse)
library(mixOmics)
library(caret)
library(pROC)
# 导入标准化后的 T1 to T3 数据
pred.all <- read.csv(file = "final/Metab_data_normalized_t.csv",stringsAsFactors = F)
pred.all <- pred.all %>%
  dplyr::select(-sample) %>%
  mutate(group=factor(group,labels = c("T1","T2","T3")))

# All proteins
# 分割数据
trainMeats <- createDataPartition(pred.norm$group, p = 2/3,list = F)
# 准备数据用于pls-da
X.all <- as.matrix(pred.all[-1])
Y.all <- pred.all$group  


#######################################
# all data

plsda.res.all <- plsda(X.all, Y.all, ncomp = 20)

set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda.all <- perf(plsda.res.all, validation = "Mfold", folds = 10, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 50) # 10 components is better
# perf.plsda.srbct$error.rate  # error rates
plot(perf.plsda.all, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

auc.plsda = mixOmics::auroc(plsda.res.all, roc.comp = 10)

plotIndiv(plsda.res.all , comp = 1:2,
          group = pred.all$group, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA on SRBCT')

# with background
background = background.predict(plsda.res.all, comp.predicted=2, dist = "max.dist") 
#optional: xlim = c(-40,40), ylim = c(-30,30))

plotIndiv(plsda.res.all, comp = 2:3,
          group = pred.all$group, ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

p.vip <- vip(plsda.res)
########################################################################################
# spls-da
# 根据PLSDA的结果，最佳ncomp=10, 方法最佳是 max.dist
# grid of possible keepX values that will be tested for each component


list.keepX <- c(1:10, seq(20, 100, 10))

### max.dist and BER
set.seed(25434564)
splsda.maxBER <- tune.splsda(X.all, Y.all, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=1

error.maxBER <- splsda.maxBER$error.rate
ncomp.maxBER <- splsda.maxBER$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maxBER <- splsda.maxBER$choice.keepX[1:ncomp.maxBER]
plot(splsda.maxBER, col = color.jet(10))

### mahalanobis.dist and BER
set.seed(25434564)
splsda.maBER <- tune.splsda(X.all, Y.all, ncomp = 10, validation = 'Mfold', folds = 10,
                             progressBar = TRUE, dist = 'mahalanobis.dist', measure = "BER",
                             test.keepX = list.keepX, nrepeat = 10) # ncomp=1

error.maBER.all <- splsda.maBER$error.rate
ncomp.maBER.all <- splsda.maBER$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maBER.all <- splsda.maBER$choice.keepX[1:ncomp.maBER]
plot(splsda.maBER, col = color.jet(10))

### max.dist and AUC
set.seed(25434564)
splsda.maxAUC <- tune.splsda(X.all, Y.all, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "AUC",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=3

error.maxAUC.all <- splsda.maxAUC$error.rate
ncomp.maxAUC.all <- splsda.maxAUC$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maxAUC.all <- splsda.maxAUC$choice.keepX[1:ncomp.maxAUC.all]
plot(splsda.maxAUC, col = color.jet(5))

### mahalanobis.dist and AUC
set.seed(25434564)
splsda.maAUC <- tune.splsda(X.all, Y.all, ncomp = 5, validation = 'Mfold', folds = 10,
                             progressBar = TRUE, dist = 'mahalanobis.dist', measure = "AUC",
                             test.keepX = list.keepX, nrepeat = 10) # ncomp=3

error.maAUC.all <- splsda.maAUC$error.rate
ncomp.maAUC.all <- splsda.maAUC$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maAUC.all <- splsda.maAUC$choice.keepX[1:ncomp.maAUC.all]
plot(splsda.maAUC, col = color.jet(5))

### centroids.dist and AUC
set.seed(25434564)
splsda.ceAUC <- tune.splsda(X.all, Y.all, ncomp = 5, validation = 'Mfold', folds = 10,
                            progressBar = TRUE, dist = 'centroids.dist', measure = "AUC",
                            test.keepX = list.keepX, nrepeat = 10) # ncomp=3

error.ceAUC.all <- splsda.ceAUC$error.rate
ncomp.ceAUC.all <- splsda.ceAUC$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.ceAUC.all <- splsda.ceAUC$choice.keepX[1:ncomp.ceAUC.all]
plot(splsda.ceAUC, col = color.jet(5))

### max.dist and overall
set.seed(25434564)
tune.splsda.maxoverall <- tune.splsda(X, Y, ncomp = 10, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "overall",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=1

error.maxoverall <- tune.splsda.maxoverall$error.rate
ncomp.maxoverall <- tune.splsda.maxoverall$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maxoverall <- tune.splsda.maxoverall$choice.keepX[1:ncomp.maxoverall]
plot(tune.splsda.maxoverall, col = color.jet(10))

########################################################################################################
set.seed(2543) 
tune.splsda.maAUC <- tune.splsda(X, Y, ncomp = 10, validation = 'Mfold', folds = 10,
                                 progressBar = TRUE, dist = 'mahalanobis.dist', measure = "AUC",
                                 test.keepX = list.keepX, nrepeat = 10) # ncomp=2
set.seed(2555)
tune.splsda.maoverall <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                 progressBar = TRUE, dist = 'mahalanobis.dist', measure = "overall",
                                 test.keepX = list.keepX, nrepeat = 10) # ncomp=1 

set.seed(45654654)
tune.splsda.maxAUC <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "AUC",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=2
set.seed(2543)
tune.splsda.maxoverall <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "overall",
                                  test.keepX = list.keepX, nrepeat = 10) # ncomp=1

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.maxAUC$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX <- tune.splsda.maxAUC$choice.keepX[1:ncomp]
plot(tune.splsda.maxAUC, col = color.jet(5))

splsda.srbct <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX)
plotIndiv(splsda.srbct, comp = c(1,2),
          group = Y, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')
auc.splsda = auroc(splsda.srbct, roc.comp = 2)

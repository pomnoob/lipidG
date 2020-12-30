library(tidyverse)
library(mixOmics)

# 导入标准化后的 T1 to T3 数据
pred.norm <- read.csv(file = "final/Metab_data_normalized_t_T1&T3.csv",stringsAsFactors = F)
pred.norm <- pred.norm %>%
  dplyr::select(-sample) %>%
  mutate(group=factor(group,labels = c("T1","T3")))


# 准备数据用于pls-da
X <- as.matrix(pred.norm[-1])
Y <- pred.norm$group  


#######################################
# all data

plsda.res <- plsda(X, Y, ncomp = 20)

set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda <- mixOmics::perf(plsda.res, validation = "Mfold", folds = 10, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 50) 

perf.plsda$error.rate

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

auc.plsda = mixOmics::auroc(plsda.res, roc.comp = 2)

plotIndiv(plsda.res , comp = 1:2,
          group = pred.norm$group, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA on SRBCT')

# with background
background = background.predict(plsda.res, comp.predicted=2, dist = "max.dist") 
#optional: xlim = c(-40,40), ylim = c(-30,30))

plotIndiv(plsda.res, comp = 1:2,
          group = pred.norm$group, ind.names = FALSE, title = "max.dist",
          legend = TRUE,  background = background)

p.vip <- vip(plsda.res)

########################################################################################
# spls-da
# 根据PLSDA的结果，最佳 ncomp=2 方法最佳是 centroids.dist
# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10, seq(20, 100, 10))
set.seed(97845613)
tune.splsda.maBER <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                    progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                    test.keepX = list.keepX, nrepeat = 50) # ncomp=3

error.maBER <- tune.splsda.maBER$error.rate
ncomp.maBER <- tune.splsda.maBER$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maBER <- tune.splsda.maBER$choice.keepX[1:ncomp.maBER] # 7 100 5
plot(tune.splsda.maBER, col = color.jet(5))



### max.dist and AUC
set.seed(25434564)
tune.splsda.maxAUC <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10,
                                  progressBar = TRUE, dist = 'max.dist', measure = "AUC",
                                  test.keepX = list.keepX, nrepeat = 50, auc = T) # ncomp=3

error.maxAUC <- tune.splsda.maxAUC$error.rate
ncomp.maxAUC <- tune.splsda.maxAUC$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX.maxAUC <- tune.splsda.maxAUC$choice.keepX[1:ncomp.maxAUC] # 20 50 6
plot(tune.splsda.maxAUC, col = color.jet(5))


# sPLS-DA
splsda <- splsda(X, Y, ncomp = ncomp.maxAUC, keepX = select.keepX.maxAUC)
plotIndiv(splsda, comp = c(1,2),
          group = Y, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(splsda, comp = c(1,3),
          group = Y, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 3')

plotIndiv(splsda, comp = c(2,3),
          group = Y, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 2 & 3')
# AUC
auc.splsda = auroc(splsda, roc.comp = 3)

# Performance
perf.splsda <- perf(splsda, validation = "Mfold", folds = 10, 
                    progressBar = FALSE, auc = TRUE, nrepeat = 50) 

perf.splsda$error.rate
plot(perf.splsda)
# break down of error rate per class is also insightful on the prediction
# of the model:
perf.splsda$error.rate.class
selectVar(splsda.srbct, comp = 1)$value
selectVar(splsda.srbct, comp = 3)$value
# stability
par(mfrow=c(1,3))
plot(perf.splsda$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.splsda$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)
plot(perf.splsda$features$stable[[3]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 3', las =2)

# 筛选出的重要蛋白质
# comp1
# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda, comp = 1)$name, 
                  names(perf.splsda$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.splsda$features$stable[[1]][ind.match])

splsda.feature.c1 <- data.frame(selectVar(splsda, comp = 1)$value, Freq)
write.csv(splsda.feature.c1,
          file = "final/SPLSDA comp1 筛选蛋白.csv",
          row.names = T)
# comp2
# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda, comp = 2)$name, 
                  names(perf.splsda$features$stable[[2]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.splsda$features$stable[[2]][ind.match])

splsda.feature.c2 <- data.frame(selectVar(splsda, comp = 2)$value, Freq)
write.csv(splsda.feature.c2,
          file = "final/SPLSDA comp2 筛选蛋白.csv",
          row.names = T)
# comp3
# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda, comp = 3)$name, 
                  names(perf.splsda$features$stable[[3]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.splsda$features$stable[[3]][ind.match])

splsda.feature.c3 <- data.frame(selectVar(splsda, comp = 3)$value, Freq)
write.csv(splsda.feature.c3,
          file = "final/SPLSDA comp3 筛选蛋白.csv",
          row.names = T)
################################################################################
# error rate per comp
# comp 1
error.comp1 <- as.data.frame(tune.splsda.maxAUC[["error.rate.all"]][["comp1"]])
error.comp2 <- as.data.frame(tune.splsda.maxAUC[["error.rate.all"]][["comp2"]])
error.comp3 <- as.data.frame(tune.splsda.maxAUC[["error.rate.all"]][["comp3"]])
write.csv(error.comp1,file = "final/splsda结果 comp1.csv",row.names = T)
write.csv(error.comp2,file = "final/splsda结果 comp2.csv",row.names = T)
write.csv(error.comp3,file = "final/splsda结果 comp3.csv",row.names = T)

## 在Excel 中整理数据，重新导入画图
error.all <- read.csv(file = "final/splsda error 结果汇总.csv")
ggplot(data = error.all)+
  geom_point(aes(x=as.factor(variable),y=mean),size=5)+
  geom_errorbar(aes(x=as.factor(variable), ymin=mean-sd,ymax=mean+sd),
                width=0.4, alpha=0.9)+
  facet_wrap(~class,scales="free",nrow = 1)

# 上面代码画出来的图不太好看，error bar相对y坐标太长

# 做boxplot图
# comp1
error.comp1 <- error.comp1 %>%
  rownames_to_column(var = "variable") %>%
  dplyr::select(1:21)

error.c1.melt <- melt(error.comp1,id.vars = "variable",
                      value.name = "AUC",
                      variable.name = "iter")
error.c1.melt$comp <- rep("Component 1",380)
# comp2
error.comp2 <- error.comp2 %>%
  rownames_to_column(var = "variable") %>%
  dplyr::select(1:21)

error.c2.melt <- melt(error.comp2,id.vars = "variable",
                      value.name = "AUC",
                      variable.name = "iter")
error.c2.melt$comp <- rep("Component 2",380)
# comp3
error.comp3 <- error.comp3 %>%
  rownames_to_column(var = "variable") %>%
  dplyr::select(1:21)

error.c3.melt <- melt(error.comp3,id.vars = "variable",
                      value.name = "AUC",
                      variable.name = "iter")
error.c3.melt$comp <- rep("Component 3",380)

# 合并数据
error.melt.all <- rbind(error.c1.melt,error.c2.melt,error.c3.melt)
error.melt.all$variable <- as.numeric(error.melt.all$variable)
error.melt.all <- error.melt.all %>%
  arrange(variable)
error.melt.all <- error.melt.all %>%
  filter(variable <=60)

error.melt.all <- error.melt.all %>%
  mutate(type=case_when(variable==20&comp=="Component 1"~"Highlight",
                        variable==30&comp=="Component 2"~"Highlight",
                        variable==6&comp=="Component 3"~"Highlight",
                        TRUE~"Normal"))


variable.sel <- ggplot(data = error.melt.all)+
  geom_boxplot(aes(x=as.factor(variable),y=AUC,fill=type))+
  scale_fill_manual(values=c("#69b3a2", "white")) +
  facet_wrap(~comp,nrow = 1)+
  theme_minimal()+
  theme(
    strip.text = element_text(size=15),
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(size=15),
    axis.text = element_text(size=15),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.position = "none"
  )+
  xlab("Number of Variables")+
  ylab("Aera Under the Curve (AUC)")
variable.sel

ggsave(filename = "final/splsda 变量AUC图.svg",
       variable.sel,
       width=14,
       height=12,
       dpi=75)

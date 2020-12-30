# 根据splsda的结果，重新作图查看各个蛋白质浓度的差异
# 准备数据
prot.orig <- read.csv(file = "final/Metab_data_processed_t.csv",stringsAsFactors = F)


# 查找重复ID
colna  <- colnames(prot.orig[3:426])
ccc <- character()

for (i in 1:length(colna )) {
  
  if (colna [i] %in% colna [-i]) {
    ccc[length(ccc)+1] <- colna [i]
  }
  
}

# 均值和SD
mean <- aggregate(prot.orig[3:426],by=list(ter=prot.orig$group),mean,na.rm=T)
sd <- aggregate(prot.orig[3:426],by=list(ter=prot.orig$group),sd,na.rm=T)

library(reshape2)
mean.m <- melt(mean,id.vars = "ter",
               variable.name = "feature",
               value.name = "mean")

sd.m <- melt(sd,id.vars = "ter",
             variable.name = "feature",
             value.name = "sd")

mean.m$sd <- sd.m$sd

# 导入之前做好的蛋白质ID转Gene ID的文件和代码
pro.id <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)

# refine the names
pro.id$gname <- toupper(pro.id$gname)
pro.id$pname <- str_replace_all(pro.id$pname,
                                "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                "")
pro.id$pname <- str_to_title(pro.id$pname)

# 查找重复ID
#data %>%
#  group_by(id1, id2) %>%
#  filter(row_number() == 1) %>%
#  ungroup()

pro.id <- pro.id %>%
  group_by(feature)%>%
  filter(row_number()==1)%>%
  ungroup

forPlot <- left_join(mean.m,pro.id,by="feature")

################################################################################
# 选择splsda模型中comp1的蛋白质作图
splsda.c1 <- read.csv(file = "final/SPLSDA comp1 筛选蛋白.csv",stringsAsFactor=F)%>%
  rename(feature=X)
c1.feature <- splsda.c1$feature
forPlot.c1 <- forPlot %>%
  filter(feature %in% c1.feature & ter != 2)
# 画图
plot.c1 <- ggplot(data=forPlot.c1)+
  geom_bar(aes(x=as.factor(ter),y=mean,fill=as.factor(ter)),stat = "identity",width = 0.5,color="black")+
  geom_errorbar(aes(x=as.factor(ter), ymin=mean,ymax=mean+sd),
                width=0.4, alpha=0.9)+
  facet_wrap(~feature,scales = "free",nrow=4)+
  scale_fill_grey(start = 0.4,end = 0.9,labels=c("Tertile 1", "Tertile 3"))+theme_minimal()+
  theme(
    strip.text = element_text(size=15),
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(size=15),
    axis.text = element_text(size=15),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.position = "bottom"
  )+
  xlab("Tertile of WAZ")+
  ylab("One Ten Thousandth of Total Protein")+
  labs(fill="Tertile of WAZ ")
plot.c1

# 选择splsda模型中comp2的蛋白质作图
splsda.c2 <- read.csv(file = "final/SPLSDA comp2 筛选蛋白.csv",stringsAsFactor=F)%>%
  rename(feature=X)%>%
  slice(1:30)
c2.feature <- splsda.c2$feature
forPlot.c2 <- forPlot %>%
  filter(feature %in% c2.feature & ter != 2)
# 画图
plot.c2 <- ggplot(data=forPlot.c2)+
  geom_bar(aes(x=as.factor(ter),y=mean,fill=as.factor(ter)),stat = "identity",width = 0.5,color="black")+
  geom_errorbar(aes(x=as.factor(ter), ymin=mean,ymax=mean+sd),
                width=0.4, alpha=0.9)+
  facet_wrap(~feature,scales = "free",nrow=5)+
  scale_fill_grey(start = 0.4,end = 0.9,labels=c("Tertile 1", "Tertile 3"))+theme_minimal()+
  theme(
    strip.text = element_text(size=15),
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(size=15),
    axis.text = element_text(size=15),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.position = "bottom"
  )+
  xlab("Tertile of WAZ")+
  ylab("One Ten Thousandth of Total Protein")+
  labs(fill="Tertile of WAZ ")
plot.c2

# 选择splsda模型中comp2的蛋白质作图
splsda.c3 <- read.csv(file = "final/SPLSDA comp3 筛选蛋白.csv",stringsAsFactor=F)%>%
  rename(feature=X)
c3.feature <- splsda.c3$feature
forPlot.c3 <- forPlot %>%
  filter(feature %in% c3.feature & ter != 2)
# 画图
plot.c3 <- ggplot(data=forPlot.c3)+
  geom_bar(aes(x=as.factor(ter),y=mean,fill=as.factor(ter)),stat = "identity",width = 0.5,color="black")+
  geom_errorbar(aes(x=as.factor(ter), ymin=mean,ymax=mean+sd),
                width=0.4, alpha=0.9)+
  facet_wrap(~feature,scales = "free",nrow=2)+
  scale_fill_grey(start = 0.4,end = 0.9,labels=c("Tertile 1", "Tertile 3"))+theme_minimal()+
  theme(
    strip.text = element_text(size=15),
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(size=15),
    axis.text = element_text(size=15),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.position = "bottom"
  )+
  xlab("Tertile of WAZ")+
  ylab("One Ten Thousandth of Total Protein")+
  labs(fill="Tertile of WAZ ")
plot.c3
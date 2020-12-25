# 2020.12.24 
# 方法最终定下来：
# 使用 Normalyzer 包对比各类 normalization 方法，结果是 quntile normalization 最优
# 选择的数据前处理方法是 Quantile Normalization + Cube root transformation + range scaling
# 数据前处理使用 MetaboAnalyst 进行
# /final 文件夹下 Metab_ 开头的文件为 MetaboAnalyst 服务器下载的处理后的数据

################################################################################
################################################################################
################################################################################
library(tidyverse)
# 先查看 twaz 各组间的人口学指标有何差异
load(file = "final/demographic across tertile of waz.Rdata")

# 整体均值
sapply(pomic.char[continousV], mean, na.rm=T)
sapply(pomic.char[continousV], sd, na.rm=T)
# 分组均值
continousV <- c("waz","ageBaby","ageMon","birthWeight","birthLength","monBMI")
char.mean <- aggregate(pomic.char[continousV],by=list(ter=pomic.char$twaz),mean,na.rm=T)
char.sd <- aggregate(pomic.char[continousV],by=list(ter=pomic.char$twaz),sd,na.rm=T)

write.csv(char.mean,file = "final/表一数据 特征平均值.csv",
          row.names = F)
write.csv(char.sd,file = "final/表一数据 特征标准差.csv",
          row.names = F)

# twaz 之间做ANOVA

fit_waz <- lm(waz~twaz,data = pomic.char)
summary(fit_waz) #2.2e-16

fit_ageBaby <- lm(ageBaby~twaz,data = pomic.char)
summary(fit_ageBaby) #p=0.002

fit_ageMon <- lm(ageMon~twaz,data = pomic.char)
summary(fit_ageMon) #p=0.128

fit_birthWeight <- lm(birthWeight~twaz,data = pomic.char)
summary(fit_birthWeight) #p=0.0002

fit_birthLength <- lm(birthLength~twaz,data = pomic.char)
summary(fit_birthLength) #p=0.086

fit_monBMI <- lm(monBMI~twaz,data = pomic.char)
summary(fit_monBMI) #p=0.366

# 查看离散型变量
# sex
table(pomic.char$sex)
prop.table(table(pomic.char$sex))
chi.sex <- xtabs(~sex+twaz,data = pomic.char)
chi.sex
prop.table(chi.sex,2)
chisq.test(chi.sex) # p=0.965
# city
table(pomic.char$city)
prop.table(table(pomic.char$city))
chi.city <- xtabs(~city+twaz,data = pomic.char)
chi.city
prop.table(chi.city,2)
chisq.test(chi.city) # p=0.016

# parity
table(pomic.char$parity)
prop.table(table(pomic.char$parity))
chi.parity <- xtabs(~parity+twaz,data = pomic.char)
chi.parity
prop.table(chi.parity,2)
chisq.test(chi.parity) # p=0.772

# education
table(pomic.char$edu)
prop.table(table(pomic.char$edu))
chi.edu <- xtabs(~edu+twaz,data = pomic.char)
chi.edu
prop.table(chi.edu,2)
chisq.test(chi.edu) # p=0.00042

################################################################################

# 对各个蛋白质做linear model
# 混杂因子：ageBaby, birthWeight, Education, City
# 使用的数据：标准化以后的蛋白质数据

# 合并数据
prot_norm <- read.csv("final/Metab_data_normalized_t.csv",
                      stringsAsFactors = F)

prot_norm <- prot_norm %>%
  dplyr::rename(id=sample)

prot.all <- left_join(prot_norm,pomic.char,by="id")

# 准备回归模型的相关参数

# 蛋白质id
prot <- prot_norm %>%
  select(-id,-group)
# 组别 
group <- prot_norm %>%
  select(group)
# 混杂因素
conf <- pomic.char %>%
  select(ageBaby,birthWeight,edu,city)

# 提取名字
prot.name <- colnames(prot)
group.name <- colnames(group)
conf.name <- colnames(conf)


# 线性模型

crntRcor.fc <- double()
crntPcor.fc <- double()
  
for (j in 1:length(prot.name)) {
  y.fc <- prot.name[j]
  x.fc <- group.name 
  cov_f <- paste(conf.name,collapse = "+")
  yhs.fc <- paste(y.fc,"~")
  xhs.fc <- paste(x.fc,"+")
    
  frma.fc <- as.formula(paste(yhs.fc,xhs.fc,cov_f))
  mod.fc <- lm(frma.fc,data = prot.all,na.action = na.exclude)
  coef.fc <- coef(mod.fc)
  names(coef.fc) <- NULL
  p.fc <- anova(mod.fc)
    
  crntRcor.fc[j] <- coef.fc[2]
  crntPcor.fc[j] <- p.fc[1,5]
    
  }
  

# 合并数据
lm.result <- data.frame(feature=prot.name,
                        p=crntPcor.fc,
                        r=crntRcor.fc)

# 校正p值
lm.result$fdr <- p.adjust(lm.result$p,method = "fdr")

# 提取fdr<0.05的数据
lm.result.sig <- lm.result %>%
  filter(fdr<0.05)

# 保存数据
write.csv(lm.result,file = "final/线性回归结果.csv",row.names=F)

# 对比一下VIP>2的结果
pls.vip <- read.csv(file = "final/Metab_plsda_vip.csv")
pls.vip <- pls.vip %>%
  select(feature=X,vip=Comp..1)

prot.vip.sig <- left_join(pls.vip,lm.result,by="feature")

prot.vip.sig.cut <- prot.vip.sig %>%
  filter(vip > 2 | fdr < 0.05)

# 导入之前做好的蛋白质ID转Gene ID的文件和代码
pro.id <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)

# refine the names
pro.id$gname <- toupper(pro.id$gname)
pro.id$pname <- str_replace_all(pro.id$pname,
                                "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                "")
pro.id$pname <- str_to_title(pro.id$pname)

prot.vip.sig.cut <- left_join(prot.vip.sig.cut,pro.id,by="feature")

write.csv(prot.vip.sig.cut,
          file = "final/VIP大于2或FDR小于0.05的蛋白质.csv",
          row.names = F)

################################################################################
################################################################################
# 根据原始浓度数据绘制图
# 需要绘制图片和列出的蛋白质
feature <- prot.vip.sig.cut$feature

# 准备数据
prot.orig <- read.csv(file = "final/Metab_data_processed_t.csv",stringsAsFactors = F)
# 均值和SD
mean <- aggregate(prot.orig[feature],by=list(ter=prot.orig$group),mean,na.rm=T)
sd <- aggregate(prot.orig[feature],by=list(ter=prot.orig$group),sd,na.rm=T)

write.csv(mean,file = "final/差异蛋白平均值.csv",row.names=F)

write.csv(sd,file = "final/差异蛋白标准差.csv",row.names=F)


################################################################################
# 在Excel 中调整数据格式，画柱状图
data.barp <- read.csv(file = "final/画图数据.csv",stringsAsFactors = F)
conc.protein <- ggplot(data=data.barp)+
  geom_bar(aes(x=as.factor(group),y=mean,fill=as.factor(group)),stat = "identity",width = 0.5,color="black")+
  geom_errorbar(aes(x=as.factor(group), ymin=mean,ymax=mean+sd),
                width=0.4, alpha=0.9)+
  facet_wrap(~feature,scales = "free",nrow=7)+
  scale_fill_grey(start = 0.4,end = 0.9,labels=c("Tertile 1", "Tertile 2", "Tertile 3"))+theme_minimal()+
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

ggsave(filename = "final/蛋白质浓度图.svg",
       conc.protein,
       width = 8,
       height = 20,
       dpi = 75)
  

###############################################################################
prot.orig.plot <- prot.orig %>%
  select(id,group,all_of(feature))

# 转置数据
library(reshape2)
prot.plot <- melt(prot.orig.plot,id.vars = c("id","group"),
                      variable.name = "feature",
                      value.name = "value")



# 作图效果不太好
ggplot(data=prot.plot)+
  geom_boxplot(aes(x=as.factor(group),y=value,
                   fill=as.factor(group)),width=0.5)+
  scale_fill_grey(start = 0.4,end = 1,
                  labels=c("Tertile 1", "Tertile 2", "Tertile 3"))+

  facet_wrap(~ feature, nrow = 5)+
  theme_minimal()+
  theme(
    panel.border  = element_rect(color = "black",
                                 fill=NA),
    axis.title = element_text(size=24),
    axis.text = element_text(size=20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20),
    legend.position = "top"
  )+
  xlab("Tertile")+
  ylab("Concentration")+
  labs(fill="Tertile of WAZ ")+
  scale_x_discrete(breaks=c(1,2,3),labels=c("T1","T2","T3"))+
  scale_y_continuous(limits = c(0,0.015),
                     breaks = c(0,0.0025,0.005,0.0075,0.010,0.015))

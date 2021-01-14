# 如果需要发好一些的杂志，需要考虑最终的验证问题
# 郑州和哈尔滨的样本没有进行过分析，可以使用这两个城市的样本进行验证工作
# 先在母乳数据库交接文件中整理出哈尔滨和郑州的蛋白组数据

# 提取其他6个城市的蛋白组名称情况
library(tidyverse)
pomic_all <- read.csv("data/pomics_all_t.csv",stringsAsFactors = F)
pomic_all <- pomic_all  %>%
  select(-group) %>%
  rename(id=sample)

pid6 <- colnames(pomic_all)
pid6 <- pid6[-1]

# 郑州数据
zz <- read.csv(file = "final/郑州蛋白组数据.csv")
zz_pid <- colnames(zz)
zz_pid <- zz_pid[-1]

# 哈尔滨数据
hrb <- read.csv(file = "final/哈尔滨蛋白组数据.csv")
hrb_pid <- colnames(hrb)
hrb_pid <- hrb_pid[-1]

# 查看在郑州和哈尔滨都检测到的蛋白质
allid <- character()

for (i in 1:length(hrb_pid)) {
  if (hrb_pid[i] %in% zz_pid) {
    allid[length(allid)+1] <- hrb_pid[i]
  } #两个城市都检测到的蛋白质只有364个
}

# 查看在郑州和哈尔滨都检测到的蛋白质在其他6个城市中是否都检测到
allid2 <- character()

for (i in 1:length(allid)) {
  if (allid[i] %in% pid6) {
    allid2[length(allid2)+1] <- allid[i]
  }
} # 结果是只有264个蛋白质重合

# 查看通过ANOVA和VIP删选的蛋白质
vip <- read.csv(file = "final/VIP大于2且FDR小于0.05的蛋白质.csv")
vipid <- vip$feature

# 查看通过ANOVA和VIP删选的蛋白质是否都在郑州和哈尔滨检测到
vipall <- character()

for (i in 1:length(vipid)) {
  if (vipid[i] %in% zz_pid) {
    vipall[length(vipall)+1] <- vipid[i]
  } #两个城市都检测到的蛋白质只有364个
}

#
splsda <- read.csv(file = "final/splsda comp123 蛋白质.csv")
splsdaid <- splsda$X

# 查看通过SPLSDA筛选的蛋白质是否都在郑州和哈尔滨检测到
splsall <- character()

for (i in 1:length(splsdaid)) {
  if (splsdaid[i] %in% zz_pid) {
    splsall[length(splsall)+1] <- splsdaid[i]
  } #两个城市都检测到的蛋白质只有364个
}

# 都含有的蛋白质进行整理
pomic_zz <- zz %>%
  select(sample, all_of(allid))

pomic_hrb <- hrb %>%
  select(sample, all_of(allid))

pomic_zh <- rbind(pomic_zz,pomic_hrb) %>%
  rename(id = sample)

################################################################################
################################################################################
# 研究思路：先合并郑州和哈尔滨的数据，用ANOVA和VIP筛选的蛋白质进行预测，
# 其中12个蛋白质中，郑州和哈尔滨都检测到的有10种
# 经过同样的标准化过程后，查看这10种蛋白质的含量有何差异
# 然后用pls模型进行预测，查看准确度

# 导入问卷数据
ques <- read.csv(file = "data/Metadata of breastmilk questionnaire.csv",stringsAsFactors = F)

# 根据问卷数据计算母亲、婴儿年龄等数据
ques$birthMon <- as.Date(paste(ques$A2A,ques$A2B,
                               ques$A2C,sep="-"))#乳母出生日期
ques$birthBaby<- as.Date(paste(ques$A5A,
                               ques$A5B,ques$A5C,sep="-"))
ques$sampleDate <- as.Date(ques$N1,
                           format="%d-%m-%Y")#样本采集日期
# 样本采集时母亲和婴儿年龄
ques <- ques %>%
  mutate(ageMon=(sampleDate-birthMon)/365,
         ageBaby=(sampleDate-birthBaby))
ques$ageBaby <- as.numeric(ques$ageBaby)
ques$ageMon <- as.numeric(ques$ageMon)
ques$id <- as.numeric(ques$id)

# 列出一下问卷数据编号代表的内容
# A4:infant gender; TB3: infant body weight; TB1:infant length
# TB6: infant head circumference; A3: parity; M1: maternal education
# M3: family annual income; M4: prepregnant maternal body weight
# M5: maternal body weight postpartum; M7: delivery mode
# TA1: maternal height; TA3: maternal body weight on site; B401: infant allergy; B5: maternal allergy history
# B6: father allergy history;B1: birth weight;B2:birth length;B3: preterm or not

quesSel <- ques %>%
  select(id, city, A4, TB3, TB1, TB6, A3, M1, M3, M4, 
         M5, M7, TA1, B401, B5, B6, ageBaby, ageMon,B1,B2,B3,TA3,B12) %>%
  rename(sex=A4,babyWeight=TB3,babyLength=TB1,
         babyHead=TB6,parity=A3,edu=M1, income=M3,
         preMonWeight=M4,postMonWeight=M5,delivery=M7,monHeight=TA1,
         monWeight=TA3,allergy=B401,MonAllergy=B5,FatAllergy=B6,
         birthWeight=B1,birthLength=B2,preterm=B3,mix=B12)

pomic <- inner_join(pomic_zh,quesSel,by="id")

# 婴儿数据
# 查看是否有身高、体重异常的样本
boxplot(pomic$babyWeight) 

# 身高无异常值
boxplot(pomic$babyLength)

# 母亲数据
summary(pomic$monHeight)
boxplot(pomic$monHeight)

# 查看母亲体重分布情况，去除异常值
summary(pomic$preMonWeight)
boxplot(pomic$preMonWeight)
pomic$preMonWeight[pomic$preMonWeight>=120] <- NA
# 产后体重
summary(pomic$postMonWeight)
boxplot(pomic$postMonWeight)

# 采样体重
summary(pomic$monWeight)
boxplot(pomic$monWeight)
pomic$monBMI <- pomic$monWeight/(pomic$monHeight/100)^2
summary(pomic$monBMI)


# 计算 z 评分
library(zscorer)
pomic <- addWGSR(data = pomic, sex = "sex", firstPart = "babyWeight",
                 secondPart = "ageBaby", index = "wfa",output = "waz")

pomic <- addWGSR(data = pomic, sex = "sex", firstPart = "babyLength",
                 secondPart = "ageBaby", index = "hfa",output = "haz")

pomic <- pomic %>%
  mutate(nBMI=case_when(monBMI>27.5~1,
                        monBMI>=18.5 & monBMI<=23~2,
                        monBMI>23 & monBMI<=27.5~3))
table(pomic$nBMI)

# 只排除早产儿
pomic.p <- pomic %>%
  filter(preterm==1)

# 去除混合喂养的
pomic.p$mix[is.na(pomic.p$mix)] <- 1
table(pomic.p$mix)
pomic.pm <- pomic.p %>%
  filter(mix != 2)

pomic.pm$twaz <- ntile(pomic.pm$waz,3)

pomic.expm <- pomic.pm %>%
  select(id,twaz,2:365) %>%
  rename(group=twaz)
pomic.expm[pomic.expm==0] <- NA

write.csv(pomic.expm,file = "final/WAZ T1toT3 all proteins no preterm not mixed feeding-zz and hrb.csv",row.names = F)

################################################################################

# 对各个蛋白质做linear model
# 混杂因子：ageBaby, birthWeight, Education, City
# 使用的数据：标准化以后的蛋白质数据

# 合并数据
prot_norm <- read.csv("final/Metabo_data_normalized_zh_t.csv",
                      stringsAsFactors = F)

pomic.char <- pomic.pm %>%
  select(id,366:392)

prot.all <- inner_join(prot_norm,pomic.char,by="id")


################################################################################

################################################################################
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

# 导入之前做好的蛋白质ID转Gene ID的文件和代码
pro.id <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)

# refine the names
pro.id$gname <- toupper(pro.id$gname)
pro.id$pname <- str_replace_all(pro.id$pname,
                                "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                "")
pro.id$pname <- str_to_title(pro.id$pname)

# 加入gene id

lm.result <- inner_join(lm.result,pro.id,by="feature")


# 提取fdr<0.05的数据
lm.result.sig <- lm.result %>%
  filter(fdr<0.05)

# 保存数据
write.csv(lm.result,file = "final/线性回归结果.csv",row.names=F)
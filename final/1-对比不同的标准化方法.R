library(tidyverse)
library(zscorer)
# 导入之前整理好的脂质组学数据
# 仍然没有哈尔滨和郑州的数据
# 因为这两个城市的样本id与检测数据没法一一对应
# 选择所有蛋白质
pomic_all <- read.csv("data/pomics_all_t.csv",stringsAsFactors = F)
pomic_all <- pomic_all  %>%
  select(-group) %>%
  rename(id=sample)

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


# 合并脂质组学和问卷数据
pomic <- inner_join(pomic_all,quesSel,by="id")
write.csv(pomic,file = "proteome for growth.csv",
          row.names = F)

# 婴儿数据
# 查看是否有身高、体重异常的样本
boxplot(pomic$babyWeight) 
# 有两个体重大于60的样本，设为NA
pomic$babyWeight[pomic$babyWeight > 60] <- NA
# 身高无异常值
boxplot(pomic$babyLength)

# 母亲数据
summary(pomic$monHeight)
boxplot(pomic$monHeight)
pomic$monHeight[pomic$monHeight<100] <- NA
# 查看母亲体重分布情况，去除异常值
summary(pomic$preMonWeight)
boxplot(pomic$preMonWeight)
pomic$preMonWeight[pomic$preMonWeight>=120] <- NA
# 产后体重
summary(pomic$postMonWeight)
boxplot(pomic$postMonWeight)
pomic$postMonWeight[pomic$postMonWeight>=120] <- NA
# 采样体重
summary(pomic$monWeight)
boxplot(pomic$monWeight)
pomic$monBMI <- pomic$monWeight/(pomic$monHeight/100)^2
summary(pomic$monBMI)

# 计算 z 评分
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

write.csv(pomic.pm, 
          file = "final/proteome data with z score.csv",
          row.names = F)

# 2021年1月14日，分析哈尔滨和郑州的样本发现重复性不太好
# 可能是北京的样本有较大偏差，去掉北京的样本进行分析

pomic.pm.exBJ <- pomic.pm %>%
  filter(city!="北京")

pomic.pm.exBJ$twaz <- ntile(pomic.pm.exBJ$waz,3)
pomic.pmpBJ <- pomic.pm.exBJ %>%
  filter(!is.na(twaz))

pomic.pm$twaz <- ntile(pomic.pm$waz,3)
pomic.pmp <- pomic.pm %>%
  filter(!is.na(twaz))
  
pomic.expm <- pomic.pmp %>%
  select(id,twaz,2:473) %>%
  rename(group=twaz)

pomic.expmBJ <- pomic.pmpBJ %>%
  select(id,twaz,2:473) %>%
  rename(group=twaz)



write.csv(pomic.expm,file = "final/WAZ T1toT3 all proteins no preterm not mixed feeding.csv",row.names = F)
write.csv(pomic.expmBJ,file = "final/WAZ T1toT3 all proteins no preterm not mixed feeding excl BJ.csv",row.names = F)

pomic.char <- pomic.pmp %>%
  select(id,475:501)
save(pomic.char,file = "final/demographic across tertile of waz.Rdata")

# 使用 Normalyzer 进行标准化
library(Normalyzer)
library(grid)
Normalyzer::normalyzer("NormalyzerDE/data for normalyzer.txt",
                       "proteome_norm3") # 无法生成 PDF


# 使用 NormalyzerDE
library(NormalyzerDE)
data_mat <- read.csv(file = "NormalyzerDE/Data matrix.csv",stringsAsFactors = F,header = F)
design_mat <- read.csv(file = "NormalyzerDE/Design matrix.csv",stringsAsFactors = F,header = F)
# 导出为txt文件
write.table(data_mat,file = "NormalyzerDE/Data matrix.txt",sep = "\t",quote = F,col.names = F,row.names = F)
write.table(design_mat,file = "NormalyzerDE/Design matrix.txt",sep = "\t",quote = F,col.names = F,row.names = F)

# 标准化数据
normalyzer("proteome_norm2",
             designPath = "NormalyzerDE/Design matrix.txt",
             dataPath = "NormalyzerDE/Data matrix.txt",
             outputDir = "NormalyzerDE/output") # 有数值低于1，无法进行后续分析


normalyzerDE("proteome_norm2",
                         designPath = "NormalyzerDE/Design matrix.txt",
                         dataPath = "NormalyzerDE/Data matrix.txt",
                         outputDir = "NormalyzerDE/output",
             comparisons = c("1-2", "1-3"))


# 使用MetaboAnalyst进行Cube root transformation 后使用 NormalyzerDE
library(NormalyzerDE)
data_mat.c <- read.csv(file = "NormalyzerDE/Data matrix cubic norm.csv",stringsAsFactors = F,header = F)
design_mat.c <- read.csv(file = "NormalyzerDE/Design matrix cubic norm.csv",stringsAsFactors = F,header = F)
# 导出为txt文件
write.table(data_mat.c,file = "NormalyzerDE/Data matrix cubic norm.txt",sep = "\t",quote = F,col.names = F,row.names = F)
write.table(design_mat.c,file = "NormalyzerDE/Design matrix cubic norm.txt",sep = "\t",quote = F,col.names = F,row.names = F)

# 标准化数据
normalyzer("proteome_norm2",
           designPath = "NormalyzerDE/Design matrix cubic norm.txt",
           dataPath = "NormalyzerDE/Data matrix cubic norm.txt",
           outputDir = "NormalyzerDE/output",
           noLogTransform = T) # 有数值低于1，无法进行后续分析


normalyzerDE("proteome_norm2",
             designPath = "NormalyzerDE/Design matrix.txt",
             dataPath = "NormalyzerDE/Data matrix.txt",
             outputDir = "NormalyzerDE/output",
             comparisons = c("1-2", "1-3"))


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
# 筛选样本
pomicSel <- pomic %>%
  # 选取月龄不小于3且非早产
  filter(ageBaby > 60 & preterm == 1)

# 只排除早产儿
pomic.p <- pomic %>%
  filter(preterm==1)



# 三分位
pomic.p$twaz <- ntile(pomic.p$waz,3)
table(pomic.p$twaz)


# 策略1
pomic.p <- pomic.p %>%
  mutate(nwaz=case_when(waz > 1~1,
                        waz < -1~3,
                        waz <= 1 & waz >= -1~2))
table(pomic.p$nwaz)
# 选出发育迟缓和正常的数据
waz1 <- pomic.p %>%
  filter(nwaz == 1 | nwaz == 2)
# 去除混合喂养的
pomic.p$mix[is.na(pomic.p$mix)] <- 1
table(pomic.p$mix)
pomic.pm <- pomic.p %>%
  filter(mix != 2)

waz1$mix[is.na(waz1$mix)] <- 1
table(waz1$mix)
waz1.m <- waz1 %>%
  filter(mix != 2)

table(waz1.m$nwaz)

# 导出1 和-1
pomic.expm <- pomic.pm %>%
  select(id,nwaz,2:473) %>%
  rename(group=nwaz)

# 导出三分位
pomic.expmt <- pomic.pm %>%
  select(id,twaz,2:473) %>%
  rename(group=twaz)

# 导出三分位T1heT3
wazT <- pomic.expmt %>%
  filter(group ==1 | group ==3)

pomic.exp <- pomic.p %>%
  select(id,nwaz,2:473) %>%
  rename(group=nwaz)

write.csv(pomic.expmt,file = "data/WAZ tertile all proteins no preterm not mixed feeding.csv",row.names = F)
write.csv(wazT,file = "data/WAZ tertile 1&3 all proteins no preterm not mixed feeding.csv",row.names = F)

write.csv(pomic.expm,file = "data/WAZ all proteins no preterm not mixed feeding.csv",row.names = F)
write.csv(pomic.exp,file = "data/WAZ all proteins no preterm.csv",row.names = F)
####### 过敏
table(pomic.p$allergy)
pomic.p$allergy <- ifelse(pomic.p$allergy==1,0,1)
p.allergy <- pomic.p %>%
  select(id,allergy,2:473) %>%
  rename(group=allergy)

write.csv(p.allergy,file = "data/allergy all proteins no preterm 1 and -1.csv",row.names = F)

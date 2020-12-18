library(tidyverse)
library(zscorer)
# 导入之前整理好的脂质组学数据
# 仍然没有哈尔滨和郑州的数据
# 因为这两个城市的样本id与检测数据没法一一对应
# 选择所有蛋白质
glycome <- read.csv("glycome/gomics.csv",stringsAsFactors = F)
glycome <- glycome  %>%
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
         M5, M7, TA1, B401, B5, B6, ageBaby, ageMon,B1,B2,B3,TA3) %>%
  rename(sex=A4,babyWeight=TB3,babyLength=TB1,
         babyHead=TB6,parity=A3,edu=M1, income=M3,
         preMonWeight=M4,postMonWeight=M5,delivery=M7,monHeight=TA1,
         monWeight=TA3,allergy=B401,MonAllergy=B5,FatAllergy=B6,
         birthWeight=B1,birthLength=B2,preterm=B3)


# 合并脂质组学和问卷数据
glycome <- inner_join(glycome,quesSel,by="id")


# 婴儿数据
# 查看是否有身高、体重异常的样本
boxplot(glycome$babyWeight) 
# 有两个体重大于60的样本，设为NA
glycome$babyWeight[glycome$babyWeight > 60] <- NA
# 身高无异常值
boxplot(glycome$babyLength)

# 母亲数据
summary(glycome$monHeight)
boxplot(glycome$monHeight)
glycome$monHeight[glycome$monHeight<100] <- NA
# 查看母亲体重分布情况，去除异常值
summary(glycome$preMonWeight)
boxplot(glycome$preMonWeight)
glycome$preMonWeight[glycome$preMonWeight>=120] <- NA
# 产后体重
summary(glycome$postMonWeight)
boxplot(glycome$postMonWeight)
glycome$postMonWeight[glycome$postMonWeight>=120] <- NA
# 采样体重
summary(glycome$monWeight)
boxplot(glycome$monWeight)
glycome$monBMI <- glycome$monWeight/(glycome$monHeight/100)^2
summary(glycome$monBMI)

# 计算 z 评分
glycome <- addWGSR(data = glycome, sex = "sex", firstPart = "babyWeight",
                 secondPart = "ageBaby", index = "wfa",output = "waz")

glycome <- addWGSR(data = glycome, sex = "sex", firstPart = "babyLength",
                 secondPart = "ageBaby", index = "hfa",output = "haz")

glycome <- glycome %>%
  mutate(nBMI=case_when(monBMI>27.5~1,
                        monBMI>=18.5 & monBMI<=23~2,
                        monBMI>23 & monBMI<=27.5~3))
table(glycome$nBMI)

# 只排除早产儿
glycome.p <- glycome %>%
  filter(preterm==1)

# 策略1
glycome.p <- glycome.p %>%
  mutate(nwaz=case_when(waz > 1~1,
                        waz < -1~2))
table(glycome.p$nwaz)
# 选出发育迟缓和正常的数据
waz1 <- glycome.p %>%
  filter(nwaz == 1 | nwaz == 2)

waz1.exp <- waz1 %>%
  select(id,nwaz,2:13) %>%
  rename(group=nwaz)

write.csv(waz1.exp,file = "glycome/WAZ HMO no preterm 1 and -1.csv",row.names = F)

####### 过敏
table(glycome.p$allergy)
glycome.p$allergy <- ifelse(glycome.p$allergy==1,0,1)
p.allergy <- glycome.p %>%
  select(id,allergy,2:13) %>%
  rename(group=allergy)

write.csv(p.allergy,file = "glycome/allergy HMO no preterm 1 and -1.csv",row.names = F)

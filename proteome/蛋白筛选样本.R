library(tidyverse)
library(zscorer)
# 导入之前整理好的脂质组学数据
# 仍然没有哈尔滨和郑州的数据
# 因为这两个城市的样本id与检测数据没法一一对应
pomic_all <- read.csv("data/pomics_20.csv",stringsAsFactors = F)
pomic_all <- pomic_all  %>%
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
# TA1: maternal height; B401: infant allergy; B5: maternal allergy history
# B6: father allergy history;B1: birth weight;B2:birth length;B3: preterm or not

quesSel <- ques %>%
  select(id, city, A4, TB3, TB1, TB6, A3, M1, M3, M4, 
         M5, M7, TA1, B401, B5, B6, ageBaby, ageMon,B1,B2,B3) %>%
  rename(sex=A4,babyWeight=TB3,babyLength=TB1,
         babyHead=TB6,parity=A3,edu=M1, income=M3,
         preMonWeight=M4,postMonWeight=M5,delivery=M7,monHeight=TA1,
         allergy=B401,MonAllergy=B5,FatAllergy=B6,
         birthWeight=B1,birthLength=B2,preterm=B3)

# 合并脂质组学和问卷数据
pomic <- inner_join(pomic_all,quesSel,by="id")

# 查看是否有身高、体重异常的样本
boxplot(pomic$babyWeight) 
# 有两个体重大于60的样本，设为NA
pomic$babyWeight[pomic$babyWeight > 60] <- NA
# 身高无异常值
boxplot(pomic$babyLength)

# 计算 z 评分
pomic <- addWGSR(data = pomic, sex = "sex", firstPart = "babyWeight",
                 secondPart = "ageBaby", index = "wfa",output = "waz")

pomic <- addWGSR(data = pomic, sex = "sex", firstPart = "babyLength",
                 secondPart = "ageBaby", index = "hfa",output = "haz")

# 筛选样本
pomicSel <- pomic %>%
  # 选取月龄不小于3且非早产
  filter(ageBaby > 60 & preterm == 1)

# 只排除早产儿
pomic.p <- pomic %>%
  filter(preterm==1)


save(pomic.p,file="data/selected proteomic samples without preterm.Rdata")
save(pomicSel,file = "data/selected proteomic samples.Rdata")

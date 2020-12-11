library(tidyverse)
library(zscorer)
# 导入之前整理好的脂质组学数据
# 仍然没有哈尔滨和郑州的数据
# 因为这两个城市的样本id与检测数据没法一一对应
pomic_all <- read.csv("data/ldomic_all.csv",stringsAsFactors = F)
pomic_all <- ldomic_all %>%
  rename(id=sample)


# 去掉NA超过20%的脂质数据

ld <- ldomic_all[c(1,2)]

for (i in 3:160) {
  
  if (sum(is.na(ldomic_all[,i]))<32){
    ld <- cbind(ld,ldomic_all[i])
  }
  
}

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
ldomic <- inner_join(ld,quesSel,by="id")
ldSel <- ldomic %>%
  filter(ageBaby > 40 & preterm == 1)

# 查看婴儿体重和身长数据是否有问题
## 体重
hist(ldSel$babyWeight)
summary(ldSel$babyWeight)
boxplot(ldSel$babyWeight)# 有一个大于60的数据点，去掉
ldSel$babyWeight[ldSel$babyWeight > 60] <- NA

## 身高
hist(ldSel$babyLength)
summary(ldSel$babyLength)
boxplot(ldSel$babyLength)



# 计算 z 评分
ldSel <- addWGSR(data = ldSel, sex = "sex", firstPart = "babyWeight",
                  secondPart = "ageBaby", index = "wfa",output = "waz")

ldSel <- addWGSR(data = ldSel, sex = "sex", firstPart = "babyLength",
                  secondPart = "ageBaby", index = "hfa",output = "haz")

boxplot(ldSel$waz)
boxplot(ldSel$haz)

save(ldSel,file = "data/selected lipidomic samples.Rdata")

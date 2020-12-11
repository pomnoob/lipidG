load(file = "data/selected lipidomic samples.Rdata")
library(tidyverse)
##############################################################################
##############################################################################

# haz分析策略
# 1. 以 haz 绝对值为标准，大于1为发育正常，小于-1为发育迟缓，n分别为44和24
# 2. 以 haz 的三分位为标准，T1和T3分别有42个样本

# 策略1
ldSel <- ldSel %>%
  mutate(nhaz=case_when(haz > 1~1,
                        haz < -1~2))
table(ldSel$nhaz)
# 选出发育迟缓和正常的数据
haz1 <- ldSel %>%
  filter(nhaz == 1 | nhaz == 2)

#t 检验
mean1_1 <- double()
mean1_2 <- double()
pvalue1 <- double()
lipid1 <- character()
for (i in 3:139) {
  # 做t检验
  t.result1 <- t.test(haz1[[i]]~haz1$nhaz,
                     na.action = na.omit)
  pvalue1[length(pvalue1)+1] <- t.result1[["p.value"]]
}

pvalue1
p.adjust(pvalue1,method ="fdr" )

##############################################################################

# 策略2
ldSel$thaz <- ntile(ldSel$haz,3)
table(ldSel$thaz)
# 选出T1和T2的数据
haz2 <- ldSel %>%
  filter(thaz == 1 | thaz == 3)
#t 检验
mean2_1 <- double()
mean2_3 <- double()
pvalue2 <- double()
lipid2 <- character()
for (i in 3:139) {
  # 做t检验
  t.result2 <- t.test(haz2[[i]]~haz2$thaz,
                      na.action = na.omit)
  pvalue2[length(pvalue2)+1] <- t.result2[["p.value"]]
}

pvalue2
p.adjust(pvalue2,method ="hommel" )

save(haz1,file="data/infant length for age strategy 1.Rdata")
save(haz2,file="data/infant length for age strategy 2.Rdata")

##############################################################################
##############################################################################

# 体重以三分位为标准进行分析
ldSel$twaz <- ntile(ldSel$waz,3)
table(ldSel$twaz)
# 选出T1和T2的数据
waz <- ldSel %>%
  filter(twaz == 1 | twaz == 3)

#t 检验
mean3_1 <- double()
mean3_3 <- double()
pvalue3 <- double()
lipid3 <- character()
for (i in 3:139) {
  # 做t检验
  t.result3 <- t.test(waz[[i]]~waz$twaz,
                      na.action = na.omit)
  pvalue3[length(pvalue3)+1] <- t.result3[["p.value"]]
}

pvalue3
p.adjust(pvalue3,method ="fdr" )

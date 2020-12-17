# 选择所有非早产儿样本

load(file = "data/selected proteomic samples without preterm.Rdata")
library(tidyverse)
##############################################################################
##############################################################################

# haz分析策略
# 1. 以 haz 绝对值为标准，大于1为发育正常，小于-1为发育迟缓，n分别为44和24
# 2. 以 haz 的三分位为标准，T1和T3分别有42个样本

# 策略1
pomic.p <- pomic.p %>%
  mutate(nhaz=case_when(haz > 1~1,
                        haz < -1~2))
table(pomic.p$nhaz)
# 选出发育迟缓和正常的数据
haz1 <- pomic.p %>%
  filter(nhaz == 1 | nhaz == 2)
# 检查各个蛋白质分布情况
write.csv(haz1,file = "data/HAZ proteome without preterm 1 and -1.csv",
          row.names = F)
#t 检验
mean1_1 <- double()
mean1_2 <- double()
pvalue1 <- double()
lipid1 <- character()
for (i in 3:195) {
  # 做t检验
  t.result1 <- t.test(haz1[[i]]~haz1$nhaz,
                     na.action = na.omit)
  pvalue1[length(pvalue1)+1] <- t.result1[["p.value"]]
}

pvalue1
p1.adj <- p.adjust(pvalue1,method ="fdr" )
p1.adj[p1.adj<0.05]
##############################################################################

# 策略2
pomic.p$thaz <- ntile(pomic.p$haz,3)
table(pomic.p$thaz)
# 选出T1和T2的数据
haz2 <- pomic.p %>%
  filter(thaz == 1 | thaz == 3)
#t 检验
mean2_1 <- double()
mean2_3 <- double()
pvalue2 <- double()
lipid2 <- character()
for (i in 3:195) {
  # 做t检验
  t.result2 <- t.test(haz2[[i]]~haz2$thaz,
                      na.action = na.omit)
  pvalue2[length(pvalue2)+1] <- t.result2[["p.value"]]
}

pvalue2
p2.adj <- p.adjust(pvalue2,method ="fdr" )
p2.adj[p2.adj<0.05]
# 导出数据
haz2.exp <- haz2 %>%
  select(id, thaz, 3:195) %>%
  rename(sample = id,
         group = thaz)
write.csv(haz2.exp,file = "data/HAZ proteome without preterm T1 and T3.csv",
          row.names = F)
##############################################################################
##############################################################################

# 体重以三分位为标准进行分析
pomic.p$twaz <- ntile(pomic.p$waz,3)
table(pomic.p$twaz)

# 选出T1和T2的数据
waz.p <- pomic.p %>%
  filter(twaz == 1 | twaz == 3)
# 保存数据
save(waz.p,file = "data/proteomic-infant weight for age without preterm.Rdata")
wazp.exp <- waz.p %>%
  select(id,twaz,3:195) %>%
  rename(group=twaz)
write.csv(wazp.exp,file = "data/WAZ proteome without preterm.csv",
          row.names = F)
#t 检验
mean3_1 <- double()
mean3_3 <- double()
pvalue3 <- double()
prot3 <- character()
colnameWAZ <- colnames(waz.p)
for (i in 3:195) {
  # 做t检验
  t.result3 <- t.test(waz.p[[i]]~waz.p$twaz,
                      na.action = na.omit)
  pvalue3[length(pvalue3)+1] <- t.result3[["p.value"]]
  mean3_1[length(mean3_1)+1] <- t.result3[["estimate"]][["mean in group 1"]]
  mean3_3[length(mean3_3)+1] <- t.result3[["estimate"]][["mean in group 3"]]
  prot3[length(prot3)+1] <- colnameWAZ[i]
}

p3.adj <- p.adjust(pvalue3,method ="fdr" )
p3.adj[p3.adj<0.05]

# 筛选出差异显著的蛋白质
protWAZ.p <- data.frame(pid=prot3,
                      mean_T1=mean3_1,
                      mean_T3=mean3_3,
                      pval=pvalue3,
                      padj=p3.adj)

# FDR 以后差异显著
protWAZpre.adj <- protWAZ.p %>%
  filter(padj<0.05) #n=90
# 原始 p 值
protWAZpre.p <- protWAZ.p %>%
  filter(pval<0.05) #n=103

# 保存数据
save(protWAZ.p,file = "data/蛋白组 WAZ T1和T3两组之间t检验-非早产.Rdata")
save(protWAZpre.adj,file = "data/FDR校正后差异显著的蛋白质 WAZ T1和T3两组之间t检验-非早产.Rdata")
save(protWAZpre.p,file = "data/不做校正差异显著的蛋白质 WAZ T1和T3两组之间t检验-非早产.Rdata")



# 以waz大于1或小于-1为标准
pomicSel <- pomicSel %>%
  mutate(nwaz=case_when(waz > 1~1,
                        waz < -1~2))
table(pomicSel$nwaz)

# 选出1和2的数据
nwaz <- pomicSel %>%
  filter(nwaz == 1 | nwaz == 2)
# 保存数据
save(nwaz,file = "data/proteomic-infant weight for age 1 and -1.Rdata")
#t 检验
mean4_1 <- double()
mean4_3 <- double()
pvalue4 <- double()
prot4 <- character()
colnameWAZ <- colnames(nwaz)
for (i in 3:195) {
  # 做t检验
  t.result4 <- t.test(nwaz[[i]]~nwaz$nwaz,
                      na.action = na.omit)
  pvalue4[length(pvalue4)+1] <- t.result4[["p.value"]]
  mean4_1[length(mean4_1)+1] <- t.result4[["estimate"]][["mean in group 1"]]
  mean4_3[length(mean4_3)+1] <- t.result4[["estimate"]][["mean in group 2"]]
  prot4[length(prot4)+1] <- colnameWAZ[i]
}

p4.adj <- p.adjust(pvalue4,method ="fdr" )
p4.adj[p4.adj<0.05]


# 筛选出差异显著的蛋白质
protWAZ.n <- data.frame(pid=prot4,
                      mean1=mean4_1,
                      mean2=mean4_3,
                      pval=pvalue4,
                      padj=p4.adj)
# FDR 以后差异显著
protWAZ.adj.n <- protWAZ.n %>%
  filter(padj<0.05) #n=18
# 原始 p 值
protWAZ.pn <- protWAZ.n %>%
  filter(pval<0.05) #n=72

# 以体重大于1或小于-1为标准
pomic.p <- pomic.p %>%
  mutate(nwaz=case_when(waz > 1~1,
                        waz < -1~2))
table(pomic.p$nwaz)
# 选出发育迟缓和正常的数据
waz2.p <- pomic.p %>%
  filter(nwaz == 1 | nwaz == 2)
#t 检验
mean5_1 <- double()
mean5_3 <- double()
pvalue5 <- double()
prot5 <- character()
colnameWAZ <- colnames(waz2.p)
for (i in 3:195) {
  # 做t检验
  t.result5 <- t.test(waz2.p[[i]]~waz2.p$nwaz,
                      na.action = na.omit)
  pvalue5[length(pvalue5)+1] <- t.result5[["p.value"]]
  mean5_1[length(mean5_1)+1] <- t.result5[["estimate"]][["mean in group 1"]]
  mean5_3[length(mean5_3)+1] <- t.result5[["estimate"]][["mean in group 2"]]
  prot5[length(prot5)+1] <- colnameWAZ[i]
}

p5.adj <- p.adjust(pvalue5,method ="fdr" )
p5.adj[p5.adj<0.05]

write.csv(waz2.p,file = "data/WAZ proteome without preterm 1 and -1.csv",
          row.names = F)
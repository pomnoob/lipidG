# 这一个版本的参数：
# 所有的蛋白质，蛋白质的筛选在MetaboAnalyst上面完成，参数是：50%，min，RSD
# 数据标准化和归一化参数：median，Log，Pareto
# 所有的参数可以在该文件夹的Rhistory中找到
library(tidyverse)
# 导入PLS，t-test，FC结果
plsVIP.all <- read.csv(file = "redo with parameter MLP/plsda_vip.csv",stringsAsFactors = F)
fc.all <- read.csv(file = "redo with parameter MLP/fold_change.csv",stringsAsFactors = F)
uni.all <- read.csv(file = "redo with parameter MLP/wilcox_rank.csv",stringsAsFactors = F)

# 数据修整
plsVIP.all <- plsVIP.all %>%
  select(X,Comp..1) %>%
  rename(feature=X,vip = Comp..1)

fc.all <- fc.all %>%
  rename(feature=X
    ,fc = Fold.Change,
         logfc=log2.FC.)

uni.all <- uni.all %>%
  rename(feature=X,
         p=p.value,
         logp=X.log10.p.) %>%
  select(-V)

result <- left_join(plsVIP.all,fc.all,by="feature")

result <- left_join(result,uni.all,by = "feature")


# 导入之前做好的蛋白质ID转Gene ID的文件和代码
pro.id <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)

# refine the names
pro.id$gname <- toupper(pro.id$gname)
pro.id$pname <- str_replace_all(pro.id$pname,
                                "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                "")
pro.id$pname <- str_to_title(pro.id$pname)

result.gene <- left_join(result,pro.id,by="feature")
result.sig <- result.gene %>%
  filter(vip>1 & (fc >= 1.5 | fc< 0.67))
# 导出CSV,用 Crytoscape进行网络分析
write.csv(pSel.gene,file = "data/all protein with vip and fc for network.csv",
          row.names = F)
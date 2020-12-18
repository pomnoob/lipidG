# 对所有蛋白都进行纳入，然后根据筛选的样本进行feature selection
# 在MetaboAnalyst online版本上做分析，下载VIP和FC的结果
plsVIP.all <- read.csv(file = "data/plsda_vip all protein.csv",stringsAsFactors = F)
fc.all <- read.csv(file = "data/fold_change all protein.csv",stringsAsFactors = F)

# 数据修整
plsVIP.all <- plsVIP.all %>%
  select(feature,Comp..1) %>%
  rename(vip = Comp..1)

fc.all <- fc.all %>%
  rename(fc = Fold.Change,
         logfc=log2.FC.)

pSel <- left_join(plsVIP.all,fc.all,by = "feature")


# 导入之前做好的蛋白质ID转Gene ID的文件和代码
pro.id <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)

# refine the names
pro.id$gname <- toupper(pro.id$gname)
pro.id$pname <- str_replace_all(pro.id$pname,
                                "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                "")
pro.id$pname <- str_to_title(pro.id$pname)

pSel.gene <- left_join(pSel,pro.id,by="feature")

# 导出CSV,用 Crytoscape进行网络分析
write.csv(pSel.gene,file = "data/all protein with vip and fc for network.csv",
          row.names = F)

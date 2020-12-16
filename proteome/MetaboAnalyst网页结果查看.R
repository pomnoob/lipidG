# 导入从MetaboAnalyst网页上下载的结果
# 对比后发现，t test的结果包含pls的结果，可以先看t test的结果
proSig <- read.csv(file = "MetaboAnalystOnLine/volcano.csv",stringsAsFactors = F)
proSigName <- proSig %>%
  select(feature)
# 导入蛋白质ID和Gene ID文件
pro.id.n <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)
pro.id.n$gname <- toupper(pro.id.n$gname)
pro.id.n$pname <- str_replace_all(pro.id.n$pname,
                                  "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                  "")
pro.id.n$pname <- str_to_title(pro.id.n$pname)
proSigName <- left_join(proSig,pro.id.n,by= "feature")

# 保存数据
write.csv(proSigName,file = "data/显著差异的蛋白质id.csv")

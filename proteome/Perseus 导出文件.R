# 导入之前做好的蛋白质ID转Gene ID的文件和代码
pro.id <- read.csv(file = "data/protein id to name.csv",stringsAsFactors = F)

# refine the names
pro.id$gname <- toupper(pro.id$gname)
pro.id$pname <- str_replace_all(pro.id$pname,
                                  "\\([[:upper:]][[:upper:][:digit:]]*\\)",
                                  "")
pro.id$pname <- str_to_title(pro.id$pname)

# 导入WAZ proteome 数据
wazPro <- read.csv(file = "data/WAZ proteome without preterm.csv")

# 转置
wazPro.t <- as.data.frame(t(wazPro))
wazPro.t <- wazPro.t[-1,]
names(wazPro.t) <- wazPro$id
# pid 作为一列
wazPro.t <- rownames_to_column(wazPro.t,var = "feature")

wazProAll <- left_join(wazPro.t,pro.id,by="feature")

write.csv(wazProAll,file = "data/Perseus usable data.csv",
          row.names = F)

write.table(wazProAll,file = "data/Perseus usable data.txt",
            sep = "\t",
            row.names = F)

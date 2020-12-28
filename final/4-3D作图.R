# 3D 作图
library(tidyverse)
library(rgl)
# 导入PLS score数据
pls.score <- read.csv(file = "final/Metab_plsda_score.csv",
                      stringsAsFactors = F)

pls.score <- pls.score %>%
  select(sample=X,c1=Comp.1,
         c2=Comp.2,
         c3=Comp.3)

# 加入分组信息
prot.orig <- read.csv(file = "final/Metab_data_normalized_t.csv",
                      stringsAsFactors = F)

prot.orig <- prot.orig %>%
  select(sample,group)

# 加入分组信息

pls.score <- inner_join(pls.score,prot.orig,by="sample")

pls.score.two <- pls.score %>%
  filter(group == 1| group == 3)

with(pls.score.two, plot3d(c1, c2, c3, 
                  type="s", col=as.numeric(group),radius = 0.1,
                  xlab = "Comp 1",
                  ylab = "Comp 2",
                  zlab = "Comp 3"))
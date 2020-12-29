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

#' Get colors for the different levels of 
#' a factor variable
#' 
#' @param groups a factor variable containing the groups
#'  of observations
#' @param colors a vector containing the names of 
#   the default colors to be used
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

cols.3d <- get_colors(pls.score.two$group,c("#1965B0","#DC050C"))

with(pls.score.two, plot3d(c1, c2, c3, 
                  type="s", col=cols.3d,radius = 0.1,box=F,
                  xlab = "Component 1",
                  ylab = "Component 2",
                  zlab = "Component 3"))

rgl.postscript("final/3d plot.svg",fmt="svg")
rgl.snapshot(filename = "final/3d plot.png")
# add legend

legend3d("right", legend = paste('Tertile', c("1","3")), col=c("#1965B0","#DC050C"),pch = 16,cex=0.9, inset=c(0.03))



proNorm <- read.csv(file = "Download/data_normalized_id.csv",stringsAsFactors = F)
proOrig <- read.csv(file = "Download/data_original_id.csv",stringsAsFactors = F)

ggplot(data=proOrig)+geom_boxplot(aes(x=factor(Label),y=P25311))
ggplot(data=proNorm)+geom_boxplot(aes(x=factor(Label),y=P25311))

library(rgl)
library(tidyverse)
x <- seq(0,2048,by=1)
y <- (80+10*x)/(80+x)
z <- (1024/7+16*x)/(1024/7+x)
data <- data.frame(x,y,z)
ggplot(data)+geom_line(aes(x,y),color="black")+
  geom_line(aes(x,z),color="red")
write.csv(data,file = "商品.csv")

library(tidyverse)
x <- seq(0,2048,by=1)

y <- (80+10*x)/(80+x)
z <- (1024/7+16*x)/(1024/7+x)
data <- data.frame(x,y,z)
ggplot(data)+geom_line(aes(x,y),color="black")+
  geom_line(aes(x,z),color="red")
write.csv(data,file = "商品.csv")

# 两段曲线，商品数量为128时，倍数=8，商品数量为10时倍数为2.5，倍数最多为16
x <- seq(0,2048,by=1)

data <- data.frame(x)
data <- data %>%
  mutate(y = case_when(x<=10~1+3/20*x,
                       x>10~(4720/11+16*(x-10))/(1888/11+(x-10))
                       )
         )
ggplot(data)+geom_line(aes(x,y),color="black")
write.csv(data,file = "商品3点.csv")
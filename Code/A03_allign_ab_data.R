library(data.table)
library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)

##### load data ####
load('Result/A01_compiled_data.rda')

##### cluster ab1 #####
ab1_wide = ab1%>%spread(key = day,value = antibody_value)
t1 = apply(ab1_wide, 1, function(x){sum(is.na(x))})
ab1_wide = ab1_wide[t1<=2,]

D = 1-cor(t(ab1_wide[,c(2,3,4,5,6)]),use = 'pairwise.complete.obs')
HC = hclust(as.dist(D),method = 'ward.D2')
CT = cutree(HC, k = 2)

ab1_wide = ab1_wide%>%
  mutate(ab0 = (`0`<0.1))%>%
  mutate(ab0 = factor(ab0,labels = c('Ab0+','Ab0-')))
ab1_long = ab1_wide%>%gather(key='day',value = 'ab',-id,-ab0)%>%
  mutate(day = as.integer(day))

pdf('Result/A03_ab_trajectory.pdf',width=6,height=5)
p = ggplot(ab1_long,aes(x=(day),y = ab,group=id))+
  geom_line()+theme_bw()
plot(p)

p = ggplot(ab1_long,aes(x=(day),y = ab,group=id,color=factor(ab0)))+
  geom_line(size = 0.2)+theme_bw()
plot(p)
dev.off()

cl1 = ab1_wide%>%dplyr::select(id,ab0)

##### correlation with demo and clinical data ######
cl1 = cl1 %>%
  inner_join(timeTo1,by='id')%>%
  inner_join(demo1%>%dplyr::select(-tx),by='id')

table(cl1$ab0,cl1$tx)
table(cl1$ab0,cl1$sex)
boxplot(age~ab0,cl1)
pdf('Result/A03_Ab0_vs_clin.pdf',5,5)
boxplot(time2prime~ab0,cl1)
boxplot(time2prog~ab0,cl1)
boxplot(time2res~ab0,cl1)
dev.off()

###### save results #####
write.csv(cl1[,1:2],'Result/A03_ab_cluster.csv',row.names = F)

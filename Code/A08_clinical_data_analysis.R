library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(survival)
library(fastcluster)
library(cluster)
library(survminer)
library(umap)
library(Rtsne)
library(enrichR)
library(DESeq2)
library(ggrepel)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
neut_ab = read.csv('Data/Antibodies/prelim %Neut at 1-50 for Pras.csv')
adj_day = read.csv('Result/A02_adj_day.csv')
base_symptom = read.csv('Data/ClinicalData/lambda_baseline_symptoms.csv')

##### get duration of symptom at day 0 #####
symp_start = base_symptom%>%
  group_by(participantId)%>%
  summarise(symp_start_day=0-max(duration,na.rm = T))%>%
  dplyr::rename(id = participantId)
symp_start$symp_start_day[!is.finite(symp_start$symp_start_day)]=NA

##### compile clinical information #####
info1=adj_day%>%group_by(id)%>%summarise(adj_day=mean(adj_day)-2.5) 
info1 = info1%>%mutate(id = as.integer(id))%>%
  inner_join(demo1,by='id')%>%
  inner_join(timeTo1[,-2],by='id')

ab1_wide = ab1%>%spread(key = day,value = antibody_value)
colnames(ab1_wide)[-1]=paste0('ab',colnames(ab1_wide)[-1])
info1 = info1%>%left_join(ab1_wide,by='id')

info1 = left_join(info1,symp1,by = 'id')
info1 = full_join(info1,symp_start,by='id')

colnames(neut_ab)=c('id','neut28','neutM4')
info1 = info1%>%full_join(neut_ab,by='id')

info1 = info1%>%
  dplyr::select(adj_day,sex,age,bmi,tx,time2prime,time2prog,time2res,
         ab0,ab5,ab14,ab28,ab120,severitycat,neut28,neutM4,symp_start_day)

info1 = info1%>%
  mutate(severitycat=factor(severitycat,
                            levels =c('Asymptomatic',
                                        "Moderate Symptomatic",
                                        "Severe")))%>%
  mutate(severitycat = as.integer(severitycat))%>%
  mutate(tx = (tx=='Lambda')*1)

##### get correlation #####
D1 = 1-(cor(info1,use = 'pairwise.complete.obs',method = 'spearman'))
pca1 = cmdscale(as.dist(D1),eig=TRUE, k=2)$points[,1:2]
colnames(pca1) = c('PC1','PC2')
pca1= cbind.data.frame(gene = colnames(info1),pca1[,1:2])
colnames(pca1)=c('gene','u1','u2')

pair1 = combn(colnames(info1),2)%>%t()%>%as.data.frame()
pair1$cor = NA
pair1$p = NA
for (i in 1:nrow(pair1)) {
  t1 = cor.test(unlist(info1[,pair1$V1[i]]),
                unlist(info1[,pair1$V2[i]]),
                use = 'pairwise.complete.obs',method = 'spearman')
  pair1$cor[i]= t1$estimate
  pair1$p[i]= t1$p.value
}
pair1 = pair1%>%filter(p<0.05)#%>%mutate(cor = cor>0)%>%
  #mutate(cor = factor(cor,labels = c('sig-','sig+')))%>%
  #mutate(group = paste0(V1,V2))%>%
  #select(-p)%>%gather(key='key',value = 'gene',V1:V2)%>%
  #inner_join(pca1,by='gene')

p = ggplot(pca1,aes(x=u1,y = u2))+
  geom_point()+theme_bw()+
  geom_label_repel(aes(label = gene))#+
  #geom_line(data=pair1,aes(x=u1,y=u2,group=group,color=cor))
plot(p)

##### plot bar plot #####
pair1 = pair1%>%
  filter(!(grepl('^ab',V1)&grepl('^ab',V2)))%>%
  filter(!(grepl('^neu',V1)&grepl('^neu',V2)))%>%
  mutate(pairs = paste(V1,V2, sep = ' vs. '))

p = ggplot(pair1,aes(x=reorder(pairs,cor), y = cor, fill = cor>0))+
  geom_bar(stat = 'identity')+coord_flip()+theme_bw()
plot(p)

p = ggplot(info1%>%filter(!is.na(sex)),
           aes(x = factor(sex), y = adj_day,color =factor(sex) ))+
  geom_boxplot()+theme_bw()
plot(p)

p = ggplot(info1%>%filter(!is.na(sex)),
           aes(x = factor(severitycat), 
               y = neutM4,color =factor(severitycat) ))+
  geom_boxplot()+theme_bw()
plot(p)

##### plot adj data #####
p = ggplot(info1%>%filter(symp_start_day>(-10)),
           aes(x=adj_day,y = symp_start_day))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)

p = ggplot(info1,aes(x=adj_day,y = time2prime))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)

plot_df = adj_day%>%group_by(id)%>%
  summarise(adj_day=mean(adj_day)-2.5)%>%
  arrange(adj_day)
p = ggplot(plot_df,
           aes(x=reorder(factor(id),adj_day),y=adj_day))+
  geom_bar(stat = 'identity')+theme_bw()+coord_flip()
plot(p)
 

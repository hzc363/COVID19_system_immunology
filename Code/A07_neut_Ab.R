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

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
neut_ab = read.csv('Data/Antibodies/prelim %Neut at 1-50 for Pras.csv')

##### correlation between IgG and neut_ab #####
neut_ab = gather(neut_ab,key='day',value = "neut_value",-id)
neut_ab = neut_ab%>%
  mutate(day=gsub('D','',day))%>%
  mutate(day = as.integer(day))

df_combine = inner_join(neut_ab,ab1,by=c('day','id'))
df_cor = df_combine%>%group_by(day)%>%
  summarise(cor=cor(neut_value,antibody_value,
                    use = 'complete',method = 'spearman'))
p = ggplot(df_combine,aes(x=antibody_value,y = neut_value))+
  geom_point()+
  facet_wrap(~factor(day),scales = 'free')+theme_bw()
plot(p)

##### symptom vs neut ab ####
df_combine = symp1%>%dplyr::select(id,severitycat)%>%
  inner_join(df_combine,by='id')

p = ggplot(df_combine,aes(y=antibody_value,x = factor(day),
                          color=severitycat))+
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))+
  theme_bw()
plot(p)

p = ggplot(df_combine,aes(y=neut_value,x = factor(day),
                          color=severitycat))+
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))+
  theme_bw()
plot(p)

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
library(lme4)
library(lmerTest)
library(emmeans)
library(gplots)
library(stringr)
#library(xCell)
library(glmnet)
library(randomForest)
##### get all results #####
severity1 = read.csv('Result/A12_severity_results.csv')
clear1 = read.csv('Result/A13_clearance_results_auc.csv')
tcell1 = read.csv('Result/A14_t_cell_results.csv')
Ab1 = read.csv('Result/A15_Ab_M4_results_new.csv')
Ab2 = read.csv('Result/A15_Ab_M7_results_new.csv')

tAll = list(severity1,clear1,tcell1,Ab1,Ab2)
f1 = function(x,y){inner_join(x,y,by='gene')}
tAll = Reduce(f1,tAll)
tAll = tAll%>%
  dplyr::select(variable = gene,
                severity_t=t.x,severity_p=p.x,
                viral_t=t.y,viral_p=p.y,
                tcell_t = t.x.x,tcell_p=p.x.x,
                Ab_28days_t = t.y.y,Ab_28days_p=p.y.y,
                Ab2_7M_t = t, Ab2_7M_p = p)%>%
  filter(grepl('olink_|GO_',variable))%>%
  filter(nchar(variable)<40)

write.csv(tAll,'Result/A25_suplement_table.csv',row.names = F)
##### bar plot for pathways viral#####
tAll=tAll%>%filter(variable!='GO_immune_response')
dfPlot = tAll%>%
  filter(grepl("GO_",variable))%>%
  mutate(variable = gsub('GO_','',variable))%>%
  mutate(variable = gsub('_',' ',variable))%>%
  group_by(variable)%>%
  summarise_all(function(x){x[1]})%>%
  top_n(10,-viral_p)

p = ggplot(dfPlot,aes(x=viral_t,y=reorder(variable,viral_t)))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A25_viral_bar_auc_GO.pdf',width =3.5,height = 3)
plot(p)
dev.off()
##### bar plot for pathways tcell#####
dfPlot = tAll%>%
  filter(grepl("GO_",variable))%>%
  mutate(variable = gsub('GO_','',variable))%>%
  mutate(variable = gsub('_',' ',variable))%>%
  group_by(variable)%>%
  summarise_all(function(x){x[1]})%>%
  top_n(10,-tcell_p)

p = ggplot(dfPlot,aes(x=tcell_t,y=reorder(variable,tcell_t)))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A25_tcell_bar_auc_GO.pdf',width = 3.5,height = 3)
plot(p)
dev.off()

##### bar plot for pathways day28 ab#####
dfPlot = tAll%>%
  filter(grepl("GO_",variable))%>%
  mutate(variable = gsub('GO_','',variable))%>%
  mutate(variable = gsub('_',' ',variable))%>%
  group_by(variable)%>%
  summarise_all(function(x){x[1]})%>%
  filter(!grepl('thym|selection',variable))%>%
  top_n(10,-Ab_28days_p)

p = ggplot(dfPlot,aes(x=Ab_28days_t,y=reorder(variable,Ab_28days_t)))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A25_day28_Ab_bar_auc_GO.pdf',width = 3.5,height = 3)
plot(p)
dev.off()

##### bar plot for pathways M7 ab#####
dfPlot = tAll%>%
  filter(grepl("GO_",variable))%>%
  mutate(variable = gsub('GO_','',variable))%>%
  mutate(variable = gsub('_',' ',variable))%>%
  group_by(variable)%>%
  summarise_all(function(x){x[1]})%>%
  filter(!grepl('glia',variable))%>%
  top_n(10,-Ab2_7M_p)

p = ggplot(dfPlot,aes(x=Ab2_7M_t,y=reorder(variable,Ab2_7M_t)))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A25_7M_Ab_bar_auc_GO.pdf',width = 3.5,height = 3)
plot(p)
dev.off()

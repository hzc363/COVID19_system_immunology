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
tcell_result = read.csv('Result/A05_tcell_result.csv')
olink_result = read.csv('Result/A04_olink_result.csv')
olink = read.csv('Result/A04_olink_with_symbol.csv')


##### prepare data #####
olink_result =olink_result%>%
  filter(p_M4<0.05|p_symp_adj<0.05)
olink = olink[,c("id",'day',olink_result$gene)]

olink = gather(olink,key = 'protein',value = 'olink_value',-id,-day)
olink=olink%>%mutate(protein=gsub('_.*','',protein))%>%
  group_by(id,day,protein)%>%
  summarise(olink_value=mean(olink_value))%>%
  as.data.frame()

tcell_result = tcell_result%>%
  filter(p_M4<0.05|p_symp_adj<0.05)
tcell = tcell[,c("id",tcell_result$gene)]
tcell = gather(tcell,key='cell',value = 'tcell_value',-id)%>%
  mutate(id = as.integer(id))
OT1 = inner_join(olink,tcell,by='id')
OT1 = OT1 %>%mutate(combo = paste(protein,cell))

##### get correlation #####
OT_result = data.frame()
for(c1 in unique(OT1$combo)){
  t1 = OT1%>%filter(combo==c1)
  LM = lm(tcell_value~day+olink_value,t1)
  LM = summary(LM)$coefficients
  t1 = t1[1,c('protein','cell')]%>%
    mutate(p = LM['olink_value','Pr(>|t|)'])%>%
    mutate(t = LM['olink_value','t value'])%>%
    mutate(b = LM['olink_value','Estimate'])
  OT_result = rbind(OT_result,t1)
}

OT_result = OT_result%>%mutate(p_adj = p.adjust(p,'fdr'))%>%
  mutate(cell = gsub('_bgsub','',cell))

##### plot heatmap #####
pdf('Result/A06_heatmap.pdf',height = 15,width = 10)
t1 = OT_result%>%filter(p_adj<0.01)
plot_df = OT_result%>%filter(protein%in%t1$protein & cell%in%t1$cell)
plot_df =plot_df%>%dplyr::select(protein,cell,t)%>%spread(key = cell,value = t)
plot_df = data.frame(plot_df[,-1],row.names = plot_df[,1])
r1 = hclust(dist(plot_df),method = 'ward.D')$order
c1 = hclust(dist(t(plot_df)),method = 'ward.D')$order

pheatmap(plot_df[r1,c1],cluster_rows=F, cluster_cols = F)

t1 = OT_result%>%filter(p_adj<0.01)
plot_df = OT_result%>%filter(protein%in%t1$protein & cell%in%t1$cell)
plot_df =plot_df%>%dplyr::select(protein,cell,p_adj)%>%
  spread(key = cell,value = p_adj)
plot_df = data.frame((plot_df[,-1]<0.05)*1,row.names = plot_df[,1])
pheatmap(plot_df[r1,c1],cluster_rows=F, cluster_cols = F)
dev.off()

##### plot heatmap subset #####
pdf('Result/A06_heatmap_sub.pdf',height = 7,width = 3)
t1 = OT_result%>%
  filter(p_adj<0.05&grepl('TNF$|INFy$|IL21$|IL10$|CD107a$',cell))
plot_df = OT_result%>%filter(protein%in%t1$protein & cell%in%t1$cell)
plot_df =plot_df%>%dplyr::select(protein,cell,t)%>%spread(key = cell,value = t)
plot_df = data.frame(plot_df[,-1],row.names = plot_df[,1])
r1 = hclust(dist(plot_df),method = 'ward.D')$order
c1 = hclust(dist(t(plot_df)),method = 'ward.D')$order

pheatmap(plot_df[r1,c1],cluster_rows=F, cluster_cols = F)

t1 = OT_result%>%
  filter(p_adj<0.05&grepl('TNF$|INFy$|IL21$|IL10$|CD107a$',cell))
plot_df = OT_result%>%filter(protein%in%t1$protein & cell%in%t1$cell)
plot_df =plot_df%>%dplyr::select(protein,cell,p_adj)%>%
  spread(key = cell,value = p_adj)
plot_df = data.frame((plot_df[,-1]<0.05)*1,row.names = plot_df[,1])
pheatmap(plot_df[r1,c1],cluster_rows=F, cluster_cols = F)
dev.off()

##### plot top correlation #####
plot_df = OT1 %>%
  filter(protein=='CCL7',cell=='SS1_bgsub_CD4p_TNF')
p = ggplot(plot_df,aes(y=tcell_value,x=olink_value,color = factor(day)))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()+
  xlab('KRT19')+ylab('SS1_bgsub_CD4p_TNF')
plot(p)

plot_df = OT1 %>%
  filter(protein=='IL6',cell=='SS1_bgsub_CD4p_TNF')
p = ggplot(plot_df,aes(y=tcell_value,x=olink_value,color = factor(day)))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()+
  xlab('IL6')+ylab('SS1_bgsub_CD4p_TNF')
plot(p)



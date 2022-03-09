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
library(corrplot)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
#t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
#       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
geneDf = read.csv('A27_top_go_genes.csv')
t1 = paste0('Gene_',geneDf$gene_symbol)

RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))|colnames(RNAseq)%in%t1]
RNAseq=RNAseq%>%mutate(Info_id = as.character(Info_id))

colnames(olink)[-(1:2)]=paste0('olink_',colnames(olink)[-(1:2)])
RNAseq = left_join(RNAseq,olink,by=c('Info_id'='id','Info_day'='day'))

colnames(ab1)[-(1:2)]=paste0('IgG_',colnames(ab1)[-(1:2)])
ab1 = ab1 %>%mutate(id = as.character(id))
RNAseq = left_join(RNAseq,ab1,by=c('Info_id'='id','Info_day'='day'))

colnames(lab1)[-(1:2)]=paste0('Lab_',colnames(lab1)[-(1:2)])
lab1 = lab1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,lab1,by=c('Info_id'='id','Info_day'='day'))

colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### correlate with t cell #####
tcell = read.csv('Data/D5 D28_booleans_AIM.csv',
                         check.names = F)%>%
  dplyr::select(id, day=days, `Tfh/ICOS+`,Stim)%>%
  filter(day==28)%>%na.omit()%>%mutate(id = as.character(id))

df1 = tcell%>%
  dplyr::select(Info_id=id,Stim,
                tcell=`Tfh/ICOS+`)%>%
  inner_join(RNAseq,by='Info_id')%>%
  mutate(Info_onset = Info_day-Info_onset)

df1$tcell[df1$tcell<0]=0

# time vs t cell
dfPlot = df1%>%select(Info_onset,tcell,Info_id)%>%
  arrange(Info_onset)%>%
  group_by(Info_id)%>%
  summarise_all(function(x){x[1]})%>%na.omit()
p = ggplot(dfPlot,aes(x=Info_onset+28,y=tcell))+
  geom_point()+geom_smooth(method = 'lm')+
  ylab('Tfh response')+
  xlab('Time since symptom onset')+
  theme_bw()
pdf("Result/A29_onset_vs_Tfh.pdf",5,5)
plot(p)
dev.off()
cor.test(dfPlot$Info_onset,dfPlot$tcell)



dfResult = data.frame()
for (i in grep('GO_|^olink_|Gene_|EVE_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c("tcell","Info_day",'Stim',
              colnames(df1)[i])]%>%na.omit()
  colnames(t1)[4]='y'
  LM1 = lm(tcell~Info_day+y+Stim,t1)
  LM1 = summary(LM1)
  
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['y','Estimate'], 
                  t= LM1$coefficients['y','t value'],
                  p = LM1$coefficients['y','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}

write.csv(dfResult,'Result/A29_tfh_results.csv',row.names = F)

dfResult = dfResult%>%
  filter(grepl('olink|GO',gene))%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)%>%
  select(-nc)%>%
  arrange(p)

write.csv(dfResult,'Result/A29_tfh_results.csv',row.names = F)


##### plot GO #####
dfPlot = dfResult%>%filter(pAdj<0.05)%>%
  filter(grepl('GO',gene))

dfPlot = dfPlot%>%
  mutate(gene = gsub('olink_|GO_','',gene))%>%
  mutate(gene = gsub('_',' ',gene))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})%>%
  top_n(10,-pAdj)

p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A29_tfh_bar_GO.pdf',width = 5,height = 3)
plot(p)
dev.off()

##### plot olink #####
dfPlot = dfResult%>%filter(pAdj<0.05)%>%
  filter(grepl('olink',gene))

dfPlot = dfPlot%>%
  mutate(gene = gsub('olink_|GO_','',gene))%>%
  mutate(gene = gsub('_.*',' ',gene))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})%>%
  top_n(10,-pAdj)

p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A29_tfh_bar_olink.pdf',width = 3,height = 3)
plot(p)
dev.off()

p = ggplot(df1, aes(y=tcell,x = `olink_MCP-3_P80098_OID00474`))+
  geom_smooth(method = 'lm')+
  geom_point(aes(color=Stim))+theme_bw()+facet_wrap(~Info_day)
pdf('Result/A29_tfh_MCP.pdf',width = 5,height = 3)
plot(p)
dev.off()

t1 = df1%>%filter(Info_day==0)
cor.test(t1$tcell,t1$`olink_MCP-3_P80098_OID00474`,method = 'spearman')

t1 = df1%>%filter(Info_day==5)
cor.test(t1$tcell,t1$`olink_MCP-3_P80098_OID00474`,method = 'spearman')


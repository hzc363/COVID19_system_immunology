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

##### combined with Ab data ######
neut_ab = read.csv('Data/Antibodies/D28_M7 Neut IgG Binding.csv')
df_combine = neut_ab%>%
  dplyr::select(Info_id=ID,D28_IgG=d28_IgG_AUC,D28_Neut,
                M7_IgG = M7_IgG_AUC,M7_Neut)%>%
  mutate(D28_IgG = log2(D28_IgG+1),
         D28_Neut=log2(D28_Neut+1),
         M7_IgG=log2(M7_IgG+1) ,
         M7_Neut=log2(M7_Neut+1))%>%
  mutate(Info_id = as.character(Info_id))
RNAseq = inner_join(df_combine,RNAseq,by='Info_id')

##### correlate with Ab #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')

plotDf = df1%>%dplyr::select(Info_id,Info_tx,D28_IgG)%>%unique()
p = ggplot(plotDf,aes(x=Info_tx,y =D28_IgG ))+
  geom_boxplot()+theme_bw()
pdf('Result/A14_tx_vs_D28_IgG.pdf',3,3)
plot(p)
dev.off()

dfResult = data.frame()
for (i in grep('^GO_|^olink_|Gene_|EVE_',colnames(df1))) { 
  t1 = df1[,c("D28_IgG",'Info_onset','olink_DDX58_O95786_OID01018',
              colnames(df1)[i])]%>%na.omit()
  colnames(t1)[4]='y'
  LM1 = lm(D28_IgG~poly(Info_onset,2)+y,t1)
  LM1 = summary(LM1)
  
  LM2 = lm(olink_DDX58_O95786_OID01018~poly(Info_onset,2)+y
           ,t1)
  LM2 = summary(LM2)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['y','Estimate'], 
                  t= LM1$coefficients['y','t value'],
                  p = LM1$coefficients['y','Pr(>|t|)'],
                  p_rig = LM2$coefficients['y','Pr(>|t|)'],
                  t_rig = LM2$coefficients['y','t value'])
  dfResult = rbind(dfResult,t1)
}

write.csv(dfResult,'Result/A15_Ab_M4_results_new.csv',row.names = F)

dfResult = dfResult%>%
  #filter(grepl('olink|EVE',gene))%>%
  filter(grepl('GO_',gene))%>%
  filter(!grepl('TBA',gene))%>%
  #mutate(type = gsub('_.*','',gene))%>%
  #group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)
  

protCL=read.csv('Result/A11_gene_CL.csv')
dfPlot = dfResult%>%filter(pAdj<0.05)
dfPlot = dfPlot%>%#left_join(protCL,by=c('gene'='prot'))%>%
  mutate(pAdjRig = p.adjust(p_rig,'fdr'))%>%
  mutate(gene = gsub('olink_|EVE_|GO_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})%>%
  group_by(t>0)%>%
  top_n(5,-pAdj)

p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  #geom_point(aes(x=max(t)+1,color = pAdj>0.05),shape = 15)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A15_Ab_M4_bar_new.pdf',width =5,height = 3)
plot(p)
dev.off()

dfPlot$t_rig[dfPlot$t_rig>100]=0
p = ggplot(dfPlot,
           aes(x=t_rig,y=reorder(gene,t)))+
  geom_point(aes(x=max(t_rig)+1,color = pAdjRig>0.05),shape = 15)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A15_Ab_M4_bar_new_rig.pdf',width =3,height = 3)
plot(p)
dev.off()

pValue = dfResult%>%
  filter(pAdj<0.01)
##### ab vs time #####
df1 = RNAseq%>%
  mutate(Info_onset = 0-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  dplyr::select(D28_IgG,Info_onset)%>%
  unique()


p = ggplot(df1,aes(x=Info_onset,y= D28_IgG))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A14_time_vs_IgG_new.pdf',3,3)
plot(p)
dev.off()
cor.test(df1$D28_IgG,df1$Info_onset)


##### plot severity vs Ab #####

df1 = RNAseq%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                                   levels =c('Asymptomatic',
                                             "Moderate Symptomatic",
                                             "Severe")))

p = ggplot(df1,aes(x=Info_severitycat,y=D28_IgG))+
  geom_boxplot()+theme_bw()+geom_jitter(width = 0.2)
pdf('Result/A15_Ab_vs_severity.pdf',width = 3,height = 3)
plot(p)
dev.off()
summary.aov(lm(D28_IgG ~ Info_severitycat, data = df1))


p = ggplot(df1,aes(x=D28_Neut,y=D28_IgG))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
pdf('Result/A15_IgG_vs_Neut.pdf',width = 3,height = 3)
plot(p)
dev.off()

cor.test(df1$D28_Neut,df1$D28_IgG,method = 'spearman',
         use = 'complete')

##### DDX58 vs Ab #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  mutate(DDX58=olink_DDX58_O95786_OID01018)%>%
  dplyr::select(Info_onset,DDX58,D28_IgG)%>%
  na.omit()%>%
  mutate(DDX58 = lm(DDX58~poly(Info_onset,2))$residuals)
p = ggplot(df1,aes(x=DDX58,y=D28_IgG))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A15_DDX58_vs_D28_IgG.pdf',3,3)
plot(p)
dev.off()

df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  mutate(IL9=`EVE_IL-9`)%>%
  dplyr::select(Info_onset,IL9,D28_IgG)%>%
  na.omit()%>%
  mutate(IL9 = lm(IL9~poly(Info_onset,2))$residuals)
p = ggplot(df1,aes(x=IL9,y=D28_IgG))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A15_IL9_vs_D28_IgG.pdf',3,3)
plot(p)
dev.off()

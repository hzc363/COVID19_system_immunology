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

##### combined with T cell data ######
t1 = "_CD107a$|_IFNy$|_IL10$|_IL21$|_TNF$" 
t1 = which(grepl(t1,colnames(tcell)) )#& grepl('Memory',colnames(tcell))
tcellSub = tcell[,c(1,2,t1)]
colnames(tcellSub)=gsub('_bgsub','',colnames(tcellSub))
cor1 = cor(tcellSub[,-c(1:2)],use = 'pairwise.complete.obs',method = 'spearman')
heatmap(cor1)


###### combined with T cell data and Ab data #####
ab28 = ab1%>%filter(day ==28)%>%dplyr::select(-day)
colnames(ab28)[2]='IgG_D28'
ab120 = ab1%>%filter(day ==120)%>%dplyr::select(-day)
colnames(ab120)[2]='IgG_M4'
tcellSub = tcellSub%>%left_join(ab28,by= 'id')%>%left_join(ab120,by= 'id')

neut_ab = read.csv('Data/Antibodies/prelim %Neut at 1-50 for Pras.csv')
neut_ab = gather(neut_ab,key='day',value = "neut_value",-id)

neut_ab = neut_ab%>%
  mutate(day=gsub('D','',day))%>%
  mutate(day = as.integer(day))%>%
  mutate(id = as.character(id))%>%
  dplyr::select(-day)
tcellSub = tcellSub%>%left_join(neut_ab,by= 'id')

cor1 = cor(tcellSub[,-c(1:2)],use = 'pairwise.complete.obs',method = 'spearman')


# concise heatmap 
tcellSmall = tcellSub[,!grepl('^MN|Memory',colnames(tcellSub))]
cor1 = cor(tcellSmall[,-c(1:2)],use = 'pairwise.complete.obs',method = 'spearman')


tcellSmall = tcellSub[,grepl('CD4',colnames(tcellSub))]
plot(tcellSmall)
tcellSmall = tcellSmall%>%
  gather(key='cell',value = 'value')
p = ggplot(tcellSmall,aes(x=value))+
  geom_histogram()+theme_bw()+
  facet_wrap(~cell,scales = 'free')
plot(p)


##### correlate with t cell #####
df1 = tcell%>%
  dplyr::select(Info_id=id,
                mCD8_CD107a=`SS1_bgsub_Memory CD4p_IFNy`)%>%
  inner_join(RNAseq,by='Info_id')%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  filter(Info_tx=='placebo')

df1$mCD8_CD107a[df1$mCD8_CD107a<0]=0

plotDf = df1%>%dplyr::select(Info_id,Info_tx,mCD8_CD107a)%>%unique()
p = ggplot(plotDf,aes(x=Info_tx,y =mCD8_CD107a ))+
  geom_boxplot()+theme_bw()

dfResult = data.frame()
for (i in grep('GO_|^olink_|Gene_|EVE_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c("mCD8_CD107a","Info_onset",'olink_DDX58_O95786_OID01018',
              colnames(df1)[i])]%>%na.omit()
  colnames(t1)[4]='y'
  LM1 = lm(mCD8_CD107a~poly(Info_onset,2)+y,t1)
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

write.csv(dfResult,'Result/A14_t_cell_results_control_only.csv',row.names = F)

dfResultAll = read.csv('Result/A14_t_cell_results.csv')
dfResultAll = dfResultAll%>%inner_join(dfResult,by='gene')%>%
  dplyr::select(gene,tAll = t.x,tControl=t.y)%>%
  mutate(type = gsub('_.*','',gene))%>%
  filter(type %in%c('GO','olink'))

p = ggplot(dfResultAll,aes(x=tAll,y = tControl))+
  geom_point()+theme_bw()+
  facet_wrap(~type,scales = 'free')+
  geom_smooth(method = 'lm')+
  xlab('t stat using full dataset')+
  ylab('t stat using control arm')
pdf('Result/A14_T cell_correlation_arms.pdf',width = 8,height = 5)  
print(p)
dev.off()

dfResultAll = dfResultAll%>%group_by(type)%>%
  summarise(cor = cor(tAll,tControl,method = 'spearman'))
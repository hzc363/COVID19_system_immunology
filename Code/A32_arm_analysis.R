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
library(PRROC)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
#t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
#       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
#geneDf = read.csv('A27_top_go_genes.csv')
#t1 = paste0('Gene_',geneDf$gene_symbol)
#RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))|colnames(RNAseq)%in%t1]
RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))]

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
df1= left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

dfResult = data.frame()
for (i in grep('^GO_|^olink_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c("Info_tx","Info_onset" ,colnames(df1)[i])]%>%na.omit()
  colnames(t1)[3]='y'
  LM1 = lm(y~poly(Info_onset,2)+Info_tx,t1)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = 0-LM1$coefficients['Info_txplacebo','Estimate'], 
                  t= 0-LM1$coefficients['Info_txplacebo','t value'],
                  p = LM1$coefficients['Info_txplacebo','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}

write.csv(dfResult,'Result/A32_arm_effect.csv',row.names = F)

dfResult = dfResult%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)%>%
  filter(!grepl('TBA',gene))%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))

p = ggplot(dfResult,aes(x=beta,y = -log10(p),color = pAdj>0.05))+
  geom_point()+theme_bw()+facet_wrap(~type,ncol = 2)
pdf('Result/A32_arm_volcano.pdf',width = 8,height = 5)
plot(p)
dev.off()

##### goodness of fit #####
dfResult = data.frame()
for (i in grep('^GO_|^olink_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c("Info_tx","Info_onset" ,colnames(df1)[i])]%>%na.omit()
  colnames(t1)[3]='y'
  
  LM0 = lm(y~1,t1)
  LM0 =BIC(LM0)
  
  LM1 = lm(y~poly(Info_onset,1),t1)
  LM1 =BIC(LM1)
  
  LM2 = lm(y~poly(Info_onset,2),t1)
  LM2 = BIC(LM2)
  
  LM3 = lm(y~poly(Info_onset,3),t1)
  LM3 = BIC(LM3)
  
  LM4 = lm(y~poly(Info_onset,4),t1)
  LM4 = BIC(LM4)
  
  t1 = data.frame(gene = colnames(df1)[i],
                  "0" = LM0-LM0,
                  "1" = LM1-LM0,
                  "2"=  LM2-LM0,
                  "3" = LM3-LM0,
                  "4" = LM4-LM0)
  dfResult = rbind(dfResult,t1)
}

pdf('Result/A32_BIC.pdf',5,5)
boxplot(dfResult[,-1],xlab="highest order polynomials",ylab='BIC')

long1 = dfResult%>%gather(key = 'order',value = 'BIC',X0:X4)%>%
  mutate(order = gsub('X','',order))
p = ggplot(long1,aes(x=order,y = BIC))+
  geom_boxplot()+
  geom_jitter(width=0.2)+theme_bw()+
  xlab('highest polynomial order')
plot(p)
dev.off()




  
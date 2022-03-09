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
library(xCell)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
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

#colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
#eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
#eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
#RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### plot survival plot for virus ######
# from onset
df1 = RNAseq%>%
  mutate(Info_time2prime = Info_time2prime-Info_onset)%>%
  mutate(Info_time2res = Info_time2res-Info_onset)%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
df1$Info_time2prime[df1$Info_time2prime<0]=0

fit = survfit(Surv(Info_time2prime,Info_prime) ~ Info_sex,
              data = df1)
fit = ggsurvplot(fit)
pdf('Result/A16_survival_prime.pdf',4,4)
fit$plot
dev.off()

fit = survfit(Surv(Info_time2res,Info_resolution) ~ Info_sex,
              data = df1)
fit = ggsurvplot(fit)
pdf('Result/A13_survival_res.pdf',4,4)
fit$plot
dev.off()


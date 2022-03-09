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
library(corrplot)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
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

#colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
#eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
#eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
#RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### combined with Ab data ######
neut_ab = read.csv('Data/Antibodies/prelim %Neut at 1-50 for Pras.csv')
neut_ab = gather(neut_ab,key='day',value = "neut_value",-id)
neut_ab = neut_ab%>%
  mutate(day=gsub('D','',day))%>%
  mutate(day = as.integer(day))%>%
  mutate(id = as.character(id))
df_combine = full_join(neut_ab,ab1,by=c('day','id'))
df_combine = df_combine%>%
  filter(day==28)%>%
  filter(!is.na(neut_value)|!is.na(IgG_antibody_value))%>%
  dplyr::select(-day)%>%
  mutate(id = as.character(id))
colnames(df_combine)=c('Info_id',"Neut_M4",'IgG_M4')
RNAseq = inner_join(df_combine,RNAseq,by='Info_id')

#### sex vs ab #####
sex1 = RNAseq[,c('Info_id','Info_sex','Info_onset')]%>%unique()
ab_sex = ab1%>%inner_join(sex1,by=c('id'='Info_id'))%>%
  mutate(Info_onset = day-Info_onset)%>%
  mutate(Info_onset_week = round(Info_onset/7))
p = ggplot(ab_sex,aes(y=IgG_antibody_value,x = Info_onset,color = factor(Info_sex)))+
  geom_point()
plot(p)

p = ggplot(ab_sex,aes(y=IgG_antibody_value,x = factor(day),color = factor(Info_sex)))+
  geom_boxplot()
plot(p)

p = ggplot(ab_sex,aes(y=IgG_antibody_value,x = factor(Info_onset_week),
                      color = factor(Info_sex)))+
  geom_boxplot()
plot(p)


##### correlate with Ab #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')

dfResult = data.frame()
for (i in grep('^GO_|^olink_|Gene_',colnames(df1))) { 
  t1 = df1[,c(2,18,i)]%>%na.omit()
  colnames(t1)[3]='y'
  LM1 = lm(Neut_M4~poly(Info_onset,2)+y,t1)
  #LM1 = lm(IgG_M4~y,t1)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['y','Estimate'], 
                  t= LM1$coefficients['y','t value'],
                  p = LM1$coefficients['y','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}

dfResult = dfResult%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))

write.csv(dfResult,'Result/A15_Ab_M4_results.csv',row.names = F)


protCL=read.csv('Result/A11_gene_CL.csv')
dfPlot = dfResult%>%filter(pAdj<0.05)
dfPlot = dfPlot%>%left_join(protCL,by=c('gene'='prot'))
p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  geom_point(aes(x=7,color = CL),shape = 15,size = 5)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A15_Ab_M4_bar.pdf',width =7,height = 7)
plot(p)
dev.off()

pValue = dfResult%>%
  filter(pAdj<0.01)
##### ab vs time #####
df1 = RNAseq%>%
  mutate(Info_onset = 30-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  dplyr::select(Neut_M4,Info_onset)%>%
  unique()


p = ggplot(df1,aes(x=Info_onset,y= Neut_M4))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A14_time_vs_neut.pdf',3,3)
plot(p)
dev.off()
##### plot trajectory #####
df1 = RNAseq%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                                   levels =c('Asymptomatic',
                                             "Moderate Symptomatic",
                                             "Severe")))%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')

df1 = df1[,c('Info_onset','IgG_M4',pValue$gene)]%>%
  gather(key='prot',value = 'value',-Info_onset,-IgG_M4)%>%
  filter(Info_onset>0 & Info_onset<20)%>%
  mutate(IgG_M4_binary=IgG_M4>median(IgG_M4,na.rm = T))%>%
  na.omit()


pdf('Result/A15_symp_timeline.pdf',width = 20,height = 20)
p = ggplot(df1,aes(x=Info_onset,y=value,color = IgG_M4_binary))+
  geom_point(size=1)+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,1),se=F)+
  facet_wrap(~prot,scales = 'free')
plot(p)

p = ggplot(df1,aes(x=IgG_M4,y=value))+
  geom_point()+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,1),se=F)+
  facet_wrap(~prot,scales = 'free')
plot(p)
dev.off()

##### plot severity vs Ab #####

df1 = RNAseq%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                                   levels =c('Asymptomatic',
                                             "Moderate Symptomatic",
                                             "Severe")))
  
p = ggplot(df1,aes(x=Info_severitycat,y=Neut_M4,color = Info_severitycat))+
  geom_boxplot()+theme_bw()+geom_jitter(width = 0.2)
pdf('Result/A15_Ab_vs_severity.pdf',width = 5,height = 3)
plot(p)
dev.off()

##### plot onset vs Ab #####

df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
  
p = ggplot(df1,aes(x=Info_onset,y=IgG_M4))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)

p = ggplot(df1,aes(x=`GO_V(D)J_recombination`,y=IgG_M4))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)


##### correlate with Ab res #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')

dfResult = data.frame()
for (i in grep('^GO_|^olink_',colnames(df1))) { 
  t1 = df1[,c(3,18,i)]%>%na.omit()
  colnames(t1)[3]='y'
  t1 = t1 %>%
    mutate(lm(y~poly(Info_onset,2))$residuals)%>%
    mutate(IgG_M4 = lm(IgG_M4~poly(Info_onset,2))$residuals)
  LM1 = lm(IgG_M4~y,t1)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['y','Estimate'], 
                  t= LM1$coefficients['y','t value'],
                  p = LM1$coefficients['y','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}

dfResult = dfResult%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))

write.csv(dfResult,'Result/A15_Ab_M4_results_res.csv',row.names = F)


protCL=read.csv('Result/A11_gene_CL.csv')
dfPlot = dfResult%>%filter(pAdj<0.05)
dfPlot = dfPlot%>%left_join(protCL,by=c('gene'='prot'))
p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  geom_point(aes(x=5,color = CL),shape = 15,size = 5)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A15_Ab_M4_bar_res.pdf',width =10,height = 12)
plot(p)
dev.off()


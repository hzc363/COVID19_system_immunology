##### load packages #####
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
library(pROC)

##### load data #####
load('Result/A01_compiled_data.rda')

favi_olink = fread("Data/Favi Data/OLINK/20212503_Jagannathan_NPX.csv",data.table = F)
colnames(favi_olink)[-1] = gsub('\\-| ','',colnames(favi_olink)[-1])
favi_clinic = fread('Data/Favi Data/Clinical Data_ Favi.csv',data.table = F)
load('Result/A20_final_models.rda')

favi_tcell = fread('Data/T cell data_favi.csv',data.table = F)
colnames(olink) = gsub('\\-| ','',colnames(olink))

##### normalize data #####
prot1 = colnames(favi_olink)[grepl('_OID',colnames(favi_olink))]
for (o1 in prot1) {
  m1 = mean(olink[,o1])
  s1 = sd(olink[,o1])
  favi_olink[o1] = scale(favi_olink[o1])*s1+m1
}

colnames(favi_olink)[-1]=paste0('olink_',(colnames(favi_olink)[-1]))

##### prepare data #####
favi_olink = favi_clinic%>%
  dplyr::select(ID =participantId,Info_age=age_c,Info_Gender=sex,
         Ab1=VN_titer28,severity=ED_day1)%>%
  inner_join(favi_olink,by='ID')%>%
  filter(!is.na(Info_age))
t1 = !is.na(favi_olink$severity)
favi_olink$severity[t1]=1
favi_olink$severity[!t1]=0

favi_olink$Ab1 = as.numeric(gsub('<','',favi_olink$Ab1 ))
favi_olink$Ab1 = log2(favi_olink$Ab1)
favi_olink$Ab1_bi = favi_olink$Ab1 >median(favi_olink$Ab1 ,na.rm = T)

favi_olink$ID = gsub('.*\\-','',favi_olink$ID)%>%as.integer()

favi_olink =favi_tcell%>%
  filter(Stim%in%c('MN','SS1'))%>%
  mutate(mCD4=`MCD4/IFNg+IL21-TNF+`+`MCD4/IFNg+IL21-TNF-`)%>%
  dplyr::select(ID=id,mCD4)%>%
  filter(!is.na(mCD4))%>%
  group_by(ID)%>%summarise(mCD4=sum(mCD4))%>%
  full_join(favi_olink ,by='ID')

favi_olink$mCD4_bi = favi_olink$mCD4 >median(favi_olink$mCD4 ,na.rm = T)

##### predict #####
pred1 = RFlist
for (r1 in names(RFlist)) {
  p1 = predict(RFlist[[r1]],favi_olink,type='prob')[,2]
  pred1[[r1]]=p1
}

boxplot(pred1$D28_IgG~favi_olink$Ab1_bi)
boxplot(pred1$Info_severitycat~favi_olink$severity)
boxplot(pred1$mCD4~favi_olink$mCD4_bi)


pdf('Result/A31_AUC_favi.pdf',height = 5,width = 7)

roc1=pROC::roc(favi_olink$Ab1_bi,pred1$D28_IgG)
plot(roc1)
print(roc1)

roc1=pROC::roc(favi_olink$severity,pred1$Info_severitycat)
plot(roc1)
print(roc1)

roc1=pROC::roc(favi_olink$mCD4_bi,pred1$mCD4)
plot(roc1)
print(roc1)

dev.off()

##### compute effect size #####
f1 = function(x,y){
  df = data.frame(x=x,y=y)
  fit = summary(lm(y~x))
  return(fit$coefficients['x','t value'])
}
long1 = favi_olink%>%
  dplyr::select(ID, Ab1 ,severity,mCD4,
         olink_PPP1R9B_Q96SB3_OID00936:olink_CSF1_P09603_OID00562)%>%
  gather(key='olink',value = 'expression',
         olink_PPP1R9B_Q96SB3_OID00936:olink_CSF1_P09603_OID00562)%>%
  gather(key='clinical',value='clinical_value',Ab1:mCD4)%>%
  group_by(clinical,olink)%>%
  summarise(t_value = f1(x=expression,y=clinical_value))

##### compare effect size #####
severity1 = read.csv('Result/A12_severity_results.csv')
clear1 = read.csv('Result/A13_clearance_results_auc.csv')
tcell1 = read.csv('Result/A14_t_cell_results.csv')
Ab1 = read.csv('Result/A15_Ab_M4_results_new.csv')
Ab2 = read.csv('Result/A15_Ab_M7_results_new.csv')

lambda1 = severity1%>%filter(grepl('olink',gene))%>%
  dplyr::select(gene,t)%>%mutate(clinical='severity')
Ab1 = Ab1%>%filter(grepl('olink',gene))%>%
  dplyr::select(gene,t)%>%mutate(clinical='Ab1')
mCD4 = tcell1 %>%filter(grepl('olink',gene))%>%
  dplyr::select(gene,t)%>%mutate(clinical='mCD4')
lambda1 = rbind(lambda1,Ab1)
lambda1 = rbind(lambda1,mCD4)


compare1 = lambda1%>%inner_join(long1,by=c('clinical','gene'='olink'))
p = ggplot(compare1,aes(x=t,y=t_value))+
  geom_point()+
  geom_smooth(method = 'lm')+theme_bw()+
  facet_wrap(~clinical,scales = 'free')
pdf('Result/A31_correlation_favi.pdf',height = 5,width = 10)
plot(p)
dev.off()

compare1%>%group_by(clinical)%>%
  summarise(cor = cor(t,t_value),p = cor.test(t,t_value)$p.value	)

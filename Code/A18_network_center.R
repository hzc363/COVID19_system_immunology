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

colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### olink data #####
RNAseq = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')

olink = colnames(RNAseq)[grep('^olink_|^EVE_|GO_',colnames(RNAseq))]
olink = expand.grid('olink_DDX58_O95786_OID01018' ,olink )
olink$Var1 = as.character(olink$Var1)
olink$Var2 = as.character(olink$Var2)
olink$p = NA; olink$t = NA
for (i in 1:nrow(olink)) {
  t1 = data.frame(y=RNAseq[,olink$Var1[i]],
                  x=RNAseq[,olink$Var2[i]],
                  time = RNAseq[,'Info_onset'])
  t1 = lm(y~poly(time,2)+x,t1)
  t1 = summary(t1)
  olink$p[i]=t1$coefficients['x','Pr(>|t|)']
  olink$t[i]=t1$coefficients['x','t value']
  
}
olink$p[olink$Var1==olink$Var2] = 1
olink$p = p.adjust(olink$p,'fdr')
connect1 = olink%>%
  mutate(sig = p<0.05)%>%
  group_by(Var1)%>%
  summarise(N = sum(sig))

t1 = olink%>%filter(grepl('DDX58',Var1))%>%
  filter(p<0.05)


##### GO data #####
RNAseq = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')

olink = colnames(RNAseq)[grep('^GO_',colnames(RNAseq))]
olink = olink[nchar(olink)<40]
olink = expand.grid('GO_RIG-I_signaling_pathway' ,olink )
olink$Var1 = as.character(olink$Var1)
olink$Var2 = as.character(olink$Var2)
olink$p = NA; olink$t=NA
for (i in 1:nrow(olink)) {
  t1 = data.frame(y=RNAseq[,olink$Var1[i]],
                  x=RNAseq[,olink$Var2[i]],
                  time = RNAseq[,'Info_onset'])
  t1 = lm(y~poly(time,2)+x,t1)
  t1 = summary(t1)
  olink$p[i]=t1$coefficients['x','Pr(>|t|)']
  olink$t[i]=t1$coefficients['x','t value']
}
olink$p[olink$Var1==olink$Var2] = 1
olink$p = p.adjust(olink$p,'fdr')
connect1 = olink%>%
  mutate(sig = p<0.05)%>%
  group_by(Var1)%>%
  summarise(N = sum(sig))

t1 = olink%>%filter(grepl('RIG-I',Var1))%>%
  filter(p<0.01)

##### rig vs go #####
dfPlot = RNAseq[,c('olink_DDX58_O95786_OID01018',
                   'GO_RIG-I_signaling_pathway',
                   'Info_onset','Info_severitycat',
                   'Info_day')]%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  na.omit()%>%
  mutate(`GO_RIG-I_signaling_pathway`=
           lm(`GO_RIG-I_signaling_pathway`~poly(Info_onset,2))$residuals)
  
p = ggplot(dfPlot,aes(x=olink_DDX58_O95786_OID01018,
                      y=`GO_RIG-I_signaling_pathway`))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
plot(p)
cor.test(dfPlot$olink_DDX58_O95786_OID01018,
    dfPlot$`GO_RIG-I_signaling_pathway`)

##### rig vs MPC-3 #####
dfPlot = RNAseq[,c('olink_DDX58_O95786_OID01018',
                   'olink_MCP-3_P80098_OID00474',
                   'Info_onset','Info_severitycat',
                   'Info_day')]%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  na.omit()%>%
  mutate(`olink_MCP-3_P80098_OID00474`=
           lm(`olink_MCP-3_P80098_OID00474`~poly(Info_onset,2))$residuals)

p = ggplot(dfPlot,aes(x=olink_DDX58_O95786_OID01018,
                      y=`olink_MCP-3_P80098_OID00474`))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
plot(p)
cor.test(dfPlot$olink_DDX58_O95786_OID01018,
         dfPlot$`olink_MCP-3_P80098_OID00474`)

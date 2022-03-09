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

##### plot survival plot ######

# from onset
df1 = RNAseq%>%
  mutate(Info_time2prime = Info_time2prime-Info_onset)%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
df1$Info_time2prime[df1$Info_time2prime<0]=0

fit = survfit(Surv(Info_time2prime,Info_prime) ~ Info_severitycat,
              data = df1)
fit = ggsurvplot(fit)

pdf('Result/A13_survival_prime.pdf',4,4)
fit$plot
dev.off()

##### correlate with viral clearence #####
df1 = RNAseq%>%
  mutate(Info_time2prime = Info_time2prime-Info_onset)%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
df1 = cbind(Info_surv = Surv(df1$Info_time2prime,df1$Info_prime),
            df1)

dfResult = data.frame()
for (i in grep('^GO_|^IgG_|^olink_|Gene_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c(1,16,i)]%>%na.omit()
  colnames(t1)[3]='y'
  LM1 = coxph(Info_surv~poly(Info_onset,2)+y,t1)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['y','coef'], 
                  t= LM1$coefficients['y','z'],
                  p = LM1$coefficients['y','Pr(>|z|)'])
  dfResult = rbind(dfResult,t1)
}

dfResult = dfResult%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))

write.csv(dfResult,'Result/A13_clearance_results.csv',row.names = F)


protCL=read.csv('Result/A11_gene_CL.csv')
dfPlot = dfResult%>%filter(pAdj<0.05)
dfPlot = dfPlot%>%left_join(protCL,by=c('gene'='prot'))
p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  geom_point(aes(x=5,color = CL),shape = 15,size = 5)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A13_prime_bar.pdf',width = 5,height = 2)
plot(p)
dev.off()

pValue = dfResult%>%
  filter(pAdj<0.01)

##### plot trajectory #####
df1 = RNAseq%>%
  mutate(Info_time2prime = Info_time2prime-Info_onset)%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
df1 = cbind(Info_surv = Surv(df1$Info_time2prime,df1$Info_prime),
            df1)

dfPlot = df1[,c('Info_id','Info_time2prime','Info_prime','Info_onset',
                pValue$gene)]%>%
  gather(key='prot',value = 'value',-Info_time2prime,
         -Info_prime,-Info_onset,-Info_id)%>%
  na.omit()%>%
  group_by(prot)%>%
  mutate(value = lm(value~poly(Info_onset,2))$residuals)%>%
  mutate(value = value>median(value,na.rm = T))


fit2 <- survfit( Surv(Info_time2prime,Info_prime) ~ value+prot, data = dfPlot )
ggsurv <- ggsurvplot(fit2, conf.int = F)
surv_pvalue(fit2)

pdf('Result/A13_prime_survival.pdf',width = 5,height = 5)
plot(ggsurv$plot+facet_wrap(~prot,scales = 'free'))
dev.off()

dfPlot = df1[,c('Info_time2prime','Info_prime','Info_onset',pValue$gene)]%>%
  gather(key='prot',value = 'value',-Info_time2prime,-Info_prime,-Info_onset)%>%
  na.omit()%>%
  mutate(cleared = (Info_onset>Info_time2prime)&(Info_prime==1))%>%
  filter(Info_prime==1)%>%
  mutate(clear_cat = Info_time2prime>median(Info_time2prime))%>%
  group_by(prot)#%>%
  #mutate(value = lm(value~poly(Info_onset,2))$residuals)

p = ggplot(dfPlot,aes(x=Info_onset,y=value,color = clear_cat))+
  geom_point()+
  facet_wrap(~prot,scales = "free")+
  theme_bw()+geom_smooth(method = 'lm')
plot(p)

##### plot trajectory lab#####
df1 = RNAseq%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                                   levels =c('Asymptomatic',
                                             "Moderate Symptomatic",
                                             "Severe")))%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')

df1 = df1[,c('Info_onset','Info_severitycat',
             'Lab_pct_neutro_v','Lab_pct_lymph_v')]%>%
  gather(key='prot',value = 'value',-Info_onset,-Info_severitycat)%>%
  filter(Info_onset>0 & Info_onset<20)%>%
  mutate(Info_onset_cat = factor(ceiling(Info_onset/7),
                                 label = c('week 1','week 2','week 3')))


pdf('Result/A12_symp_timeline_lab.pdf',width = 7,height = 3)
p = ggplot(df1,aes(x=Info_onset,y=value,color = Info_severitycat))+
  geom_point(size=1)+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,1),se=F)+
  facet_wrap(~prot,scales = 'free',ncol = 2)
plot(p)
dev.off()

##### plot boxplot #####
df1 = RNAseq[,c('Info_severitycat',pValue$gene)]%>%
  gather(key='prot',value = 'value',-Info_severitycat)

pdf('Result/A12_symp_boxplot.pdf',width = 4,height = 12)
p = ggplot(df1,aes(x=Info_severitycat,y=value,color = Info_severitycat))+
  geom_boxplot(outlier.shape = NA)+theme_bw()+
  geom_jitter(width = 0.2,alpha=0.3)+
  facet_wrap(~prot,scales = 'free',ncol = 1)
plot(p)
dev.off()

df1 = RNAseq[,c('Info_severitycat','Lab_pct_neutro_v','Lab_pct_lymph_v')]%>%
  gather(key='prot',value = 'value',-Info_severitycat)

pdf('Result/A12_symp_boxplot_lab.pdf',width = 7,height = 3)
p = ggplot(df1,aes(x=Info_severitycat,y=value,color = Info_severitycat))+
  geom_boxplot(outlier.shape = NA)+theme_bw()+
  geom_jitter(width = 0.2,alpha=0.3)+
  facet_wrap(~prot,scales = 'free',ncol = 2)
plot(p)
dev.off()

##### severity vs viral clearance #####
severity1 = read.csv('Result/A12_severity_results.csv')
clear1 = read.csv('Result/A13_clearance_results.csv')

dfPlot = inner_join(severity1,clear1,by='gene')
dfPlot$sig = "other"
dfPlot$sig[dfPlot$pAdj.x<0.05 & dfPlot$pAdj.y<0.05] = "both"
dfPlot$sig[dfPlot$pAdj.x<0.05 & dfPlot$pAdj.y>0.05] = "severe only"


p = ggplot(dfPlot,aes(t.x,t.y))+
  geom_point(aes(color=sig))+
  ylim(-5,5)+xlim(-5,5)+
  theme_bw()+
  xlab('t stat with severity')+
  xlab('t stat with viral clearance')
  #geom_smooth(method = 'lm')
plot(p)

##### heapmap #####
heatPlot = dfPlot%>%
  filter(pAdj.x<0.05)%>%
  filter(grepl('olink',gene))%>%
  arrange(t.x,t.y)%>%
  mutate(gene=gsub('olink_','',gene))%>%
  mutate(gene=gsub('_.*','',gene))

t1 = as.matrix(heatPlot[,c('t.x','t.y')])
rownames(t1)=heatPlot$gene
p1 = as.matrix(heatPlot[,c('pAdj.x','pAdj.y')])
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
pdf('Result/A13_heatmap_symp_vs_clearance.pdf',width =10,height = 5)
corrplot::corrplot(t(t1),is.corr = F,sig.level = 0.05,
                   p.mat = t(p1),insig = "label_sig",
                   cl.pos = "b",cl.ratio = 1, col = rev(col2(50)))
dev.off()




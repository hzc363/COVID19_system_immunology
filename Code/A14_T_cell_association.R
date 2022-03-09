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
pdf('Result/A14_Tcell_correlation.pdf',7,7)
corrplot(cor1,order = 'hclust')
dev.off()

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
pdf('Result/A14_Tcell_Ab_correlation.pdf',10,10)
corrplot(cor1,order = 'hclust')
dev.off()

# concise heatmap 
tcellSmall = tcellSub[,!grepl('^MN|Memory',colnames(tcellSub))]
cor1 = cor(tcellSmall[,-c(1:2)],use = 'pairwise.complete.obs',method = 'spearman')
pdf('Result/A14_Tcell_Ab_correlation_small.pdf',5,5)
corrplot(cor1,order = 'hclust',hclust.method = 'ward.D2')
dev.off()

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
  filter(Info_severitycat!='Asymptomatic')

df1$mCD8_CD107a[df1$mCD8_CD107a<0]=0

plotDf = df1%>%dplyr::select(Info_id,Info_tx,mCD8_CD107a)%>%unique()
p = ggplot(plotDf,aes(x=Info_tx,y =mCD8_CD107a ))+
  geom_boxplot()+theme_bw()
pdf('Result/A14_tx_vs_tcell.pdf',3,3)
plot(p)
dev.off()

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

write.csv(dfResult,'Result/A14_t_cell_results.csv',row.names = F)

dfResult = dfResult%>%
  #filter(grepl('olink|EVE',gene))%>%
  filter(grepl('GO_',gene))%>%
  #mutate(type = gsub('_.*','',gene))%>%
  #group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)



protCL=read.csv('Result/A11_gene_CL.csv')
dfPlot = dfResult%>%filter(pAdj<0.05)
dfPlot = dfPlot%>%#left_join(protCL,by=c('gene'='prot'))%>%
  mutate(pAdjRig = p.adjust(p_rig,'fdr'))%>%
  filter(!grepl('TBA',gene))%>%
  mutate(gene = gsub('olink_|EVE_|GO_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})%>%
  top_n(10,-pAdj)

p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  #geom_point(aes(x=max(t)+1,color = pAdj>0.05),shape = 15)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A14_t_cell_bar.pdf',width = 5,height = 3)
plot(p)
dev.off()

dfPlot$t_rig[dfPlot$t_rig>100]=0
p = ggplot(dfPlot,
           aes(x=t_rig,y=reorder(gene,t)))+
  geom_point(aes(x=max(t_rig)+1,color = pAdjRig>0.05),shape = 15)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A14_t_cell_bar_rig.pdf',width =5,height = 3)
plot(p)
dev.off()


pValue = dfResult%>%
  filter(pAdj<0.01)

##### t cell vs time #####
df1 = tcell%>%
  dplyr::select(Info_id=id,mCD8_CD107a=`MN_bgsub_CD8p_CD107a`)%>%
  inner_join(RNAseq,by='Info_id')%>%
  mutate(Info_onset = 30-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  dplyr::select(mCD8_CD107a,Info_onset)%>%
  unique()

df1$mCD8_CD107a[df1$mCD8_CD107a<0]=0

p = ggplot(df1,aes(x=Info_onset,y= mCD8_CD107a))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A14_time_vs_cd8.pdf',3,3)
plot(p)
dev.off()

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

##### severity vs CD8+ T cell #####
df1 = tcell%>%
  dplyr::select(Info_id=id,mCD8_CD107a=`SS1_bgsub_Memory CD8p_CD107a`)%>%
  inner_join(RNAseq,by='Info_id')%>%
  dplyr::select(Info_severitycat,mCD8_CD107a)%>%
  unique()
df1$mCD8_CD107a[df1$mCD8_CD107a<0]=0


p = ggplot(df1,aes(x=Info_severitycat,y = mCD8_CD107a))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+theme_bw()
pdf('Result/A14_mCD8_CD107a.pdf',width =3,height = 3)
plot(p)
dev.off()

kruskal.test(mCD8_CD107a ~ Info_severitycat, data = df1)

summary.aov(lm(mCD8_CD107a ~ Info_severitycat, data = df1))

##### DDX58 vs t cell #####
df1 = tcell%>%
  dplyr::select(Info_id=id,
                tcell=`SS1_bgsub_Memory CD4p_IFNy`)%>%
  inner_join(RNAseq,by='Info_id')%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  mutate(DDX58=olink_DDX58_O95786_OID01018)%>%
  dplyr::select(Info_onset,DDX58,tcell)%>%
  na.omit()%>%
  mutate(DDX58 = lm(DDX58~poly(Info_onset,2))$residuals)
p = ggplot(df1,aes(x=DDX58,y=tcell))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A15_DDX58_vs_tcell.pdf',3,3)
plot(p)
dev.off()

df1 = tcell%>%
  dplyr::select(Info_id=id,
                tcell=`SS1_bgsub_Memory CD4p_IFNy`)%>%
  inner_join(RNAseq,by='Info_id')%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  mutate(IL9=`EVE_IL-9`)%>%
  dplyr::select(Info_onset,IL9,tcell)%>%
  na.omit()%>%
  mutate(IL9 = lm(IL9~poly(Info_onset,2))$residuals)
p = ggplot(df1,aes(x=IL9,y=tcell))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A15_IL9_vs_tcell.pdf',3,3)
plot(p)
dev.off()

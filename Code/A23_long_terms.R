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
library(glmnet)
library(randomForest)
library(PRROC)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
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

longterm1 = fread('Data/ClinicalData/outcomes_def2_w_groups-2.csv',
                  data.table = F)%>%
  mutate(participant_id = as.character(participant_id))%>%
  select(Info_id=participant_id,LTgroup = `Long COVID`,def2_m4)
RNAseq = right_join(longterm1,RNAseq,
                   by=c('Info_id'))


##### correlate with long term category #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(LTgroup%in%c('Consistent long COVID', 'Consistent no long COVID',
                      'Resolving'))%>%
  mutate(LTgroup= factor(LTgroup))%>%
  mutate(LTgroup=relevel(LTgroup,ref = "Consistent long COVID"))

dfResult = data.frame()
for (i in grep('Gene_|olink_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c('LTgroup',"Info_onset",colnames(df1)[i])]%>%na.omit()
  colnames(t1)[3]='y'
  #t1$y = log2(t1$y)
  t1 = t1%>%na.omit()
  if(nrow(t1)<50){next}
  LM1 = lm(y~poly(Info_onset,2)+LTgroup,t1)
  LM2 = lm(y~poly(Info_onset,2),t1)
  aov1 = anova(LM1, LM2)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  p_anova = aov1$`Pr(>F)`[2],
                  t_long= LM1$coefficients['LTgroupConsistent no long COVID','t value'],
                  b_long= LM1$coefficients['LTgroupConsistent no long COVID','Estimate'],
                  p_long = LM1$coefficients['LTgroupConsistent no long COVID','Pr(>|t|)'],
                  t_resolve= LM1$coefficients['LTgroupResolving','t value'],
                  p_resolve = LM1$coefficients['LTgroupResolving','Pr(>|t|)'],
                  b_resolve= LM1$coefficients['LTgroupResolving','Estimate'])
  dfResult = rbind(dfResult,t1)
}


dfResult = dfResult%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj_anova = p.adjust( p_anova,'fdr'))%>%
  mutate(pAdj_long = p.adjust( p_long,'fdr'))%>%
  mutate(pAdj_resolve = p.adjust( p_resolve,'fdr'))%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)%>%select(-nc)
write.csv(dfResult,'Result/A23_long_symp_result.csv',row.names = F)

##### plot heatmap #####
dfPlot = dfResult%>%filter(pAdj_anova<0.05)
dfPlot = RNAseq[,c('LTgroup',dfPlot$gene)]
dfPlot = gather(dfPlot,key = 'gene',value='expression',-LTgroup)%>%
  filter(LTgroup%in%c('Consistent long COVID', 'Consistent no long COVID',
                      'Resolving'))
p = ggplot(dfPlot,aes(x=LTgroup,y= expression,color = LTgroup))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(width = 0.3)+
  theme_bw()+facet_wrap(~gene,scales = 'free')+coord_flip()
pdf('Result/A23_boxplots.pdf',width = 12,height=6)
plot(p)
dev.off()

##### volcano plot #####
dfPlot = dfResult%>%filter(type %in%c('Gene','GO','olink'))
df_text = dfResult%>%filter(pAdj_long<0.05)
p = ggplot(dfPlot,aes(x=-b_long, y = -log10(p_long)))+
  geom_point()+facet_wrap(~type,scales = 'free_x')+
  geom_text(data=df_text,aes(x=-b_long, y = -log10(p_long),label=gene))+
  theme_bw()+xlab('beta (long covid vs no long covid')+ylab('-log10(p)')
pdf('Result/A23_vocano_long_vs_no.pdf',height = 5,width = 10)
plot(p)
dev.off()


df_text = dfResult%>%filter(pAdj_resolve<0.05)
p = ggplot(dfPlot,aes(x=-b_resolve, y = -log10(p_resolve)))+
  geom_point()+facet_wrap(~type,scales = 'free_x')+
  geom_text(data=df_text,aes(x=-b_resolve, y = -log10(p_resolve),label=gene))+
  theme_bw()+xlab('beta (long covid vs resolving)')+ylab('-log10(p)')
pdf('Result/A23_vocano_long_vsresolving.pdf',height = 5,width = 10)
plot(p)
dev.off()

##### correlate with long term category #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  mutate(LTgroup= factor(LTgroup))%>%
  mutate(LTgroup=relevel(LTgroup,ref = "Consistent long COVID"))

dfResult = data.frame()
for (i in grep('Gene_|olink_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c('def2_m4',"Info_onset",colnames(df1)[i])]%>%na.omit()
  colnames(t1)[3]='y'
  #t1$y = log2(t1$y)
  t1 = t1%>%na.omit()
  if(nrow(t1)<50){next}
  LM1 = lm(y~poly(Info_onset,2)+def2_m4,t1)
  LM2 = lm(y~poly(Info_onset,2),t1)
  aov1 = anova(LM1, LM2)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  p_anova = aov1$`Pr(>F)`[2],
                  t_m4= LM1$coefficients['def2_m4','t value'],
                  b_m4= LM1$coefficients['def2_m4','Estimate'],
                  p_m4 = LM1$coefficients['def2_m4','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}


dfResult = dfResult%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj_anova = p.adjust( p_anova,'fdr'))%>%
  mutate(pAdj_m4 = p.adjust( p_m4,'fdr'))%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)%>%select(-nc)
write.csv(dfResult,'Result/A23_M4_result.csv',row.names = F)


##### volcano plot #####
dfPlot = dfResult%>%filter(type %in%c('Gene','GO','olink'))%>%
  mutate(fdr = pAdj_m4)%>%as.data.frame()

df_text = dfPlot%>%filter(pAdj_m4<0.05)%>%
  mutate(gene = gsub('olink_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))

p = ggplot(dfPlot,aes(x=-b_m4, y = -log10(p_m4),color= fdr<0.05))+
  geom_point()+facet_wrap(~type,scales = 'free_x')+
  geom_text(data=df_text,aes(x=-b_m4, y = -log10(p_m4)+0.2,label=gene))+
  theme_bw()+xlab('beta (M4 PASC vs control)')+ylab('-log10(p)')+
  scale_color_manual(values=c("black", "red"))
pdf('Result/A23_M4.pdf',height = 5,width = 10)
plot(p)
dev.off()


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
library(glmnet)
library(randomForest)

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

##### correlate with onset #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')

# remove random effect
dfTractory = data.frame()
days = data.frame(Info_onset=0:15)
for (i in grep('^GO_|^IgG_|^Lab_|^olink_|^xcell_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c(15,i)]%>%na.omit()
  colnames(t1)[2]='y'
  #if(n_distinct(t1$y)>1){t1$y = scale(t1$y)}
  LM1 = lm(y~poly(Info_onset,2),t1)
  LM2 = lm(y~1,t1)
  LM3 = lm(y~Info_onset,t1)
  
  pValue = anova(LM1,LM2)$`Pr(>F)`[2]
  pValue2 = anova(LM1,LM3)$`Pr(>F)`[2]
  
  pred1 = predict(object=LM1,newdata=days)
  t1 = data.frame(gene = colnames(df1)[i], 
                  onset= days,
                  predValue = pred1,
                  p = pValue,
                  p2 = pValue2,
                  mean1 = mean(t1$y),
                  var1 = var(t1$y))
  dfTractory = rbind(dfTractory,t1)
}
write.csv(dfTractory,'Result/A11_trajectory.csv',row.names = F)

dfTractory = dfTractory%>%dplyr::select(-mean1,-var1,-p2)
pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  filter(!grepl('TBA',gene))%>%
  unique()%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)

dfPlot = dfTractory%>%
  filter(gene %in% pValue$gene)%>%
  mutate(Info_onset = str_pad(as.character(Info_onset), 2, pad = "0",side = 'left'))%>%
  mutate(day=paste0('day',Info_onset))%>%
  dplyr::select(-p,-Info_onset)%>%
  spread(key = day,value=predValue)%>%
  mutate(gene = gsub('olink_|GO_|EVE_','',gene))%>%
  mutate(gene = gsub('_Q.*|_P.*|_O.*','',gene))%>%
  mutate(gene = gsub('_',' ',gene))

dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

geneHC = hclust(dist(dfPlot),method="ward.D")
geneCL = cutree(geneHC,4)%>%
  factor(labels = c('red','green','blue','orange'))%>%
  as.character()

pdf('Result/A11_heatmap_onset.pdf',width = 10,height = 10)
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          RowSideColors= geneCL[geneHC$order],
          col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'row',hclustfun=f1)
dev.off()


##### plot cytokines #####
protCL = data.frame(prot= rownames(dfPlot),CL = geneCL)
write.csv(protCL,'Result/A11_gene_CL.csv',row.names = F)

protCL = protCL%>%filter(grepl('olink',prot))
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')
df1 = df1[,c('Info_onset',protCL$prot)]%>%
  gather(key='prot',value = 'value',-Info_onset)
df1 = df1%>%inner_join(protCL,by='prot')

p = ggplot(df1,aes(x=Info_onset,y=value,color=CL))+
  geom_point()+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  facet_wrap(~prot,scales = 'free',ncol = 2)+
  scale_color_manual(values=gsub('.*_','',unique(protCL$CL)))
pdf('Result/A11_prot.pdf',width = 6,height = 10)
plot(p)
dev.off()
##### plot cluster mean #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')

dfPlot = data.frame(cluster = paste0("cluster_",geneCL),gene = rownames(dfPlot))
t1 = df1[,c('Info_id',"Info_onset",dfPlot$gene)]%>%
  gather(key='gene',value = 'value',-Info_onset,-Info_id)
dfPlot = dfPlot%>%inner_join(t1,by='gene')

dfPlot = dfPlot%>%group_by(gene)%>%
  mutate(value = scale(value))%>%
  group_by(cluster,Info_onset,Info_id)%>%
  na.omit()%>%
  summarise(value = mean(value))%>%
  filter(Info_onset>=0 & Info_onset<=15)

p = ggplot(dfPlot,aes(x=Info_onset,y=value,color = cluster, group=cluster))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm',formula = y~poly(x,2))+
  scale_color_manual(values=gsub('.*_','',unique(dfPlot$cluster)))+
  facet_wrap(~cluster,ncol = 1)
pdf('Result/A11_mean_trajectory.pdf',width = 5,height = 8)
plot(p)
dev.off()

##### major cells #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  filter(Info_onset<=15)

df1 = df1%>%
  dplyr::select(Info_onset,`xcell_B-cells`,
                #`xcell_Basophils`,xcell_Eosinophils,xcell_DC,
                `xcell_CD4+ T-cells`,`xcell_CD8+ T-cells`,
                xcell_Macrophages,xcell_Monocytes,
                `xcell_NK cells`,xcell_Neutrophils)%>%
  gather(key='cells',value = 'value',-Info_onset)


p = ggplot(df1,aes(x=Info_onset,y=value))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm',formula = y~poly(x,2))+
  facet_wrap(~cells,scales = 'free',ncol = 2)
pdf('Result/A11_xcell.pdf',width = 6,height = 10)
plot(p)
dev.off()

##### plot neutrophil counts ####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  dplyr::select(Lab_pct_neutro_v,xcell_Neutrophils)
p = ggplot(df1,aes(Lab_pct_neutro_v,xcell_Neutrophils))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
pdf('Result/A11_neutrophils.pdf',3,3)
plot(p)
dev.off()
cor.test(df1$Lab_pct_neutro_v,df1$xcell_Neutrophils)





plotDf= df1[,c("Info_onset","Lab_pct_lymph_v")]%>%na.omit()
p = ggplot(plotDf,aes(x=Info_onset,y=Lab_pct_lymph_v))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  theme_bw()
plot(p)

plotDf= df1[,c("Info_onset","Lab_pct_neutro_v")]%>%na.omit()
p = ggplot(plotDf,aes(x=Info_onset,y=Lab_pct_neutro_v))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  theme_bw()
plot(p)

plotDf= df1[,c("Info_onset","xcell_CD4+ T-cells")]%>%na.omit()
p = ggplot(plotDf,aes(x=Info_onset,y=`xcell_CD4+ T-cells`))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,1))+
  theme_bw()
plot(p)

plotDf= df1[,c("Info_onset","xcell_B-cells")]%>%na.omit()
p = ggplot(plotDf,aes(x=Info_onset,y=`xcell_B-cells`))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  theme_bw()
plot(p)

plotDf= df1[,c("Info_onset","GO_B_cell_activation")]%>%na.omit()
p = ggplot(plotDf,aes(x=Info_onset,y=GO_B_cell_activation))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  theme_bw()
plot(p)

plotDf= df1[,c("Lab_pct_neutro_v","xcell_Neutrophils")]%>%na.omit()
p = ggplot(plotDf,aes(x=Lab_pct_neutro_v,y=xcell_Neutrophils))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,1))+
  theme_bw()
plot(p)




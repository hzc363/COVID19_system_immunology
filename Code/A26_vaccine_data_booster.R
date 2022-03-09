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

##### load data #####
vaccine1 = read.csv('Data/Olink_Pfizer vaccine_Prabhu_BP Lab.csv')
vaccine1 = vaccine1%>%select(-QC_Warning)
dfTractory = data.frame()
days = data.frame(Day=0:28)
for (i in 6:ncol(vaccine1)) { 
  t1 = vaccine1[,c('Day',colnames(vaccine1)[i])]%>%
    na.omit()%>%filter(Day<=21)
  colnames(t1)[2]='y'
  LM1 = lm(y~poly(Day,2),t1)
  LM2 = lm(y~1,t1)
  pValue1 = anova(LM1,LM2)$`Pr(>F)`[2]
  pred1 = predict(object=LM1,newdata=days%>%filter(Day<=21))
  
  t1 = vaccine1[,c('Day',colnames(vaccine1)[i])]%>%
    na.omit()%>%filter(Day>=21)
  colnames(t1)[2]='y'
  LM1 = lm(y~poly(Day,2),t1)
  LM2 = lm(y~1,t1)
  pValue2 = anova(LM1,LM2)$`Pr(>F)`[2]
  pred2 = predict(object=LM1,newdata=days%>%filter(Day>=21))
  
  predAll = c(pred1[-length(pred1)],
              (pred1[length(pred1)]+pred2[1])/2,pred2[-1])
  t1 = data.frame(gene = colnames(vaccine1)[i], 
                  Day= days$Day,
                  predValue = predAll,
                  p = min(pValue1,pValue2))
  dfTractory = rbind(dfTractory,t1)
}

dfTractory = dfTractory%>%
  mutate(gene = gsub('olink_|GO_|EVE_','',gene))%>%
  mutate(gene = gsub('_Q.*|_P.*|_O.*','',gene))%>%
  mutate(gene = gsub('_|\\-| ','',gene))%>%
  mutate(gene = toupper(gene))

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)

sigVaccine = pValue$gene%>%unique()

dfPlot = dfTractory%>%
  filter(gene %in% pValue$gene)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(day=paste0('day',Day ))%>%
  dplyr::select(-p,-Day )%>%
  spread(key = day,value=predValue)

dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

geneHC = hclust(dist(dfPlot),method="ward.D2")
geneCL = cutree(geneHC,4)%>%
  factor(labels = c('red','green','blue','orange'))%>%
  as.character()

pdf('Result/A25_heatmap_onset.pdf',width = 10,height = 10)
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          RowSideColors= geneCL[geneHC$order],
          col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'row',hclustfun=f1)
dev.off()
#####infection trajectory #####
infection1 = read.csv('Result/A11_trajectory.csv')
infection1 = infection1%>%
  filter(grepl('olink',gene))%>%
  mutate(gene = gsub('olink_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  mutate(gene = gsub('_|\\-| ','',gene))%>%
  mutate(gene = toupper(gene))

sigInfect= infection1%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)
sigInfect = sigInfect$gene%>%unique()

intersectGenes = intersect(dfTractory$gene,infection1$gene)
sigInfect = sigInfect[sigInfect%in% intersectGenes]
sigVaccine = sigVaccine[sigVaccine %in% intersectGenes]
intersectSig = intersect(sigInfect,sigVaccine)

##### plot scatter ####
t1 = dfTractory%>%select(gene,Day,predValue)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(Day=paste0('Vac',Day))%>%
  group_by(gene)%>%
  mutate(predValue=scale(predValue))
         
t2 = infection1%>%select(gene,Day=Info_onset,predValue)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(Day=paste0('Infect',Day))%>%
  group_by(gene)%>%
  mutate(predValue=scale(predValue))
         
combined = rbind(t1,t2)%>%
  filter(gene%in%union(sigInfect,sigVaccine))

dfPlot = combined%>%
  group_by(Day,gene)%>%
  summarise(predValue=mean(predValue))%>%
  spread(value = predValue,key=Day)%>%as.data.frame()

dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

geneHC = hclust(dist(dfPlot),method="ward.D2")
#geneCL = cutree(geneHC,4)%>%
#  factor(labels = c('red','green','blue','orange'))%>%
#  as.character()



pdf('Result/A24_heatmap_combined.pdf',width = 10,height = 6)
sig1 = 'white' 
sig1[rownames(dfPlot)%in%sigInfect]='black'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)

sig1 = 'white' 
sig1[rownames(dfPlot)%in%sigVaccine]='blue'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)
dev.off()

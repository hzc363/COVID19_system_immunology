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

##### get correlation with day shift #####
dayShift = data.frame()
for (i in 0:30) {
  shiftedInfection = infection1%>%mutate(Info_onset=Info_onset+i)
  combined = inner_join(dfTractory,shiftedInfection,
                        by=c('gene','Day'='Info_onset'))
  combined = combined%>%group_by(gene)%>%
    mutate(predValue.x=scale(predValue.x),
           predValue.y=scale(predValue.y))%>%
    #filter(p.y<0.05)%>%
    filter(gene%in%sigInfect | gene%in%sigVaccine)
  t1= data.frame(
    cor1 = cor(combined$predValue.x,combined$predValue.y),
    days = length(unique(combined$Day)),
    shift = i,ngene = length(unique(combined$gene))
  )
  dayShift = rbind(dayShift,t1)
}
dayShift = dayShift%>%filter(days>10)
p = ggplot(dayShift,aes(y=cor1,x=shift))+geom_bar(stat = 'identity')
plot(p)

intersectGenes = intersect(dfTractory$gene,infection1$gene)
sigInfect = sigInfect[sigInfect%in% intersectGenes]
sigVaccine = sigVaccine[sigVaccine %in% intersectGenes]
intersectSig = intersect(sigInfect,sigVaccine)

##### plot scatter ####
shiftedInfection = infection1%>%mutate(Info_onset=Info_onset+17)
combined = inner_join(dfTractory,shiftedInfection,
                      by=c('gene','Day'='Info_onset'))
combined = combined%>%group_by(gene)%>%
  mutate(predValue.x=scale(predValue.x),
         predValue.y=scale(predValue.y))
dfPlot =combined%>%
  select(gene,Day,vaccineVar = predValue.x,
         infectionVar = predValue.y)%>%
  filter(gene%in%union(sigInfect,sigVaccine))
p = ggplot(dfPlot,aes(x=infectionVar,y=vaccineVar,color=gene)) +
  geom_point()
plot(p)

##### plot heatmap #####
dfPlot =combined%>%
  select(gene,Day,vaccineVar = predValue.x,
         infectionVar = predValue.y)%>%
  filter(gene%in%union(sigInfect,sigVaccine))

dfPlot = dfPlot%>%
  gather(key = 'type',value = 'pred',vaccineVar:infectionVar)
dfPlot = dfPlot%>%mutate(Day=paste0(type,'_Day',Day))%>%
  group_by(gene,Day,type)%>%summarise(pred=mean(pred))
dfPlot = dfPlot%>%select(-type)%>%spread(value = pred,key=Day)%>%as.data.frame()

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

pdf('Result/A24_heatmap_combined.pdf',width = 10,height = 10)
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          #RowSideColors= geneCL[geneHC$order],
          col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)
dev.off()

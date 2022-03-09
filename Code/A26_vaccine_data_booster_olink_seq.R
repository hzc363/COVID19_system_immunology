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
vaccine1 = fread('Result/A24_RNA_seq_go.csv',data.table = F)
dfTractory = data.frame()
#days = data.frame(Info_Day=0:28)
days = data.frame(Info_Day=c(0:7,21:28))

for (i in grep('GO_|olink',colnames(vaccine1))) { 
  t1 = vaccine1[,c('Info_Day',colnames(vaccine1)[i])]%>%
    na.omit()
  colnames(t1)[2]='y'
  #t1$y = scale(t1$y)
  LM1 = lm(y~factor(Info_Day),t1)
  LM2 = lm(y~1,t1)
  pValue1 = anova(LM1,LM2)$`Pr(>F)`[2]
  pred1 = approx(x=t1$Info_Day, y=t1$y,xout=days$Info_Day,method = "linear")
  
  t1 = data.frame(gene = colnames(vaccine1)[i], 
                  Day= days$Info_Day,
                  predValue = pred1$y,
                  p = pValue1,
                  mean1 = mean(t1$y),
                  var1 = var(t1$y))
  dfTractory = rbind(dfTractory,t1)
}

dfTractory = dfTractory%>%
#  mutate(gene = gsub('olink_|GO_|EVE_','',gene))%>%
  mutate(gene = gsub('_Q[0-9].*|_P[0-9].*|_O[0-9].*','',gene))%>%
  mutate(gene = gsub('\\-','',gene))%>%
  mutate(gene = toupper(gene))

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(type= gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)%>%
  filter(grepl('OLINK',gene))

sigVaccine = pValue$gene%>%unique()

#####infection trajectory #####
infection1 = read.csv('Result/A11_trajectory.csv')
infection1 = infection1%>%
 #filter(grepl('olink',gene))%>%
 #mutate(gene = gsub('olink_','',gene))%>%
 mutate(gene = gsub('_P[0-9].*|_O[0-9].*|_Q[0-9].*','',gene))%>%
 mutate(gene = gsub('\\-|','',gene))%>%
 mutate(gene = toupper(gene))

sigInfect= infection1%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)%>%
  filter(grepl('OLINK',gene))

sigInfect = sigInfect$gene%>%unique()

intersectGenes = intersect(dfTractory$gene,infection1$gene)
sigInfect = sigInfect[sigInfect%in% intersectGenes]
sigVaccine = sigVaccine[sigVaccine %in% intersectGenes]
intersectSig = intersect(sigInfect,sigVaccine)

##### plot ####
t1 = dfTractory%>%
  #mutate(predValue=(predValue-mean1)/var1)%>%
  select(gene,Day,predValue)%>%
  mutate(boost = Day<21)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(Day=paste0('Vac',Day))%>%
  group_by(gene)%>%
  mutate(predValue=scale(predValue))%>%
  ungroup()%>%select(-boost)

t2 = infection1%>%
  #mutate(predValue=(predValue-mean1)/var1)%>%
  select(gene,Day=Info_onset,predValue)%>%
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
dfPlot = as.matrix(dfPlot)#%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

col1=grepl('Infec',colnames(dfPlot))|
           (as.integer(gsub('Vac|Infect','',colnames(dfPlot)))<21)
geneHC = hclust(dist(dfPlot),
                method="ward.D2")
#geneCL = cutree(geneHC,4)%>%
#  factor(labels = c('red','green','blue','orange'))%>%
#  as.character()

pdf('Result/A24_heatmap_combined.pdf',width = 10,height = 6)
sig1 = 'white' 
sig1[rownames(dfPlot)%in%sigInfect]='black'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          labRow=gsub('OLINK_','',rownames(dfPlot)[geneHC$order]),
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)

sig1 = 'white' 
sig1[rownames(dfPlot)%in%sigVaccine]='blue'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          labRow=gsub('OLINK_','',rownames(dfPlot)[geneHC$order]),
        
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)
dev.off()

##### plot t cell correlated genes #####
geneDf = read.csv('Result/A25_suplement_table.csv')
geneDf = geneDf%>%
  filter(grepl('olink_',variable))%>%
  mutate(variable = gsub('_P[0-9].*|_O[0-9].*|_Q[0-9].*','',variable))%>%
  mutate(variable = gsub('\\-|','',variable))%>%
  mutate(variable = gsub(' ','',variable))%>%
  mutate(variable = toupper(variable))%>%
  filter(variable%in%intersectGenes)%>%
  arrange(desc(tcell_t))

dfPlot = dfTractory%>%
  select(gene,Day,predValue)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(Day=paste0('Vac',Day))%>%
  group_by(gene)%>%
  mutate(predValue=scale(predValue))%>%
  ungroup()%>%
  filter(gene%in%geneDf$variable[1:10])%>%
  group_by(Day,gene)%>%
  summarise(predValue=mean(predValue))%>%
  spread(value = predValue,key=Day)%>%as.data.frame()

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  filter(grepl('OLINK',gene))%>%
  filter(gene%in%dfPlot$gene)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)

dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)#%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

col1=grepl('Infec',colnames(dfPlot))|
  (as.integer(gsub('Vac|Infect','',colnames(dfPlot)))<21)
geneHC = hclust(dist(dfPlot),
                method="ward.D2")

pdf('Result/A24_heatmap_t_cell_associate.pdf',width = 6,height = 4)
sig1 = 'white' 
sig1[rownames(dfPlot)%in%pValue$gene]='black'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          labRow=gsub('OLINK_','',rownames(dfPlot)[geneHC$order]),
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)
dev.off()



##### plot Ab correlated genes #####
geneDf = read.csv('Result/A25_suplement_table.csv')
geneDf = geneDf%>%
  filter(grepl('olink_',variable))%>%
  mutate(variable = gsub('_P[0-9].*|_O[0-9].*|_Q[0-9].*','',variable))%>%
  mutate(variable = gsub('\\-|','',variable))%>%
  mutate(variable = gsub(' ','',variable))%>%
  mutate(variable = toupper(variable))%>%
  filter(variable%in%intersectGenes)%>%
  arrange(desc(Ab_28days_t))

dfPlot = dfTractory%>%
  select(gene,Day,predValue)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(Day=paste0('Vac',Day))%>%
  group_by(gene)%>%
  mutate(predValue=scale(predValue))%>%
  ungroup()%>%
  filter(gene%in%geneDf$variable[1:10])%>%
  group_by(Day,gene)%>%
  summarise(predValue=mean(predValue))%>%
  spread(value = predValue,key=Day)%>%as.data.frame()

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  filter(grepl('OLINK',gene))%>%
  filter(gene%in%dfPlot$gene)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)

dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)#%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

col1=grepl('Infec',colnames(dfPlot))|
  (as.integer(gsub('Vac|Infect','',colnames(dfPlot)))<21)
geneHC = hclust(dist(dfPlot),
                method="ward.D2")

pdf('Result/A24_heatmap_Ab_associate.pdf',width = 6,height = 4)
sig1 = 'white' 
sig1[rownames(dfPlot)%in%pValue$gene]='black'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          labRow=gsub('OLINK_','',rownames(dfPlot)[geneHC$order]),
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)
dev.off()



##### plot severity correlated genes #####
geneDf = read.csv('Result/A25_suplement_table.csv')
geneDf = geneDf%>%
  filter(grepl('olink_',variable))%>%
  mutate(variable = gsub('_P[0-9].*|_O[0-9].*|_Q[0-9].*','',variable))%>%
  mutate(variable = gsub('\\-|','',variable))%>%
  mutate(variable = gsub(' ','',variable))%>%
  mutate(variable = toupper(variable))%>%
  filter(variable%in%intersectGenes)%>%
  arrange(desc(severity_t))

dfPlot = dfTractory%>%
  select(gene,Day,predValue)%>%
  mutate(Day = str_pad(as.character(Day ), 2, pad = "0",side = 'left'))%>%
  mutate(Day=paste0('Vac',Day))%>%
  group_by(gene)%>%
  mutate(predValue=scale(predValue))%>%
  ungroup()%>%
  filter(gene%in%geneDf$variable[1:10])%>%
  group_by(Day,gene)%>%
  summarise(predValue=mean(predValue))%>%
  spread(value = predValue,key=Day)%>%as.data.frame()

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  filter(grepl('OLINK',gene))%>%
  filter(gene%in%dfPlot$gene)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  filter(pAdj<0.05)

dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)#%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()

col1=grepl('Infec',colnames(dfPlot))|
  (as.integer(gsub('Vac|Infect','',colnames(dfPlot)))<21)
geneHC = hclust(dist(dfPlot),
                method="ward.D2")

pdf('Result/A24_heatmap_severity_associate.pdf',width = 6,height = 4)
sig1 = 'white' 
sig1[rownames(dfPlot)%in%pValue$gene]='black'
heatmap.2(dfPlot[geneHC$order,],Rowv=F, Colv = F,
          labRow=gsub('OLINK_','',rownames(dfPlot)[geneHC$order]),
          RowSideColors= sig1[geneHC$order],col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'none',hclustfun=f1)
dev.off()

##### severity vs ab ######
geneDf = read.csv('Result/A25_suplement_table.csv')
geneDf = geneDf%>%
  mutate(variable = toupper(variable))%>%
  filter(grepl('OLINK_',variable))%>%
  mutate(variable = gsub('_P[0-9].*|_O[0-9].*|_Q[0-9].*','',variable))%>%
  mutate(variable = gsub('\\-|','',variable))%>%
  mutate(variable = gsub(' ','',variable))%>%
  filter(variable%in%intersectGenes)%>%
  mutate(vaccineSig = variable%in% sigVaccine)%>%
  mutate(variable = gsub('OLINK_','',variable))

t1 = geneDf%>%
  filter(vaccineSig==TRUE)%>%
  filter((Ab_28days_t<(-2.5)&tcell_t<(-2.5)))#|(Ab_28days_t>3&tcell_t>3)

p = ggplot(geneDf,aes(x=tcell_t,y=Ab_28days_t,color = vaccineSig))+
  geom_point()+
  ggrepel::geom_label_repel(data = t1,aes(label=variable),size =2)+
  scale_color_manual(values=c("grey", "red"))+theme_bw()

pdf('Result/A26_negative_TB.pdf',width = 5, height = 4)  
plot(p)
dev.off()
plot(geneDf[,c('severity_t','Ab_28days_t','tcell_t')])

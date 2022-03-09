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

##### prepare data #####
RNAseq = cbind('sample_id'=paste0('S',1:nrow(RNAseq)),RNAseq)
expr1 = RNAseq[,-c(1:3)]%>%as.matrix()%>%t()%>%round()
expr1 <- vst(expr1)
#expr1 = apply(expr1, 2, dense_rank)
#expr1 = apply(expr1, 2, function(x)(x/max(x)))
RNAseq[,-c(1:3)] = t(expr1)

##### join with info #####
info1=RNAseq[,1:3]
info1 = info1%>%mutate(id = as.integer(id))%>%
  inner_join(demo1,by='id')%>%
  inner_join(timeTo1[,-2],by='id')%>%
  mutate(day = as.factor(day))

ab1_wide = ab1%>%spread(key = day,value = antibody_value)
ab1_wide = ab1_wide%>%
  mutate(ab0 = (`0`<0.1))%>%
  mutate(ab0 = factor(ab0,labels = c('Ab0_pos','Ab0_neg')))%>%
  dplyr::select(id,ab0)
info1 = info1%>%left_join(ab1_wide,by='id')%>%filter(!is.na(ab0))
info1 = info1%>%left_join(sympOnset,by='id')
info1 = info1%>%left_join(symp1[,c('id','severitycat')],by='id')

RNAseq = info1%>%dplyr::select(-id,-day)%>%
  inner_join(RNAseq,by=c('sample_id'))

##### PCA of samples ####
expr1 = RNAseq[,-c(1:18)]
t1 = apply(expr1,2, var)%>%order(decreasing = T)
pca1 = prcomp(expr1[,t1[1:500]])$x
pca1 = cbind.data.frame(RNAseq[,1:18],pca1[,1:2])
pca1 = pca1 %>%
  #mutate(group = paste(ab0,day))
  mutate(group = paste('day',day))

pdf('Result/A02_PCA_days.pdf',width = 5,height = 4)
library("viridis")
t1 = pca1%>%filter(!is.na(onset))%>%
  mutate(day=(day-onset))
t1$day[t1$day<0]=0;t1$day[t1$day>15]=15
mid1 = median(t1$day)
p = ggplot(t1,aes(x=PC1,y = PC2, color=day))+
  geom_point()+theme_bw()+
  scale_color_gradient2(midpoint = 7.5, low = "#4374B6", mid = "#FFFFBC",
                        high = "#D92E1D", space = "Lab" )
#scale_color_gradientn(colours = viridis(10))
plot(p)

p = ggplot(t1,aes(x=PC1,y = PC2, color=tx))+
  geom_point()+theme_bw()
plot(p)
dev.off()

pca_center = pca1%>%
  group_by(group)%>%
  summarise(PC1=mean(PC1),PC2=mean(PC2))%>%
  mutate(type = 'Center',sample_id=NA)

pca_df = pca1%>%dplyr::select(sample_id,group, PC1,PC2)%>%
  mutate(type = 'Sample')
pca_df = rbind.data.frame(pca_center,pca_df)


p = ggplot(pca1,aes(x=PC1,y = PC2, color=group))+
  geom_point()
plot(p)

pdf('Result/A02_PCA_Ab0.pdf',width = 5,height = 4)
p = ggplot(pca_df,aes(x=PC1,y = PC2, color=group,size = (type=='Center')))+
  geom_point()+theme_bw()
plot(p)
dev.off()

pdf('Result/A02_PCA_tx.pdf',width = 5,height = 4)
pca_center = pca1%>%
  filter(day==5)%>%
  group_by(tx)%>%
  summarise(PC1=mean(PC1),PC2=mean(PC2))%>%
  mutate(type = 'Center',sample_id=NA)

pca_df = pca1%>%filter(day==5)%>%dplyr::select(sample_id,tx, PC1,PC2)%>%
  mutate(type = 'Sample')
pca_df = rbind.data.frame(pca_center,pca_df)

p = ggplot(pca_df,aes(x=PC1,y = PC2, color=tx,size = (type=='Center')))+
  geom_point()+theme_bw()
plot(p)
dev.off()

pdf('Result/A02_PCA_onset.pdf',width = 5,height = 4)
p = ggplot(pca1,aes(x=PC1,y = (day-onset)))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)
dev.off()

##### Correlation with PCs #####
subDf = pca1[,c('sex','age','bmi','race',
                'tx','ab0',"severitycat",'PC1','PC2','onset','day')]%>%
  filter(severitycat!='Asymptomatic')%>%
  mutate(onsetTime = day-onset)%>%
  dplyr::select(-day,-onset)

resultDf = expand_grid(c('PC1','PC2'),
                       colnames(subDf)[!grepl('PC',colnames(subDf))])
colnames(resultDf)=c('PC','Clin')

resultDf$p = NA
resultDf$r2 = NA
for (i in 1:nrow(resultDf)) {
  iDf = data.frame(x = subDf[,resultDf$Clin[i]],
                   y = subDf[,resultDf$PC[i]])
  LM = lm(y~x,iDf)
  sumLM = summary(LM)
  resultDf$r2[i] = sumLM$r.squared
  aovLM = summary.aov(LM)
  resultDf$p[i]=aovLM[[1]]['x','Pr(>F)']
}
plotDf = resultDf%>%dplyr::select(-p)%>%
  spread(key=Clin,value = r2)%>%as.data.frame()
plotDf = data.frame(plotDf[,-1],row.names =plotDf[,1] )

pdf('Result/A02_cor_with_PC.pdf',5,5)
corrplot::corrplot(as.matrix(plotDf)%>%t(),is.corr = F,col = 'blue',
                   hclust.method = 'ward.D')
dev.off()

p = ggplot(subDf,aes(x = onsetTime, y = PC1))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A02_onset_vs_PC1.pdf',width=5,height=5)
plot(p)
dev.off()

##### get PC1 interactions#####
subDf = pca1[,c('sex','age','bmi','race',
                'tx','ab0',"severitycat",'PC1','PC2','onset','day')]%>%
  filter(severitycat!='Asymptomatic')%>%
  mutate(onsetTime = day-onset)%>%
  dplyr::select(-day,-onset)

resultDf = expand_grid(c('PC1'),
                       colnames(subDf)[!grepl('PC|onset',colnames(subDf))])
colnames(resultDf)=c('PC','Clin')

resultDf$p = NA
resultDf$r2 = NA
for (i in 1:nrow(resultDf)) {
  iDf = data.frame(x = subDf[,resultDf$Clin[i]],
                   y = subDf[,resultDf$PC[i]],
                   onset = subDf[,'onsetTime'])%>%na.omit()
  LM1 = lm(y~x+onset,iDf)
  LM2 = lm(y~onset,iDf)
  sumLM1 = summary(LM1)
  sumLM2 = summary(LM2)
  resultDf$r2[i] =sumLM1$r.squared-sumLM2$r.squared
  aovLM = summary.aov(LM1)
  resultDf$p[i]=aovLM[[1]]['x','Pr(>F)']
}
plotDf = resultDf%>%dplyr::select(-p)%>%
  spread(key=Clin,value = r2)%>%as.data.frame()
plotDf = data.frame(plotDf[,-1],row.names =plotDf[,1] )

pdf('Result/A02_interact_with_PC.pdf',5,5)
corrplot::corrplot(as.matrix(plotDf)%>%t(),is.corr = F,col = 'blue',
                   hclust.method = 'ward.D')
dev.off()

p = ggplot(subDf,aes(x = onsetTime, y = PC1,color = factor(sex)))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A02_onset_vs_PC1_sex.pdf',width=5,height=5)
plot(p)
dev.off()

subDf = subDf%>%mutate(PC1_res = lm(PC1~onsetTime)$residuals)
p = ggplot(subDf,aes(x = severitycat, y = PC1_res))+
  geom_boxplot()+geom_jitter(width=0.2)+theme_bw()
pdf('Result/A02_onset_vs_PC1_res_sex.pdf',width=5,height=5)
plot(p)
dev.off()

summary(lm(severitycat~PC1_res+sex,subDf))
t1=table(subDf$sex,subDf$severitycat)
fisher.test(t1)



##### plot onset day #####
plotDf = RNAseq%>%filter(severitycat!='Asymptomatic')%>%
  mutate(dayFactor = as.factor(day))%>%
  mutate(time2onset = day-onset)%>%
  group_by(dayFactor,time2onset)%>%
  mutate(N = 1:n())

pdf('Result/A02_onset_count.pdf',width = 6, height = 5)
p = ggplot(plotDf ,aes(x = time2onset, y = N,color =dayFactor)) +            
  theme_bw()+
  geom_point(position = position_dodge(width = 0.7))
plot(p)
dev.off()

plotDf = ab1%>%
  mutate(id = as.character(id))%>%
  right_join(plotDf,by=c('id','day')) 

p = ggplot(plotDf ,aes(x = time2onset, y = antibody_value)) +            
  theme_bw()+
  geom_point(position = position_dodge(width = 0.7))+
  geom_smooth(method = 'lm')
pdf('Result/A02_ab_time.pdf',width = 3, height = 3)
plot(p)
dev.off()


p = ggplot(plotDf ,aes(x = (day-onset), fill =dayFactor)) +            
  geom_histogram(alpha = 0.5, position = "identity",bins=20)+
  theme_bw()+
  geom_point(data = plotDf%>%filter(day==0),
             aes(x = (day-onset),y = -1,color = dayFactor))+
  geom_point(data = plotDf%>%filter(day==5),
             aes(x = (day-onset),y = -2,color = dayFactor))
pdf('Result/A02_adj_day.pdf',width = 6, height = 5)
plot(p)
dev.off()


plotDf0 = plotDf%>%filter(day==0)
p = ggplot(plotDf0,
           aes(x=ab0,y = day-onset,color = ab0))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  theme_bw()

pdf('Result/A02_ab0_vs_onset.pdf',width = 4,height = 3)
plot(p)
t.test((day-onset)~ab0,plotDf0)
dev.off()

p = ggplot(ab1%>%filter(day==0),aes(x = antibody_value)) +            
  geom_histogram(bins=29)+theme_bw()

plot(p)


#### get adjusted day #####
pca_df = pca_df%>%filter(type == 'Sample')
pca_df$day = 0
pca_df$day[pca_df$group%in%c("Ab0_neg 5","Ab0_pos 0")] = 5
pca_df$day[pca_df$group%in%c("Ab0_pos 5")] = 10

LM1 = lm(day ~ PC1+PC2, pca_df)
pred1 = predict(LM1,pca_df)
pca_df = pca_df%>%mutate(adj_day = pred1)

adj_day = pca_df%>%dplyr::select(sample_id,adj_day)
adj_day = pca1%>%dplyr::select(sample_id,id,day)%>%
  inner_join(adj_day,by='sample_id')
write.csv(adj_day,'Result/A02_adj_day.csv',row.names = F)

adj_day = adj_day%>%mutate(id = as.integer(id))%>%
  inner_join(info1[,c('id','ab0')],by='id')%>%
  unique()%>%
  group_by(ab0,id)%>%
  mutate(adj_day=adj_day-day)%>%
  summarise(adj_day=mean(adj_day))
  
p = ggplot(adj_day, aes(x = adj_day, fill =ab0)) +            
  geom_histogram(alpha = 0.5, position = "identity",bins=20)+
  theme_bw()+
  geom_point(data = adj_day%>%filter(ab0=='Ab0_pos'),aes(x = adj_day,y = -1,color = ab0))+
  geom_point(data = adj_day%>%filter(ab0=='Ab0_neg'),aes(x = adj_day,y = -2,color = ab0))
pdf('Result/A02_adj_day.pdf',width = 6, height = 5)
plot(p)
dev.off()

##### plot symptom onset day #####
sympOnset = fread('Data/ClinicalData/Lambda_symptom_onset.csv',
                  data.table = F)
sympOnset = sympOnset%>%
  inner_join(adj_day,by=c('participantId'='id'))%>%
  mutate(id=participantId)%>%
  dplyr::select(id,ab0,adj_day,symptom_onset)%>%
  mutate(symptom_onset=0-symptom_onset)
sympOnset = sympOnset%>%
  full_join(unique(info1[,c('id','time2prime','time2res')]),by='id')%>%
  gather(key='type',value = 'unadjusted',symptom_onset:time2res)

sympOnset$adjusted = NA

t1 = sympOnset$type=='symptom_onset'
sympOnset$adjusted[t1] = sympOnset$unadjusted[t1]-sympOnset$adj_day[t1]

t1 = sympOnset$type=='time2prime'
sympOnset$adjusted[t1] = sympOnset$unadjusted[t1]+sympOnset$adj_day[t1]

t1 = sympOnset$type=='time2res'
sympOnset$adjusted[t1] = sympOnset$unadjusted[t1]+sympOnset$adj_day[t1]

sympOnset = sympOnset%>%
  gather(key='value_type',value = 'value',unadjusted:adjusted)%>%
  na.omit()

p = ggplot(sympOnset, aes(y = value, x =ab0,color = ab0)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  theme_bw()+
  facet_wrap(value_type~type,scales = 'free')

pdf('Result/A02_adjsted_days.pdf',width =8, height = 4)
plot(p)
dev.off()

t1=sympOnset%>%
  filter(value_type=='adjusted'&type=='symptom_onset')%>%
  dplyr::select(value)%>%unlist()
t1 = c(t1,t1+5)%>%hist()

sympOnset%>%group_by(value_type,type)%>%
  summarise(p = wilcox.test(value~ab0)$p.value)

##### get go score ####
expr1 = RNAseq[,-c(1:18)]
rownames(expr1)=RNAseq$sample_id
expr1 = t(expr1)
expr1 = cbind.data.frame(gene_symbol = rownames(expr1),expr1)

go_df = go_df%>%group_by(go_id)%>%mutate(N=n())%>%
  filter(N>=5)%>%as.data.frame()
go1 = inner_join(go_df,expr1,by="gene_symbol")
go1 = go1%>%dplyr::select(-go_id,-gene_id,-gene_symbol)%>%
  group_by(Term)%>%summarise_all(mean)
go1 = data.frame(go1[,-1],row.names =go1$Term )
go1 = t(go1)
go1 = cbind.data.frame(sample_id = rownames(go1),go1)
colnames(go1)[-1]=paste0('GO_',colnames(go1)[-1])
go1 = go1[-1,]

colnames(RNAseq)[2:18]=paste0('Info_',colnames(RNAseq)[2:18])
colnames(RNAseq)[19:ncol(RNAseq)]=paste0('Gene_',colnames(RNAseq)[19:ncol(RNAseq)])

RNAseq = RNAseq%>%
  inner_join(go1,by='sample_id')

write.csv(RNAseq,'Result/A02_RNA_seq_go.csv',row.names = F)

##### get xcell score #####
xcell = RNAseq[,c(1,grep('^Gene',colnames(RNAseq)))]
xcell = data.frame(xcell[,-1],row.names = xcell[,1])
xcell = t(xcell)
rownames(xcell) = gsub('^Gene_','',rownames(xcell))
xcell = xCellAnalysis(xcell,save.raw = T,
                      file.name='Result/A02_xcell'	)

# use raw score
#xcell = fread('Result/A02_xcell_RAW.txt', data.table = F)
#xcell = data.frame(xcell[,-1],row.names = xcell[,1])

xcell = t(xcell)
xcell = cbind.data.frame(sample_id=rownames(xcell),xcell)
colnames(xcell)[-1]=paste0('xcell_',colnames(xcell)[-1])

RNAseq = RNAseq%>%
  inner_join(xcell,by='sample_id')

write.csv(RNAseq,'Result/A02_RNA_seq_derived.csv',row.names = F)

##### plot interferon #####
plot_df = RNAseq[,grepl('response_to_interferon|^id$|^day|^ab0',
                        colnames(RNAseq))]
plot_df = gather(plot_df,key='var',value = 'value',-id,-day,-ab0)

plot_df = symp1%>%dplyr::select(id,severitycat)%>%
  mutate(id = as.character(id))%>%
  na.omit()%>%inner_join(plot_df,by=c('id'))%>%
  mutate(severitycat = factor(severitycat,levels =c('Asymptomatic',
                                        "Moderate Symptomatic",
                                        "Severe")))%>%
  mutate(group = paste(ab0,day))

p = ggplot(plot_df,aes(color=severitycat,y=value,x = factor(day)))+
  geom_boxplot()+facet_wrap(~var,scales = 'free')+theme_bw()
plot(p)

##### get top 10000 features #####
#expr1 = RNAseq[,-c(1:18)]
#t1 = apply(expr1,2, var)%>%order(decreasing = T)
#expr1 = expr1[,t1[1:10000]]
#RNAseq = cbind(RNAseq[,1:18],expr1)

##### correlate with adj_day #####
df1 = sympOnset%>%
  dplyr::select(id,adj_day)%>%
  unique()%>%
  mutate(id = as.character(id))%>%
  inner_join(RNAseq,by=c('id'))%>%
  mutate(adj_day = adj_day+day)%>%
  filter(!is.na(adj_day))

# remove random effect
dfTractory = data.frame()
days = data.frame(adj_day=round(min(df1$adj_day)):round(max(df1$adj_day)))
for (i in 18:ncol(df1)) { 
  t1 = df1[,c(1,2,i)]%>%na.omit()
  colnames(t1)[3]='y'
  LM1 = lm(y~adj_day+I(adj_day^2),t1)
  LM2 = lm(y~1,t1)
  pValue = anova(LM1,LM2)$`Pr(>F)`[2]
  pred1 = predict(object=LM1,newdata=days)
  t1 = data.frame(gene = colnames(df1)[i], 
                  adj_day= days,
                  predValue = pred1,
                  p = pValue)
  dfTractory = rbind(dfTractory,t1)
}

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(pAdj = p.adjust(p,'BH'))%>%
  filter(pAdj<0.01)

dfPlot = dfTractory%>%
  filter(gene %in% pValue$gene)%>%
  mutate(day=paste0('day',adj_day))%>%
  dplyr::select(-p,-adj_day)%>%
  spread(key = day,value=predValue)
dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D2")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()
heatmap.2(dfPlot,Rowv=T, Colv = F,
          #RowSideColors=t2$col,ColSideColors = col_color,
          col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'row',hclustfun=f1)

plotDf= df1[,c(1,2,grep('CXCR3',colnames(df1)))]%>%na.omit()
p = ggplot(plotDf,aes(x=adj_day,y=CXCR3))+
  geom_point()+geom_line(aes(group=id))+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  theme_bw()
plot(p)
##### correlate with onset #####
df1 = RNAseq%>%
  mutate(onset = day-onset)%>%
  filter(!is.na(onset))%>%
  filter(severitycat!='Asymptomatic')

# remove random effect
dfTractory = data.frame()
days = data.frame(onset=0:15)
for (i in grep('^GO_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c(15,i)]%>%na.omit()
  colnames(t1)[2]='y'
  LM1 = lm(y~poly(onset,2),t1)
  LM2 = lm(y~1,t1)
  pValue = anova(LM1,LM2)$`Pr(>F)`[2]
  pred1 = predict(object=LM1,newdata=days)
  t1 = data.frame(gene = colnames(df1)[i], 
                  onset= days,
                  predValue = pred1,
                  p = pValue)
  dfTractory = rbind(dfTractory,t1)
}

pValue = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(pAdj = p.adjust(p,'BH'))%>%
  filter(pAdj<0.01)

dfPlot = dfTractory%>%
  filter(gene %in% pValue$gene)%>%
  mutate(onset = str_pad(as.character(onset), 2, pad = "0",side = 'left'))%>%
  mutate(day=paste0('day',onset))%>%
  dplyr::select(-p,-onset)%>%
  spread(key = day,value=predValue)
dfPlot = data.frame(dfPlot[,-1],row.names = dfPlot[,1])
col1 = colnames(dfPlot)
dfPlot = as.matrix(dfPlot)%>%apply(1, scale)%>%t()
colnames(dfPlot)=col1

f1 <- function(x) {hclust(x, method="ward.D2")}
t1=colorRampPalette(brewer.pal(7, "RdYlBu"))(25)%>%rev()
pdf('Result/A02_heatmap_GO_onset.pdf',width = 10,height = 10)
heatmap.2(dfPlot,Rowv=T, Colv = F,
          #RowSideColors=t2$col,ColSideColors = col_color,
          col= t1,
          density.info="none", trace="none", dendrogram=c("none"), 
          scale = 'row',hclustfun=f1)
dev.off()
plotDf= df1[,c(16,15,grep('CXCR3',colnames(df1)))]%>%na.omit()
p = ggplot(plotDf,aes(x=onset,y=CXCR3))+
  geom_point()+geom_line(aes(group=id))+
  geom_smooth(method = 'lm',formula = y~poly(x,2))+
  theme_bw()
plot(p)

##### correlate with same day Ab  #####
# get correlation
df1 = ab1%>%mutate(id = as.character(id))%>%
  inner_join(RNAseq,by=c('id','day'))

df_result = data.frame()
for (i in 18:ncol(df1)) {
  t1 = df1[,c(2,3,17,i)]
  colnames(t1)[4]='y'
  t1 = lm(antibody_value~y+day+ab0,t1)
  t1 = summary(t1)
  t1 = data.frame(gene = colnames(df1)[i], 
                  t_ab= t1$coefficients['y','t value'],
                  p_ab = t1$coefficients['y','Pr(>|t|)'])
  df_result = rbind(df_result,t1)
}
df_result = df_result%>%arrange(p_ab)%>%
  mutate(GO = grepl('GO_',gene))%>%group_by(GO)%>%
  mutate(p_ab_adj = p.adjust(p_ab,method = 'fdr'))

##### correlate with M4 ab #####
df1 = ab1%>%filter(day == 120)%>%dplyr::select(-day)%>%
  mutate(id = as.character(id))%>%
  na.omit()%>%inner_join(RNAseq,by=c('id'))

df_result$t_M4 = NA
df_result$p_M4 = NA

for (i in 1:nrow(df_result)) {
  t1 = df1[,c('antibody_value','day','ab0',df_result$gene[i])]
  colnames(t1)[4]='y'
  t1 = lm(antibody_value~y+day+ab0,t1)
  t1 = summary(t1)
  df_result$t_M4[i] = t1$coefficients['y','t value']
  df_result$p_M4[i] = t1$coefficients['y','Pr(>|t|)']
}
df_result = df_result%>%arrange(p_M4)%>%
  mutate(GO = grepl('GO_',gene))%>%group_by(GO)%>%
  mutate(p_M4_adj = p.adjust(p_M4,method = 'fdr'))


##### correlate with symptom #####
df1 = symp1%>%dplyr::select(id,severitycat)%>%
  mutate(id = as.character(id))%>%
  na.omit()%>%inner_join(RNAseq,by=c('id'))%>%
  mutate(severitycat = factor(severitycat,
                              levels =c('Asymptomatic',
                                        "Moderate Symptomatic",
                                        "Severe")))%>%
  mutate(severitycat = as.integer(severitycat))

df_result$t_symp = NA
df_result$p_symp = NA

for (i in 1:nrow(df_result)) {
  t1 = df1[,c('severitycat','day','ab0',df_result$gene[i])]
  colnames(t1)[4]='y'
  t1 = lm(y~severitycat+day+ab0,t1)
  t1 = summary(t1)
  df_result$t_symp[i] = t1$coefficients['severitycat','t value']
  df_result$p_symp[i] = t1$coefficients['severitycat','Pr(>|t|)']
}
df_result = df_result%>%arrange(p_symp)%>%
  mutate(GO = grepl('GO_',gene))%>%group_by(GO)%>%
  mutate(p_symp_adj = p.adjust(p_symp,method = 'fdr'))

##### define category #####
df_result$cat_ab = 'non-sig'
df_result$cat_ab[df_result$p_ab_adj<0.05&df_result$t_ab>0] = 'sig+'
df_result$cat_ab[df_result$p_ab_adj<0.05&df_result$t_ab<0] = 'sig-'

df_result$cat_M4 = 'non-sig'
df_result$cat_M4[df_result$p_M4_adj<0.05&df_result$t_M4>0] = 'sig+'
df_result$cat_M4[df_result$p_M4_adj<0.05&df_result$t_M4<0] = 'sig-'

df_result$cat_symp = 'non-sig'
df_result$cat_symp[df_result$p_symp_adj<0.05&df_result$t_symp>0] = 'sig+'
df_result$cat_symp[df_result$p_symp_adj<0.05&df_result$t_symp<0] = 'sig-'

##### PCA of genes #####
df_sig = df_result%>%
  filter(!(cat_ab=='non-sig'&cat_M4=='non-sig'&cat_symp=='non-sig'))

expr1 = RNAseq[,colnames(RNAseq)%in%df_sig$gene]%>%
  apply(2,rank)%>%t()

pca_fit = prcomp(expr1)
pca1= cbind.data.frame(gene = rownames(expr1),pca_fit$x[,1:2])
colnames(pca1)=c('gene','u1','u2')

df_sig = df_sig%>%inner_join(pca1,by='gene')
p = ggplot(pca1,aes(x=u1,y = u2,color=!grepl("_",gene)))+
  geom_point()+theme_bw()#+ylim(-10,10)+xlim(-15,15)
plot(p)

##### cluster ######
# cluster
expr1 = RNAseq[,colnames(RNAseq)%in%df_sig$gene]%>%
  apply(2,rank)%>%t()
D = dist(expr1)
HC = hclust(D,method = 'ward.D')

# find cluter N
sil = c()
for (i in 2:10) {
  CT = cutree(HC,k = i)
  ss <- silhouette(CT, D)
  sil = c(sil,mean(ss[, 3]))
}
optimal_N = which.max(sil)+1

# get cluster info
t1 = cutree(HC,k = 4)
t1 = data.frame(cluster = t1,gene =names(t1) )
df_sig = left_join(df_sig,t1,by='gene')

##### plot #####
t1 = RNAseq[,-c(1:16)]
t1 = t1[,!colnames(t1)%in%df_sig$gene]%>%
  apply(2,rank)
pred1 = predict(pca_fit,t(t1))
t1 = data.frame(gene = colnames(t1),
                cat_ab='non-sig',
                cat_M4='non-sig',
                cat_symp='non-sig',
                u1=pred1[,1],u2=pred1[,2],
                cluster=0)

plot_df = df_sig%>%
  dplyr::select(gene,cat_ab,cat_M4,cat_symp,u1,u2,cluster)%>%
  rbind(t1)%>%
  gather(key = 'outcome',value = 'cat',cat_ab,cat_M4,cat_symp,cluster)#%>%
#filter(cat!='non-sig')

plot_df$alpha = 1
plot_df$alpha[plot_df$cat%in%c('non-sig','0')]=0.01
plot_df = plot_df%>%
  mutate(cat=factor(cat,levels=c('sig+','sig-','non-sig',
                                 '1','2','3','4','0')))%>%
  arrange(desc(cat))

p = ggplot(plot_df,aes(x=u1,y = u2,color=cat,alpha=alpha))+
  geom_point()+theme_bw()+facet_wrap(~outcome,nrow=2)
plot(p)


##### enrichR #####
sub_df = plot_df%>%filter(cat!='non-sig'&outcome!='cluster')%>%
  mutate(group = paste(outcome,cat))
f1 = function(x){
  t1 = sub_df$gene[sub_df$group==x]
  t1 = enrichr(t1, 'GO_Biological_Process_2018')[[1]]
  return(cbind(group = x,t1))
  }
sub_list = lapply(unique(sub_df$group), f1)
sub_df = do.call(rbind,sub_list)
sub_df = sub_df%>%filter(Adjusted.P.value<0.05)

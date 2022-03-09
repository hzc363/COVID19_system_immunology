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
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
RNAseq = RNAseq[,grepl('olink|Info',colnames(RNAseq))]
RNAseq=RNAseq%>%mutate(Info_id = as.character(Info_id))

colnames(olink)[-(1:2)]=paste0('olink_',colnames(olink)[-(1:2)])
RNAseq = inner_join(RNAseq,olink,by=c('Info_id'='id','Info_day'='day'))

##### join with clinical data #####
neut_ab = read.csv('Data/Antibodies/D28_M7 Neut IgG Binding.csv')
df_combine = neut_ab%>%
  dplyr::select(Info_id=ID,D28_IgG=d28_IgG_AUC,D28_Neut,
                M7_IgG = M7_IgG_AUC,M7_Neut)%>%
  mutate(D28_IgG = log2(D28_IgG+1),
         D28_Neut=log2(D28_Neut+1),
         M7_IgG=log2(M7_IgG+1) ,
         M7_Neut=log2(M7_Neut+1))%>%
  mutate(Info_id = as.character(Info_id))
RNAseq = right_join(df_combine,RNAseq,by='Info_id')

tcell = read.csv('Data/Tcell_Data/Lambda ICS D28_combination cytokines_bgsub_freq of parent.csv')
RNAseq = tcell%>%
  filter(Stim%in%c('MN','SS1'))%>%
  dplyr::select(Info_id=id,mCD4)%>%
  group_by(Info_id)%>%summarise(mCD4=sum(mCD4))%>%
  right_join(RNAseq,by='Info_id')

RNAseq = opAUC%>%mutate(id = as.character(id))%>%
  right_join(RNAseq,by=c('id'='Info_id'))

colnames(RNAseq) = gsub('\\-| ','',colnames(RNAseq))
txt1='olink|mCD4|^auc|M7_IgG|D28_IgG|Info_age|Info_sex|Info_severitycat'
RNAseq = RNAseq[,grepl(txt1,colnames(RNAseq))]

##### binaries y #####
RNAseq$Info_severitycat = factor(grepl('Sever',RNAseq$Info_severitycat),
                                 labels = c('low','high'))
RNAseq$auc = factor(RNAseq$auc>median(RNAseq$auc,na.rm = T),
                    labels = c('low','high'))
RNAseq$mCD4 = factor(RNAseq$mCD4>median(RNAseq$mCD4,na.rm = T),
                     labels = c('low','high'))
RNAseq$D28_IgG = factor(RNAseq$D28_IgG>median(RNAseq$D28_IgG,na.rm = T),
                        labels = c('low','high'))
RNAseq$M7_IgG = factor(RNAseq$M7_IgG>median(RNAseq$M7_IgG,na.rm = T),
                       labels = c('low','high'))

RNAseq = RNAseq%>%dplyr::select(-olink_IL6_P05231_OID00482)
RNAseq = as.data.frame(RNAseq)

##### get 30 with highest var #####
t1 = colnames(RNAseq)[grepl('olink',colnames(RNAseq))]
xlist = c('Info_sex','Info_age',t1)
yList = c("auc","mCD4","D28_IgG","M7_IgG","Info_severitycat")
comboY = RNAseq%>%na.omit()
t1 = apply(comboY[,yList],1,paste0,collapse='_')
comboY = cbind(comboY,y = t1)

dfTrain = comboY[,c('y',xlist)]
RF = randomForest(factor(y)~.,dfTrain)
t1 = data.frame(var = rownames(RF$importance))
t1$gini = RF$importance[,1]
imp1 = t1 %>%mutate(gini = gini/mean(gini))%>%
  arrange(desc(gini))

##### cross validation #####
predDf = expand.grid(1:15,1:nrow(RNAseq),yList)
colnames(predDf)=c("N","leaveOut",'clin')
predDf$y = NA
predDf$pred = NA
predDf$clin = as.character(predDf$clin)
for (i in 1:nrow(predDf)) {
  dfAll = RNAseq[,c(predDf$clin[i],xlist)]
  colnames(dfAll)[1]="y"
  dfAll = na.omit(dfAll)
  dfTrain = dfAll[-predDf$leaveOut[i],]
  set.seed(1)
  high1 = (dfTrain$y=='high')
  s1 = sample(which(high1),sum(!high1),replace = T)
  s1 = c(which(!high1),s1)
  dfTrain = dfTrain[s1,]
  impSub = imp1#%>%filter(y==predDf$clin[i])
  v1 = unique(c('y',impSub$var[1:predDf$N[i]]))
  RF = randomForest(y~.,dfTrain[,v1])
  predDf$pred[i]=predict(RF,dfAll[predDf$leaveOut[i],-1,drop=F],
                         type='prob')[,2]
  predDf$y[i] = as.character(dfAll[predDf$leaveOut[i],1])
}

auc1 = predDf%>%na.omit()%>%group_by(N,clin)%>%
  summarise(auc=pROC::roc(y,pred)$auc[[1]])

p = ggplot(auc1,aes(N, auc,color = clin,group=clin))+
  geom_point()+theme_bw()+geom_line()+ylim(0.5,0.9)
pdf('Result/A20_auc_vs_N.pdf',width = 5,height=4)
plot(p)
dev.off()

save(auc1, file = 'Result/A20_auc1.rda')

##### final models #####
load('Result/A20_auc1.rda')
finalImp = data.frame()
RFlist = list()
auc10 = auc1%>%filter(N<=10)%>%
  arrange(desc(auc))%>%group_by(clin)%>%
  summarise(N = N[1],auc = max(auc))%>%
  mutate(data='olink+Demo')
yList = auc10$clin
N = auc10$N

for (i in 1:length(yList)) {
  impSub = imp1%>%filter(y==yList[i])
  v1 = unique(c(yList[i],impSub$var[1:N[i]]))
  dfTrain = RNAseq[,v1]
  colnames(dfTrain)[1]="y"
  dfTrain = na.omit(dfTrain)
  
  set.seed(1)
  high1 = (dfTrain$y=='high')
  s1 = sample(which(high1),sum(!high1),replace = T)
  s1 = c(which(!high1),s1)
  dfTrain = dfTrain[s1,]
  
  RF = randomForest(factor(y)~.,dfTrain)
  t1 = data.frame(var = rownames(RF$importance))
  t1$gini = RF$importance[,1]
  t1$y = yList[i]
  finalImp = finalImp%>%rbind(t1)
  RFlist[[yList[i]]]=RF
}

save(RFlist, file = 'Result/A20_final_models.rda')

##### plot feature importance #####
plotDf = finalImp%>%spread(key = y,value = gini,fill = 0)%>%
  mutate(var = gsub('Info_|olink_','',var))%>%
  mutate(var = gsub('_.*','',var))
plotDf = data.frame(plotDf[,-1],row.names = plotDf[,1])
plotDf = plotDf[,c('Info_severitycat','mCD4','D28_IgG','M7_IgG')]
plotDf = plotDf[apply(plotDf, 1, sum)>0,]
r1 = hclust(dist(plotDf))$order
c1 = hclust(dist(t(plotDf)))$order
p1 = as.matrix(plotDf[r1,c1]==0)
pdf('Result/A20_importance.pdf',10,10)
t1 = as.matrix(plotDf[r1,c1])
corrplot::corrplot(t1,is.corr = F,
                   col = 'blue',tl.srt = 45,p.mat = p1,
                   insig = "pch")
dev.off()

##### age and gender #####
predDf = expand.grid(1:nrow(RNAseq),yList)
colnames(predDf)=c("leaveOut",'clin')
predDf$y = NA
predDf$pred = NA
predDf$clin = as.character(predDf$clin)
for (i in 1:nrow(predDf)) {
  dfAll = RNAseq[,c(predDf$clin[i],c('Info_age','Info_sex'))]
  colnames(dfAll)[1]="y"
  dfAll = na.omit(dfAll)
  dfTrain = dfAll[-predDf$leaveOut[i],]
  set.seed(1)
  high1 = (dfTrain$y=='high')
  s1 = sample(which(high1),sum(!high1),replace = T)
  s1 = c(which(!high1),s1)
  dfTrain = dfTrain[s1,]
  RF = randomForest(y~.,dfTrain)
  predDf$pred[i]=predict(RF,dfAll[predDf$leaveOut[i],-1,drop=F],
                         type='prob')[,2]
  predDf$y[i] = as.character(dfAll[predDf$leaveOut[i],1])
}

aucDemo = predDf%>%na.omit()%>%group_by(clin)%>%
  summarise(auc=roc(y,pred)$auc[[1]])%>%
  mutate(data = 'Demo')

aucAll = rbind(aucDemo,auc10%>%dplyr::select(-N))
p = ggplot(aucAll,aes(x = clin, y = auc, fill = data))+
  geom_bar(stat = 'identity',position = 'dodge')+
  coord_cartesian(ylim = c(0.5, 0.85)) +
  theme_bw()
pdf('Result/A20_auc_bar.pdf',width = 5,height = 4)
plot(p)
dev.off()

##### validation #####
olinkValid = read.csv('Data/Olink_validation.csv',check.names = F)
olinkValid = olinkValid%>%
  filter(Group!=0)%>%
  mutate(y = Group>1)%>%
  mutate(Info_age = Age)%>%
  mutate(Info_sex = (Gender=='F')*1)
olinkValid = olinkValid%>%na.omit()
colnames(olinkValid) = gsub('\\-| ','',colnames(olinkValid))

pred = predict(RFlist$Info_severitycat,olinkValid,type='prob')[,1]
roc1=pROC::roc(olinkValid$y,pred)
pr1 <- pr.curve(pred[olinkValid$y], pred[!olinkValid$y],curve = TRUE );
print(roc1)
plot(roc1)

p = data.frame(prob = pred, y = olinkValid$y)
p = ggplot(p,aes(x=y,y = pred))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+theme_bw()
pdf('Result/A20_test_box.pdf',3,3)
plot(p)
dev.off()

##### validation 2 #####
olinkValid = read.csv('Data/OLINK Proteomics/Olink_Science_BP Lab.csv',
                      check.names = F)
nameDf1 = finalImp%>%filter(y=='Info_severitycat')
nameDf1 = nameDf1%>%
  mutate(prot = gsub('olink_','',var))%>%
  mutate(prot = gsub('_.*','',prot))

nameDf2 = data.frame(varValid = colnames(olinkValid))%>%
  mutate(prot = gsub('\\-| ','',varValid))

nameDf = left_join(nameDf1,nameDf2,by='prot')

olinkValid = olinkValid%>%
  filter(Group!=0)%>%
  mutate(y = Group>1)%>%
  mutate(Info_age = Age)%>%
  mutate(Info_sex = (Gender=='F')*1)
olinkValid = olinkValid%>%na.omit()
colnames(olinkValid) = gsub('\\-| ','',colnames(olinkValid))

pred = predict(RFlist$Info_severitycat,olinkValid,type='prob')[,1]
roc1=pROC::roc(olinkValid$y,pred)
pr1 <- pr.curve(pred[olinkValid$y], pred[!olinkValid$y],curve = TRUE );
print(roc1)
plot(roc1)

p = data.frame(prob = pred, y = olinkValid$y)
p = ggplot(p,aes(x=y,y = pred))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+theme_bw()
pdf('Result/A20_test_box.pdf',3,3)
plot(p)
dev.off()

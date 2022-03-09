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

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
#t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
#       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
#geneDf = read.csv('A27_top_go_genes.csv')
#t1 = paste0('Gene_',geneDf$gene_symbol)
#RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))|colnames(RNAseq)%in%t1]
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

colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### correlate with severity #####
df1 = RNAseq%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                              levels =c('Asymptomatic',
                                        "Moderate Symptomatic",
                                        "Severe")))%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  mutate(Info_severitycat = as.integer(Info_severitycat))

dfResult = data.frame()
for (i in grep('^GO_|^IgG_|^Lab_|^olink_|^xcell_|Gene_|EVE_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c(15,16,i)]%>%na.omit()
  colnames(t1)[3]='y'
  LM1 = lm(y~poly(Info_onset,2)+Info_severitycat,t1)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['Info_severitycat','Estimate'], 
                  t= LM1$coefficients['Info_severitycat','t value'],
                  p = LM1$coefficients['Info_severitycat','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}

write.csv(dfResult,'Result/A12_severity_results.csv',row.names = F)

dfResult = dfResult%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)%>%
  filter(!grepl('TBA',gene))%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))

protCL=read.csv('Result/A11_gene_CL.csv')
dfPlot = dfResult%>%filter(round(pAdj,2)<=0.05)
dfPlot = dfPlot%>%
  mutate(gene = gsub('olink_|GO_|EVE_','',gene))%>%
  mutate(gene = gsub('_Q.*|_P.*|_O.*','',gene))%>%
  mutate(gene = gsub('_',' ',gene))%>%
  left_join(protCL,by=c('gene'='prot'))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})

p = ggplot(dfPlot,
           aes(x=t,y=reorder(gene,t)))+
  geom_point(aes(x=6,color = CL),shape = 15,size = 5)+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A12_symp_bar.pdf',width = 7,height = 7)
plot(p)
dev.off()

##### rig_I vs symptom #####
p = ggplot(RNAseq,aes(x=Info_severitycat,y=olink_DDX58_O95786_OID01018))+
  geom_boxplot()+geom_jitter(width = 0.1,height = 0)+theme_bw()
pdf('Result/A15_DDX58_vs_severity.pdf',4,4)
plot(p)
dev.off()

##### plot trajectory #####
df1 = RNAseq%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                                   levels =c('Asymptomatic',
                                             "Moderate Symptomatic",
                                             "Severe")))%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')

t1 = colnames(df1)[grepl('IP\\-10|olink_MCP\\-1|PRDX5|CXCL11|MCP\\-2',
                         colnames(df1))]
df1 = df1[,c('Info_onset','Info_severitycat',t1)]%>%
  gather(key='prot',value = 'value',-Info_onset,-Info_severitycat)%>%
  filter(Info_onset>0 & Info_onset<20)%>%
  mutate(Info_onset_cat = factor(ceiling(Info_onset/7),
                                 label = c('week 1','week 2','week 3')))


pdf('Result/A12_symp_timeline.pdf',width = 4.5,height = 8)
p = ggplot(df1,aes(x=Info_onset,y=value,
                   color = grepl('M',Info_severitycat)))+
  geom_point(size=1)+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,1),se=F)+
  facet_wrap(~prot,scales = 'free',ncol = 1)
plot(p)

p = ggplot(df1,aes(x=Info_onset,y=value))+
  geom_point()+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,2),se=F)+
  facet_wrap(~prot,scales = 'free',ncol = 1)
plot(p)
dev.off()

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


pdf('Result/A12_symp_timeline_lab.pdf',width = 7,height = 2.5)
p = ggplot(df1,aes(x=Info_onset,y=value,color = grepl('M',Info_severitycat)))+
  geom_point(size=1)+theme_bw()+
  geom_smooth(method = 'lm',formula = y~poly(x,1),se=F)+
  facet_wrap(~prot,scales = 'free',ncol = 2)
plot(p)
dev.off()

##### plot boxplot #####
t1 = colnames(RNAseq)[grepl('IP\\-10|olink_MCP\\-1|PRDX5|CXCL11|MCP\\-2',
                         colnames(RNAseq))]
df1 = RNAseq[,c('Info_severitycat',t1)]%>%
  gather(key='prot',value = 'value',-Info_severitycat)

pdf('Result/A12_symp_boxplot.pdf',width = 4.5,height = 8)
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

##### prediction glm #####
dfResultSub = dfResult%>%
  filter(grepl('olink',gene))%>%
  arrange(p)

df1 = RNAseq[,grep('Info_age|Info_sex|^olink_|severity',colnames(RNAseq))]
df1 = na.omit(df1)
df1 = df1%>%
  dplyr::rename(y=Info_severitycat)%>%
  mutate(y=grepl('Severe',y))

# get features that's not correlated with age and gender
dfCor = data.frame()
for (i in grep('^olink',colnames(df1))) {
  t1 = lm(df1[,i]~df1[,'Info_age']+df1[,'Info_sex'])
  t1 = predict(t1,df1[,c('Info_age','Info_sex')])
  t1 = cor(t1,df1[,i])
  t1 = data.frame(var =colnames(df1)[i], cor = t1)
  dfCor = rbind(dfCor, t1)
}

# get 20 with highest var
t1 = apply(df1[,-(1:3)], 2, var)%>%sort(decreasing = T)
df1 = df1[,c('y',"Info_age",'Info_sex',names(t1)[1:20])]

# find lambdas
set.seed(1)
s1 = sample(which(df1$y==TRUE),length(df1$y)-sum(df1$y),replace = T)
s1 = c(which(df1$y==FALSE),s1)
t1 = df1[s1,]
l1 = glmnet(x=as.matrix(t1[,-1]),y = t1$y)

# get the optimal lambda
predDf = expand.grid(l1$lambda,1:nrow(df1))
colnames(predDf)=c("lambda","leaveOut")
predDf$y = df1$y[predDf$leaveOut]
predDf$pred = NA
predDf$N = NA
for (i in 1:nrow(predDf)) {
  t1 = df1[-predDf$leaveOut[i],]
  set.seed(1)
  s1 = sample(which(t1$y==TRUE),length(t1$y)-sum(t1$y),replace = T)
  s1 = c(which(t1$y==FALSE),s1)
  t1 = t1[s1,]
  t1 = glmnet(x=as.matrix(t1[,-1]),y = t1$y,family = 'binomial',
              lambda = predDf$lambda[i])
  predDf$pred[i]=predict(t1,newx=as.matrix(df1[predDf$leaveOut[i],-1,drop=F]))
  predDf$N[i] = sum(t1$beta[,1]!=0)
}
auc1 = predDf%>%group_by(lambda)%>%
  summarise(auc=pROC::roc(y,pred)$auc[[1]],
            N = mean(N))

p = ggplot(auc1,aes(N, auc))+
  geom_point()+theme_bw()
pdf('Result/A12_auc_vs_N.pdf',width = 3,height = 3)
plot(p)
dev.off()

# fit final model
set.seed(1)
s1 = sample(which(df1$y==TRUE),length(df1$y)-sum(df1$y),replace = T)
s1 = c(which(df1$y==FALSE),s1)
t1 = df1[s1,]
l1 = glmnet(x=as.matrix(t1[,-1]),y = t1$y,
            lambda = 0.015)#auc1$lambda[which.max(auc1$auc)])
var1 = rownames(l1$beta)[l1$beta[,1]!=0]
print(length(var1))

# fit glm
l1 = glm(y~.,t1[,c('y',var1)],family = 'binomial')
l2 = glm(y~.,t1[,c('y','Info_age','Info_sex')],family = 'binomial')

plotDf = data.frame(gene = gsub('`','',rownames(summary(l1)$coefficients)[-1]),
                    z = summary(l1)$coefficients[-1,'z value'],
                    p = summary(l1)$coefficients[-1,'Pr(>|z|)'])
p = ggplot(plotDf,aes(y=reorder(gene,z), x = z))+
  geom_bar(stat = 'identity')+
  theme_bw()
pdf('Result/A12_BarPlot.pdf',5,5)
plot(p)
dev.off()

##### training error #####
pred = predict(l1,df1,type='response')
roc1=pROC::roc(df1$y,pred)
pr1 <- pr.curve(pred[df1$y], pred[!df1$y],curve = TRUE );
print(roc1)

pred = predict(l2,df1,type='response')
roc2=pROC::roc(df1$y,pred)
pr2 <- pr.curve(pred[df1$y], pred[!df1$y],curve = TRUE );
print(roc2)


pdf('Result/A12_ROC_train.pdf',3,3)
plot(roc1)
plot(roc2)
dev.off()


pdf('Result/A12_PR_curve_train.pdf',height=3,width=3.5)
plot(pr1)
plot(pr2)
dev.off()

##### validation #####
olinkValid = read.csv('Data/Olink_validation.csv',check.names = F)
olinkValid = olinkValid%>%
  filter(Group!=0)%>%
  mutate(y = Group>1)%>%
  mutate(Info_age = Age)%>%
  mutate(Info_sex = (Gender=='F')*1)
olinkValid = olinkValid%>%na.omit()

pred = predict(l1,olinkValid,type='response')
roc1=pROC::roc(olinkValid$y,pred)
pr1 <- pr.curve(pred[olinkValid$y], pred[!olinkValid$y],curve = TRUE );
print(roc1)

pred = predict(l2,olinkValid,type='response')
roc2=pROC::roc(olinkValid$y,pred)
pr2 <- pr.curve(pred[olinkValid$y], pred[!olinkValid$y],curve = TRUE );
print(roc2)


pdf('Result/A12_ROC.pdf',3,3)
plot(roc1)
plot(roc2)
dev.off()


pdf('Result/A12_PR_curve.pdf',height=3,width=3.5)
plot(pr1)
plot(pr2)
dev.off()



##### random forest prepare #####

dfResultSub = dfResult%>%
  filter(grepl('olink',gene))%>%
  arrange(p)

df1 = RNAseq[,grep('Info_age|Info_sex|^olink_|severity',colnames(RNAseq))]
df1 = na.omit(df1)
df1 = df1%>%
  dplyr::rename(y=Info_severitycat)%>%
  mutate(y=grepl('Severe',y))

# get 30 with highest var
t1 = apply(df1[,-(1:3)], 2, var)%>%sort(decreasing = T)
df1 = df1[,c('y',"Info_age",'Info_sex',names(t1)[1:30])]
df1 = data.frame(df1)

set.seed(1)
s1 = sample(which(df1$y==TRUE),length(df1$y)-sum(df1$y),replace = T)
s1 = c(which(df1$y==FALSE),s1)
t1 = data.frame(df1[s1,])

l1 = randomForest(factor(y)~.,t1)

imp1 = data.frame(var = rownames(l1$importance),
                  gini = l1$importance[,1])%>%
  arrange(desc(gini))

imp1 = data.frame(var = rownames(l1$importance),
                  gini = l1$importance[,1])%>%
  arrange(desc(gini))
p= ggplot(imp1,aes(y=reorder(var,gini),x=gini))+
  geom_bar(stat = 'identity')+
  theme_bw()
plot(p)

##### leave one out RF #####
predDf = expand.grid(1:nrow(imp1),1:nrow(df1))
colnames(predDf)=c("N","leaveOut")
predDf$y = df1$y[predDf$leaveOut]
predDf$pred = NA
for (i in 1:nrow(predDf)) {
  t1 = df1[-predDf$leaveOut[i],]
  set.seed(1)
  s1 = sample(which(t1$y==TRUE),length(t1$y)-sum(t1$y),replace = T)
  s1 = c(which(t1$y==FALSE),s1)
  t1 = t1[s1,]
  v1 = unique(c('y',imp1$var[1:predDf$N[i]]))
  t1 = randomForest(factor(y)~.,t1[,v1])
  predDf$pred[i]=predict(t1,df1[predDf$leaveOut[i],-1,drop=F],type='prob')[,2]
}

auc1 = predDf%>%group_by(N)%>%
  summarise(auc=pROC::roc(y,pred)$auc[[1]])

p = ggplot(auc1,aes(N, auc))+
  geom_point()+theme_bw()

pdf('Result/A12_auc_vs_N_RF.pdf',width = 3,height = 3)
plot(p)
dev.off()

set.seed(1)
s1 = sample(which(df1$y==TRUE),length(df1$y)-sum(df1$y),replace = T)
s1 = c(which(df1$y==FALSE),s1)
t1 = data.frame(df1[s1,])

v1 = unique(c('y',imp1$var[1:8]))
l1 = randomForest(factor(y)~.,t1[,v1])
l2 = randomForest(factor(y)~.,t1[,c('y','Info_age','Info_sex')])

imp2 = data.frame(var = rownames(l1$importance),
                  gini = l1$importance[,1])%>%
  arrange(desc(gini))%>%
  mutate(var = gsub('olink_','',var))%>%
  mutate(var = gsub('_.*','',var))

p= ggplot(imp2,aes(y=reorder(var,gini),x=gini))+
  geom_bar(stat = 'identity')+
  theme_bw()
pdf('Result/A12_bar_RF.pdf',width = 2.5,height = 4)
plot(p)
dev.off()
##### training error #####
pred = predict(l1,df1,type='prob')[,1]
roc1=pROC::roc(df1$y,pred)
pr1 <- pr.curve(pred[df1$y], pred[!df1$y],curve = TRUE );
print(roc1)

pred = predict(l2,df1,type='prob')[,1]
roc2=pROC::roc(df1$y,pred)
pr2 <- pr.curve(pred[df1$y], pred[!df1$y],curve = TRUE );
print(roc2)


pdf('Result/A12_ROC_train_RF.pdf',3,3)
plot(roc1)
plot(roc2)
dev.off()


pdf('Result/A12_PR_curve_train_RF.pdf',height=3,width=3.5)
plot(pr1)
plot(pr2)
dev.off()

##### validation RF #####
olinkValid = read.csv('Data/Olink_validation.csv',check.names = F)
olinkValid = olinkValid%>%
  filter(Group!=0)%>%
  mutate(y = Group>1)%>%
  mutate(Info_age = Age)%>%
  mutate(Info_sex = (Gender=='F')*1)
olinkValid = olinkValid%>%na.omit()%>%data.frame()

pred = predict(l1,olinkValid,type='prob')[,1]
roc1=pROC::roc(olinkValid$y,pred)
pr1 <- pr.curve(pred[olinkValid$y], pred[!olinkValid$y],curve = TRUE );
print(roc1)

pred = predict(l2,olinkValid,type='prob')[,1]
roc2=pROC::roc(olinkValid$y,pred)
pr2 <- pr.curve(pred[olinkValid$y], pred[!olinkValid$y],curve = TRUE );
print(roc2)


pdf('Result/A12_ROC_RF.pdf',3,3)
plot(roc1)
plot(roc2)
dev.off()


pdf('Result/A12_PR_curve_RF.pdf',height=3,width=3.5)
plot(pr1)
plot(pr2)
dev.off()

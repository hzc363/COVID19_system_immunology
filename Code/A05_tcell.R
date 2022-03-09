library(data.table)
library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library("pheatmap")
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(survival)
require(GO.db)
require(biomaRt)
library(org.Hs.eg.db)
library(fastcluster)
library(cluster)
library(survminer)
library(gpairs)

##### load data #####
# ab1,demo1,eve1,olink,RNAseq,timeTo1, tcell
load('Result/A01_compiled_data.rda')

##### get data #####
tcell = cbind('samples'=1:nrow(tcell),tcell)
info1=tcell[,1:3]
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
info1 = info1%>%left_join(symp1,by='id')
  
expr1 = info1%>%dplyr::select(samples)%>%inner_join(tcell,by='samples')
expr1 =  t(as.matrix(expr1[,-(1:3)]))
colnames(expr1)=info1$samples

##### PCA for samples#####
t1 = apply(expr1, 2, function(x){sum(!is.na(x))})
expr1_t = expr1[,colnames(expr1)!='120']
info1_t = info1[info1$samples!='120',]
D1 = 1-cor(expr1_t,use = 'pairwise.complete.obs')#,method = 'spearman')
pca1 = cmdscale(as.dist(D1),eig=TRUE, k=2)$points[,1:2]
colnames(pca1) = c('PC1','PC2')
pca1 = cbind.data.frame(info1_t,pca1)
pca1 = pca1 %>%mutate(group = paste(ab0,day))

pca_center = pca1%>%
  group_by(group)%>%
  summarise(PC1=mean(PC1),PC2=mean(PC2))%>%
  mutate(type = 'Center')
pca_df = pca1%>%dplyr::select(group, PC1,PC2)%>%
  mutate(type = 'Sample')
pca_df = rbind.data.frame(pca_center,pca_df)


p = ggplot(pca1,aes(x=PC1,y = PC2))+
  geom_point()
plot(p)

p = ggplot(pca1,aes(x=PC1,y = PC2, color=severitycat))+
  geom_point()
plot(p)

p = ggplot(pca1,aes(x=PC1,y = PC2, color=ab0))+
  geom_point()
plot(p)

##### PCA for variables#####
expr1_t = expr1
expr1_t[expr1_t<0]=0
t1 = apply(expr1_t, 1, sd,na.rm=T)
expr1_t = expr1_t[t1>0&rownames(expr1_t)!='SS1_bgsub_Memory CD8p_CD107apIFNypIL10nIL21pTNFn',]
D1 = 1-cor(t(expr1_t),use = 'pairwise.complete.obs',method = 'spearman')
pca1 = cmdscale(as.dist(D1),eig=TRUE, k=2)$points[,1:2]
colnames(pca1) = c('PC1','PC2')
pca1 = as.data.frame(pca1)
pca1$gene = rownames(pca1)
p = ggplot(pca1,aes(x=PC1,y = PC2))+
  geom_point()
plot(p)

##### direct correlate with antibody titer genes #####
expr1 = t(expr1)
expr1 = cbind.data.frame(id = info1$id,day=info1$day,
                         ab0 = info1$ab0,expr1)

df1 = ab1%>%mutate(day = factor(day))%>%
  inner_join(expr1,by=c('id','day'))
df_cor = data.frame()
for (i in 6:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,3],method = 'spearman')
  t1 = data.frame(gene = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_cor = rbind(df_cor,t1)
}
df_cor = df_cor%>%arrange(p)%>%
  mutate(p_adj = p.adjust(p,method = 'fdr'))
write.csv(df_cor,'Result/A05_cor_ab_gene.csv',row.names = F)

plot_df = inner_join(pca1,df_cor,by='gene')

plot_df$cat = 'non-sig'
plot_df$cat[plot_df$p_adj<0.05&plot_df$cor>0] = 'Sig cor>0'
plot_df$cat[plot_df$p_adj<0.05&plot_df$cor<0] = 'Sig cor<0'

p = ggplot(plot_df,aes(x=PC1,y = PC2,color=cat))+
  geom_point()
plot(p)

plot_df = plot_df%>%filter(p_adj<0.05)

# plot
plot(df1[,c('MN_bgsub_CD4p_IL21','SS1_bgsub_CD4p_TNF','antibody_value')])
t1= df1[,c('MN_bgsub_CD4p_IL21','SS1_bgsub_CD4p_TNF','antibody_value')]%>%
  apply( 2, rank)%>%as.data.frame()
plot(t1)

##### correlate with M4 ab #####
df1 = ab1%>%filter(day == 120)%>%na.omit()
t1 = expr1
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day.x, -day.y)

df_M4 = data.frame()
for (i in 5:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,2])
  t1 = data.frame(gene = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_M4 = rbind(df_M4,t1)
}
df_M4 = df_M4%>%arrange(p)%>%mutate(p_adj = p.adjust(p,method = "fdr"))
write.csv(df_M4,'Result/A05_cor_with_M4_Ab.csv',row.names = F)

plot_df = inner_join(pca1,df_M4,by='gene')
plot_df$cat = 'non-sig'
plot_df$cat[plot_df$p<0.01&plot_df$cor>0] = 'Sig cor>0'
plot_df$cat[plot_df$p<0.01&plot_df$cor<0] = 'Sig cor<0'

p = ggplot(plot_df,aes(x=PC1,y = PC2,color=cat))+
  geom_point()
plot(p)

plot_df = plot_df%>%filter(p<0.01)

# plot
plot(df1[,c('MN_bgsub_CD4p_IL21','SS1_bgsub_CD4p_TNF','antibody_value')])
t1= df1[,c('MN_bgsub_CD4p_IL21','SS1_bgsub_CD4p_TNF','antibody_value')]%>%
  apply( 2, rank)%>%as.data.frame()
plot(t1)

##### correlate with symptom severity #####
df1 = symp1%>%dplyr::select(id,severitycat)%>%na.omit()
t1 = expr1
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day)
df1 = cbind(severity_score = 1,df1)
df1$severity_score[df1$severitycat=='Moderate Symptomatic']=2
df1$severity_score[df1$severitycat=='Severe']=3


df_symp = data.frame()
for (i in 5:ncol(df1)) {
  t1 = df1[,c(1:4,i)]
  colnames(t1)[5]="y"
  t1 = lm(y~severity_score+ab0,t1)
  t1 = summary(t1)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  t = t1['severity_score','t value'],
                  p = t1['severity_score','Pr(>|t|)'])
  df_symp = rbind(df_symp,t1)
}
df_symp = df_symp%>%arrange(p)%>%mutate(p_adj = p.adjust(p,method = "fdr"))
write.csv(df_symp,'Result/A05_cor_with_symp.csv',row.names = F)

plot_df = inner_join(pca1,df_symp,by='gene')
plot_df$cat = 'non-sig'
plot_df$cat[plot_df$p_adj<0.05&plot_df$t>0] = 'Sig beta>0'
plot_df$cat[plot_df$p_adj<0.05&plot_df$t<0] = 'Sig beta<0'

p = ggplot(plot_df,aes(x=PC1,y = PC2,color=cat))+
  geom_point()
plot(p)

plot_df = plot_df%>%filter(p_adj<0.05)

##### correlate with symptom cluster #####
df1 = symp1%>%
  dplyr::select(id,symptomcluster)%>%
  na.omit()%>%
  filter(symptomcluster!='')

t1 = expr1
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day)
df1 = cbind(severity_score = 1,df1)
df1$severity_score[df1$symptomcluster=='B']=2

df_symp = data.frame()
for (i in 5:ncol(df1)) {
  t1 = df1[,c(1:4,i)]
  colnames(t1)[5]="y"
  t1 = lm(y~severity_score+ab0,t1)
  t1 = summary(t1)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  t = t1['severity_score','t value'],
                  p = t1['severity_score','Pr(>|t|)'])
  df_symp = rbind(df_symp,t1)
}
df_symp = df_symp%>%arrange(p)%>%mutate(p_adj = p.adjust(p,method = "fdr"))
write.csv(df_symp,'Result/A05_cor_with_symp.csv',row.names = F)

plot_df = inner_join(pca1,df_symp,by='gene')
plot_df$cat = 'non-sig'
plot_df$cat[plot_df$p_adj<0.05&plot_df$t>0] = 'Sig beta>0'
plot_df$cat[plot_df$p_adj<0.05&plot_df$t<0] = 'Sig beta<0'

p = ggplot(plot_df,aes(x=PC1,y = PC2,color=cat))+
  geom_point()
plot(p)

plot_df = plot_df%>%filter(p<0.01)
##### associate with clinical info #####

tcell = tcell%>%
  mutate(id = as.integer(id))
plot_df = info1%>%inner_join(tcell,by='id')%>%
  dplyr::select(sex, age, bmi, race, ethnic, tx, time2prime,ab0,
                time2prog,time2res,severitycat,
                symptomcluster,MN_bgsub_CD4p_IL21,SS1_bgsub_CD4p_TNF)

plot_m = model.matrix(~.,plot_df)
cor1 = cor(plot_m,method = 'spearman',use = 'pairwise.complete.obs')
cor1[upper.tri(cor1, diag = T)] = 100

cor1 = cbind.data.frame(var1=rownames(cor1),cor1)%>%
  gather(key="var2",value='cor',-1)%>%
  filter(cor<1)%>%
  mutate(abs_cor = abs(cor))

ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=age))+
  geom_point()
ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=factor(race)))+
  geom_point()
ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=factor(ethnic)))+
  geom_point()
ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=severitycat))+
  geom_point()
ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=symptomcluster))+
  geom_point()
ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=factor(tx)))+
  geom_point()
ggplot(plot_df,aes(x=MN_bgsub_CD4p_IL21,y=SS1_bgsub_CD4p_TNF,color=factor(ab0)))+
  geom_point()



plot(plot_df$MN_bgsub_CD4p_IL21~plot_df$age)
cor.test(plot_df$MN_bgsub_CD4p_IL21,plot_df$age)

boxplot(plot_df$MN_bgsub_CD4p_IL21~plot_df$race)
summary.aov(lm(MN_bgsub_CD4p_IL21~race,plot_df))

boxplot(plot_df$MN_bgsub_CD4p_IL21~plot_df$ethnic)
boxplot(plot_df$MN_bgsub_CD4p_IL21~plot_df$severitycat)
boxplot(plot_df$MN_bgsub_CD4p_IL21~plot_df$symptomcluster)

plot(plot_df$SS1_bgsub_CD4p_TNF~plot_df$age)
boxplot(plot_df$SS1_bgsub_CD4p_TNF~plot_df$race)
boxplot(plot_df$SS1_bgsub_CD4p_TNF~plot_df$ethnic)
boxplot(plot_df$SS1_bgsub_CD4p_TNF~plot_df$severitycat)
boxplot(plot_df$SS1_bgsub_CD4p_TNF~plot_df$symptomcluster)





##### combined results #####
olink_sig = result_df%>%filter(stage_p_adj<0.05)%>%
  dplyr::select(from=gene,p = stage_p_adj)%>%
  mutate(direction = 0,to='stage',relation='associated with',
         dataset = 'olink',from_type = 'protein',to_type = 'clinical')%>%
  dplyr::select(from, to,relation,direction,p,from_type,to_type,dataset)

t1 = result_df_BP%>%
  mutate(stage_p_adj = p.adjust(stage_p,method = 'fdr'))%>%
  filter(stage_p_adj<0.05)%>%
  dplyr::select(from=gene,p = stage_p_adj)%>%
  mutate(direction = 0,to='stage',relation='associated with',
         dataset = 'olink',from_type = 'protein_GO',to_type = 'clinical')%>%
  dplyr::select(from, to,relation,direction,p,from_type,to_type,dataset)
olink_sig = rbind(olink_sig,t1)

t1 = df_cor%>%
  filter(p_adj<0.05)%>%
  dplyr::select(from=gene,p = p_adj,direction=cor)%>%
  mutate(direction = sign(direction),to='Current Ab',relation='associated with',
         dataset = 'olink',from_type = 'protein',to_type = 'Ab')%>%
  dplyr::select(from, to,relation,direction,p,from_type,to_type,dataset)
olink_sig = rbind(olink_sig,t1)

t1 = df_cor_BP%>%
  filter(p_adj<0.05)%>%
  dplyr::select(from=pathway,p = p_adj,direction=cor)%>%
  mutate(direction = sign(direction),to='Current Ab',relation='associated with',
         dataset = 'olink',from_type = 'protein_GO',to_type = 'Ab')%>%
  dplyr::select(from, to,relation,direction,p,from_type,to_type,dataset)
olink_sig = rbind(olink_sig,t1)

t1 = df_M4%>%
  mutate(p_adj = p.adjust(p,method = 'fdr'))%>%
  filter(p_adj<0.05)%>%
  dplyr::select(from=gene,p = p_adj,direction=cor)%>%
  mutate(direction = sign(direction),to='M4 Ab',relation='associated with',
         dataset = 'olink',from_type = 'protein',to_type = 'Ab')%>%
  dplyr::select(from, to,relation,direction,p,from_type,to_type,dataset)
olink_sig = rbind(olink_sig,t1)







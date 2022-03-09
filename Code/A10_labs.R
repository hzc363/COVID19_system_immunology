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
library(ggrepel)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')

##### prepare data #####
RNAseq = cbind('sample_id'=paste0('S',1:nrow(lab1)),lab1)
t1 = apply(RNAseq, 1, function(x){sum(is.na(x))})
RNAseq = RNAseq[t1<10,]

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
RNAseq = info1%>%dplyr::select(-id,-day)%>%
  inner_join(RNAseq,by=c('sample_id'))%>%
  filter(day%in% c(0,5))

##### PCA of samples ####
#RNAseq = RNAseq%>%filter(sample_id!='S120')
#info1 = info1%>%filter(sample_id!='S120')

expr1 = t(RNAseq[,-c(1:16)])%>%apply(2,as.numeric)
D1 = 1-cor(expr1,use = 'pairwise.complete.obs',method = 'spearman')
pca1 = cmdscale(as.dist(D1),eig=TRUE, k=2)$points[,1:2]
colnames(pca1) = c('PC1','PC2')
pca1 = cbind.data.frame(RNAseq[,1:16],pca1[,1:2])
pca1 = pca1 %>%mutate(group = paste(ab0,day))

pca_center = pca1%>%
  group_by(group)%>%
  summarise(PC1=mean(PC1),PC2=mean(PC2))%>%
  mutate(type = 'Center')
pca_df = pca1%>%dplyr::select(group, PC1,PC2)%>%
  mutate(type = 'Sample')
pca_df = rbind.data.frame(pca_center,pca_df)

p = ggplot(pca1,aes(x=PC1,y = PC2, color=group))+
  geom_point()
plot(p)

p = ggplot(pca_df,aes(x=PC1,y = PC2, color=group,size = (type=='Center')))+
  geom_point()
plot(p)

##### correlate with same day Ab  #####
# get correlation
df1 = ab1%>%
  mutate(day = as.character(day))%>%
  inner_join(RNAseq,by=c('id','day'))

df_result = data.frame()
for (i in 18:ncol(df1)) {
  t1 = df1[,c(2,3,17,i)]
  colnames(t1)[4]='y'
  t1 = lm(antibody_value~y+ab0,t1)
  t1 = summary(t1)
  if(!'y'%in% rownames(t1$coefficients)){next}
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
  na.omit()%>%inner_join(RNAseq,by=c('id'))

df_result$t_M4 = NA
df_result$p_M4 = NA

for (i in 1:nrow(df_result)) {
  t1 = df1[,c('antibody_value','day','ab0',df_result$gene[i])]
  colnames(t1)[4]='y'
  t1 = lm(antibody_value~y+ab0,t1)
  #t1 = lm(antibody_value~y,t1)
  t1 = summary(t1)
  if(!'y'%in% rownames(t1$coefficients)){next}
  df_result$t_M4[i] = t1$coefficients['y','t value']
  df_result$p_M4[i] = t1$coefficients['y','Pr(>|t|)']
}
df_result = df_result%>%arrange(p_M4)%>%
  mutate(GO = grepl('GO_',gene))%>%group_by(GO)%>%
  mutate(p_M4_adj = p.adjust(p_M4,method = 'fdr'))


##### correlate with symptom #####
df1 = symp1%>%dplyr::select(id,severitycat)%>%
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
  t1 = lm(severitycat~y+ab0,t1)
  t1 = summary(t1)
  if(!'y'%in% rownames(t1$coefficients)){next}
  df_result$t_symp[i] = t1$coefficients['y','t value']
  df_result$p_symp[i] = t1$coefficients['y','Pr(>|t|)']
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
df_sig = df_result#%>%
#filter(!(cat_ab=='non-sig'&cat_M4=='non-sig'&cat_symp=='non-sig'))

expr1 = RNAseq[,colnames(RNAseq)%in%df_sig$gene]#%>%
#apply(2,rank)%>%t()
t1 = apply(expr1, 2, sd,na.rm=T)
expr1 = expr1[,t1>0&colnames(expr1)!='SS1_bgsub_Memory CD8p_CD107apIFNypIL10nIL21pTNFn']
D1 = 1-cor(expr1,use = 'pairwise.complete.obs',method = 'spearman')
pca1 = cmdscale(as.dist(D1),eig=TRUE, k=2)$points[,1:2]
colnames(pca1) = c('PC1','PC2')
#pca_fit = prcomp(expr1)
pca1= cbind.data.frame(gene = colnames(expr1),pca1[,1:2])
colnames(pca1)=c('gene','u1','u2')

df_sig = df_sig%>%inner_join(pca1,by='gene')
p = ggplot(pca1,aes(x=u1,y = u2,color=!grepl("_",gene)))+
  geom_point()+theme_bw()#+ylim(-10,10)+xlim(-15,15)
plot(p)

##### plot #####
#t1 = RNAseq[,-c(1:16)]
#t1 = t1[,!colnames(t1)%in%df_sig$gene]%>%
#  apply(2,rank)
#pred1 = predict(pca_fit,t(t1))
#t1 = data.frame(gene = colnames(t1),
#                #cat_ab='non-sig',
#                cat_M4='non-sig',
#                cat_symp='non-sig',
#               u1=pred1[,1],u2=pred1[,2])

plot_df = df_sig%>%as.data.frame()%>%
  dplyr::select(gene,#cat_ab,
                cat_M4,cat_symp,u1,u2)%>%
  #rbind(t1)%>%
  gather(key = 'outcome',value = 'cat',#cat_ab,
         cat_M4,cat_symp)#%>%
#filter(cat!='non-sig')

plot_df$alpha = 1
plot_df$alpha[plot_df$cat%in%c('non-sig','0')]=0.01
plot_df = plot_df%>%
  mutate(cat=factor(cat,levels=c('sig+','sig-','non-sig')))%>%
  arrange(desc(cat))

p = ggplot(plot_df,aes(x=u1,y = u2,color=cat,alpha=alpha))+
  geom_point()+
  #geom_label_repel(data=plot_df%>%filter(cat!='non-sig'),
  #aes(x=u1,y = u2,color=cat, label = gene),size=1)+
  theme_bw()+
  facet_wrap(~outcome,nrow=1)
pdf('Result/A10_PCA.pdf',width = 15, height = 7)
plot(p)
dev.off()

# bar plots
plot_df1 = df_sig%>%as.data.frame()%>%
  dplyr::select(gene,t_M4,t_symp)%>%
  gather(key = 'outcome',value = 't',t_M4,t_symp)

plot_df2 = df_sig%>%as.data.frame()%>%
  dplyr::select(gene,cat_M4,cat_symp)%>%
  gather(key = 'outcome',value = 'cat',cat_M4,cat_symp)

plot_df = plot_df1[plot_df2$cat!='non-sig',]%>%
  arrange(desc(t))

p = ggplot(plot_df,aes(x=t,y = reorder(gene,t),fill=(t<0)))+
  geom_bar(stat = 'identity')+theme_bw()+
  facet_wrap(~outcome,nrow=1,scales = 'free')
plot(p)

##### save #####
write.csv(df_result,'Result/A10_tcell_result.csv',row.names = F)


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
library(org.Hs.eg.db)


##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')

##### check data #####
plot(olink$IL10_P22301_OID00528,olink$IL10_P22301_OID00993)
cor(olink$IL10_P22301_OID00528,olink$IL10_P22301_OID00993)
plot(olink$IL6_P05231_OID00482,olink$IL6_P05231_OID00947)
cor(olink$IL6_P05231_OID00482,olink$IL6_P05231_OID00947)

##### prepare data #####
RNAseq = cbind('sample_id'=paste0('S',1:nrow(olink)),olink)

prot1 = data.frame(id = colnames(RNAseq)[-c(1:3)])
prot1 = prot1%>%
  mutate(uniprot = gsub('_OID.*','',id))%>%
  mutate(uniprot = gsub('.*_','',uniprot))

#uniprot to id
x <- org.Hs.egUNIPROT
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
id2uniprot = lapply(1:length(xx), function(i){
  data.frame(gene_id = names(xx)[i],uniprot = xx[[i]])
})
id2uniprot  = do.call(rbind,id2uniprot)

# id to symbol
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
id2symbol = lapply(1:length(xx), function(i){
  data.frame(gene_id = names(xx)[i], symbol= xx[[i]])
})
id2symbol  = do.call(rbind,id2symbol)

# uniprot to symbol
uniprot2symbol = inner_join(id2uniprot,id2symbol,by='gene_id')%>%
  group_by(uniprot)%>%summarise_all(function(x){x[1]})

prot1 = prot1%>%
  left_join(uniprot2symbol,by='uniprot')
colnames(RNAseq)[-c(1:3)]=prot1$symbol

##### combine duplicated genes #####
t1 = colnames(RNAseq)%>%table()
t1 = names(t1)[t1>1]
m_df = data.frame('a'=rep(NA,nrow(RNAseq)))
col1 = c()
for(c1 in t1){
  i = which(colnames(RNAseq)==c1)
  col1 = c(col1,i)
  m1 = RNAseq[,i]%>%apply(1, mean)
  m_df = cbind(m_df,m1)
}
m_df = m_df[,-1]
colnames(m_df)=t1
RNAseq = RNAseq[,-col1]%>%cbind(m_df)

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
RNAseq = RNAseq%>%filter(sample_id!='S180')
info1 = info1%>%filter(sample_id!='S180')

expr1 = RNAseq[,-c(1:18)]
t1 = apply(expr1,2, var)%>%order(decreasing = T)
pca1 = prcomp(expr1[,t1])$x
pca1 = cbind.data.frame(RNAseq[,1:18],pca1[,1:2])
pca1 = pca1 %>%
  #mutate(group = paste(ab0,day))
  mutate(group = paste("day",day))

# find bridge IDs
t1 = pca1%>%filter(day==0)%>%arrange(sex,PC1)
t1 = t1[(0:7)*13+2,c("sample_id","sex","age"  ,"severitycat",
                     "id", "day", "PC1", "PC2" )]
print(t1%>%arrange(PC1))


pdf('Result/A04_PCA_days.pdf',width = 5,height = 4)
library("viridis")
t1 = pca1%>%filter(!is.na(onset))%>%
  mutate(day=(day-onset))
t1$day[t1$day<0]=0;t1$day[t1$day>15]=15
mid1 = median(t1$day)
p = ggplot(t1,aes(x=PC1,y = PC2, color=day))+
  geom_point()+theme_bw()+
  scale_color_gradient2(midpoint = 7.5, low = "#4374B6", mid = "#FFFFBC",
                        high = "#D92E1D", space = "Lab" )
plot(p)

p = ggplot(t1,aes(x=PC1,y = PC2, color=tx))+
  geom_point()+theme_bw()
plot(p)
dev.off()

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
  geom_point()+theme_bw()
pdf('Result/A04_PCA_Ab0.pdf',width = 5,height = 4)
plot(p)
dev.off()

pdf('Result/A04_PCA_tx.pdf',width = 5,height = 4)
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

pdf('Result/A04_cor_with_PC.pdf',5,5)
corrplot::corrplot(as.matrix(plotDf)%>%t(),is.corr = F,col = 'blue',
                   hclust.method = 'ward.D')
dev.off()

p = ggplot(subDf,aes(x = onsetTime, y = PC2))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
pdf('Result/A04_onset_vs_PC1.pdf',width=5,height=5)
plot(p)
dev.off()

##### get go score ####
# expr1 = RNAseq[,-c(1:16)]
# rownames(expr1)=RNAseq$sample_id
# expr1 = t(expr1)
# expr1 = cbind.data.frame(gene_symbol = rownames(expr1),expr1)
# 
# go_df = go_df%>%
#   mutate(in_olink = gene_symbol%in%expr1$gene_symbol)%>%
#   group_by(go_id)%>%mutate(N=sum(in_olink),Perc = sum(in_olink)/n())%>%
#   filter(N>=2&Perc>=0.2)%>%filter(in_olink==T)%>%as.data.frame()
# 
# go1 = inner_join(go_df,expr1,by="gene_symbol")
# go1 = go1%>%dplyr::select(-go_id,-gene_id,-gene_symbol)%>%
#   group_by(Term)%>%summarise_all(mean)
# go1 = data.frame(go1[,-1],row.names =go1$Term )
# go1 = t(go1)
# go1 = cbind.data.frame(sample_id = rownames(go1),go1)
# colnames(go1)[-1]=paste0('GO_',colnames(go1)[-1])
# 
# RNAseq = RNAseq%>%
#   inner_join(go1,by='sample_id')

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

df1 = symp1%>%dplyr::select(id,severitycat)%>%
  mutate(id = as.character(id))%>%
  na.omit()%>%inner_join(RNAseq,by=c('id'))%>%
  mutate(severitycat = factor(severitycat,
                              levels =c('Asymptomatic',
                                        "Moderate Symptomatic",
                                        "Severe")))
p = ggplot(df1,aes(x = factor(day),y=IFNG,color=severitycat))+
  geom_boxplot()+theme_bw()
plot(p)
##### define category #####
df_result$cat_ab = 'non-sig'
df_result$cat_ab[df_result$p_ab<0.05&df_result$t_ab>0] = 'sig+'
df_result$cat_ab[df_result$p_ab<0.05&df_result$t_ab<0] = 'sig-'

df_result$cat_M4 = 'non-sig'
df_result$cat_M4[df_result$p_M4<0.05&df_result$t_M4>0] = 'sig+'
df_result$cat_M4[df_result$p_M4<0.05&df_result$t_M4<0] = 'sig-'

df_result$cat_symp = 'non-sig'
df_result$cat_symp[df_result$p_symp<0.05&df_result$t_symp>0] = 'sig+'
df_result$cat_symp[df_result$p_symp<0.05&df_result$t_symp<0] = 'sig-'

##### PCA of genes #####
df_sig = df_result%>%
  #filter(!(cat_ab=='non-sig'&cat_M4=='non-sig'&cat_symp=='non-sig'))%>%
  filter(GO==F)

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
#t1 = RNAseq[,-c(1:16)]
# t1 = t1[,!colnames(t1)%in%df_sig$gene]%>%
#   apply(2,rank)
# pred1 = predict(pca_fit,t(t1))
# t1 = data.frame(gene = colnames(t1),
#                 #cat_ab='non-sig',
#                 cat_M4='non-sig',
#                 cat_symp='non-sig',
#                 u1=pred1[,1],u2=pred1[,2])

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
  geom_point()+theme_bw()+facet_wrap(~outcome,nrow=1)
plot(p)

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

##### enrichR #####
sub_df = plot_df%>%
  filter(cat!='non-sig'&cat!='cat_ab'&outcome!='cluster')%>%
  mutate(group = paste(outcome,cat))
f1 = function(x){
  t1 = sub_df$gene[sub_df$group==x]
  t1 = enrichr(t1, 'GO_Biological_Process_2018')[[1]]
  return(cbind(group = x,t1))
}
sub_list = lapply(unique(sub_df$group), f1)
sub_df = do.call(rbind,sub_list)
sub_df = sub_df%>%
  mutate(Overlap = gsub('\\/.*','',Overlap))%>%
  mutate(Overlap = as.integer(Overlap))%>%
  filter(Overlap>2)%>%
  filter(!grepl("^positive reg|^regulation of|^negative reg",Term))%>%
  mutate(Adjusted.P.value = p.adjust(P.value,'fdr'))%>%
  filter(Adjusted.P.value<0.1)

#### save #####
write.csv(df_result,'Result/A04_olink_result.csv',row.names = F)
write.csv(RNAseq,'Result/A04_olink_with_symbol.csv',row.names = F)

  

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



##### load data #####
# ab1,demo1,eve1,olink,RNAseq,timeTo1
load('Result/A01_compiled_data.rda')

##### get deseq2 obj #####
RNAseq = cbind('samples'=1:nrow(RNAseq),RNAseq)
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

expr1 = info1%>%dplyr::select(samples)%>%inner_join(RNAseq,by='samples')
expr1 =  round(t(as.matrix(expr1[,-(1:3)])))
colnames(expr1)=info1$samples

dds <- DESeqDataSetFromMatrix(countData = expr1 ,
                              colData =info1,design = ~ day+tx+ab0)
keep <- rowSums(counts(dds)) >= 30
dds <- dds[keep,]

vsd <- vst(dds, blind=TRUE)

##### PCA #####
pca1= plotPCA(vsd, intgroup=c("ab0", "day"),ntop = 500)
pca_center = pca1$data%>%
  group_by(group)%>%
  summarise(PC1=mean(PC1),PC2=mean(PC2))%>%
  mutate(type = 'Center')
pca_df = pca1$data%>%dplyr::select(group, PC1,PC2)%>%
  mutate(type = 'Sample')
pca_df = rbind.data.frame(pca_center,pca_df)

p = ggplot(pca1$data,aes(x=PC1,y = PC2, color=group))+
  geom_point()
plot(p)

p = ggplot(pca_df,aes(x=PC1,y = PC2, color=group,size = (type=='Center')))+
  geom_point()
plot(p)

##### Clustering #####
vsd_m = assay(vsd)
colnames(vsd_m)= paste(info1$id,info1$day,sep = '_')
# D = 1-cor(t(vsd_m))
# HC = fastcluster::hclust(as.dist(D),method = 'ward.D2')
# sil = c()
# for (i in (1:10)*100) {
#   CT = cutree(HC,k = i)
#   ss <- silhouette(CT, D)
#   sil = c(sil,mean(ss[, 3]))
# }
# optimal_N = which.max(sil)+1
#plot(HC)

#### correlations #####
COR = cor(vsd_m,method = 'spearman')
COR = COR%>%as.data.frame()%>%
  mutate(id1 = rownames(.))%>%
  gather(key='id2',value = 'cor',-id1)

t1 = info1%>%dplyr::select(id,tx)%>%mutate(id = as.character(id))
COR= COR%>%
  filter(id1!=id2)%>%
  mutate(day1=gsub('.*_','',id1))%>%
  mutate(day2=gsub('.*_','',id2))%>%
  mutate(id1=gsub('_.*','',id1))%>%
  mutate(id2=gsub('_.*','',id2))%>%
  inner_join(t1,by=c('id1'='id'))%>%
  inner_join(t1,by=c('id2'='id'))

t1 = read.csv('Result/A03_ab_cluster.csv')
t1$id = as.character(t1$id)
COR=COR%>%inner_join(t1,by=c('id1'='id'))%>%
  inner_join(t1,by=c('id2'='id'))

boxplot(cor~(id1==id2),COR,outline=FALSE)
boxplot(cor~(day1==day2),COR,outline=FALSE)
boxplot(cor~((day1==day2)&(cl.x==cl.y)),COR,outline=FALSE)

t1 = COR%>%filter((day1==5)&(day2==5))
boxplot(cor~(tx.x==tx.y),t1,outline=FALSE)

##### differential expression #####
expr1 = t(vsd_m)
t1 = gsub('_.*','',rownames(expr1))%>%as.integer()
t2 = gsub('.*_','',rownames(expr1))%>%as.integer()
expr1 = cbind.data.frame(id = t1,day=t2,expr1)
expr1 = inner_join(unique(info1[,c('id','ab0')]),expr1,by='id')
expr1 = cbind.data.frame(stage =NA,expr1)
expr1$stage[expr1$day==0&expr1$ab0=='Ab0_neg']='stage1'
expr1$stage[expr1$day==5&expr1$ab0=='Ab0_neg']='stage2'
expr1$stage[expr1$day==0&expr1$ab0=='Ab0_pos']='stage2'
expr1$stage[expr1$day==5&expr1$ab0=='Ab0_pos']='stage3'
expr1 = expr1%>%na.omit()

# result_df = data.frame()
# for (i in 5:ncol(expr1)) {
#   t1 = expr1[,c(1,2,i)]
#   colnames(t1)[3]='Y'
#   LM1 = lmer(Y~stage+(1|id),t1,REML = F)
#   EM1 = emmeans(LM1,'stage')
#   CT1 = contrast(EM1, 'tukey')
#   AN1 = anova(LM1)
#   EM1 = summary(EM1)
#   CT1 = summary(CT1)
#   t1 = data.frame(
#     gene = colnames(expr1)[i],
#     stage_p = AN1$`Pr(>F)`,
#     stage1_mean = EM1[1,'emmean'],
#     stage2_mean = EM1[2,'emmean'],
#     stage3_mean = EM1[3,'emmean'],
#     stage1_se = EM1[1,'SE'],
#     stage2_se = EM1[2,'SE'],
#     stage3_se = EM1[3,'SE'],
#     stage1v2_beta = CT1[1,'estimate'],
#     stage1v3_beta = CT1[2,'estimate'],
#     stage2v3_beta = CT1[3,'estimate'],
#     stage1v2_se = CT1[1,'SE'],
#     stage1v3_se = CT1[2,'SE'],
#     stage2v3_se = CT1[3,'SE'],
#     stage1v2_p = CT1[1,'p.value'],
#     stage1v3_p = CT1[2,'p.value'],
#     stage2v3_p = CT1[3,'p.value']
#   )
#   result_df = rbind(result_df,t1)
# }
# write.csv(result_df,'Result/A02_DE_table.csv',row.names = F)

result_df = read.csv('Result/A02_DE_table.csv')

##### heatmap #####
plot_df = result_df%>%
  arrange(stage_p)%>%
  dplyr::select(gene,stage1_mean,stage2_mean,stage3_mean)
plot_df = plot_df[1:30,]%>%as.data.frame()
plot_df = data.frame(plot_df[,-1],row.names = plot_df[,1])
colnames(plot_df)=gsub('_.*','',colnames(plot_df))

c1 = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
#pdf("Result/C07_lm_heatmap.pdf",width = 5, height = 10)
pheatmap(as.matrix(plot_df),  scale="row",color = c1)
result_df = result_df%>%
  mutate(stage_p_adj = p.adjust(stage_p,method = 'fdr'))

##### go analysis #####

filter1 <- as.list(GOBPANCESTOR)
filter1 = sapply(filter1, function(x){'GO:0002376'%in% x})
filter1 = names(filter1)[filter1]

m_list = lapply(filter1, function(x){
  allegs = mget(x, org.Hs.egGO2ALLEGS, ifnotfound = NA)[[1]]
  if(is.na(allegs[1])){return(NA)}
  return(unlist(mget(allegs,org.Hs.egSYMBOL)))
})
filter1 = data.frame(go_id = filter1)
tbl=toTable(GOTERM)
tbl = tbl[,c(1,3)]%>%mutate(Term=gsub(" ",'_',Term))%>%unique()
tbl = left_join(filter1,tbl, by = 'go_id')
names(m_list)=tbl$Term
t1 = sapply(m_list, function(x){sum(!is.na(x))})
m_list = m_list[t1>=5]
m_list = m_list[!grepl('^regulation_of|negative_regu|positive_regu',names(m_list))]

#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") 
#gene.data <- getBM(attributes=c('hgnc_symbol', 'go_id'),
#                   filters = 'go', values = filter1, mart = ensembl)
#gene.data = gene.data%>%filter(go_id%in%filter1)
#tbl=toTable(GOTERM)
#tbl = tbl[,c(1,3)]%>%mutate(Term=gsub(" ",'_',Term))
#gene.data  = gene.data %>%inner_join(tbl,by='go_id')
#m_list = gene.data %>% split(x = .$hgnc_symbol, f = .$Term)

#m_df = msigdbr(species = "Homo sapiens")
#m_df = m_df%>%filter(gs_subcat=="BP")
#m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

expr_m = expr1[,-(1:4)]%>%t()
gsea1 = apply(expr_m, 2, function(x){
  t1=fgsea(pathways = m_list, stats = rank(x),3,scoreType = "pos")
  t1 = t1%>%dplyr::select(pathway,ES)
})

gsea1 = Reduce(function(x,y){inner_join(x,y,by="pathway")},gsea1)
gsea1 = as.data.frame(gsea1)
gsea1 = data.frame(gsea1[,-1],row.names = gsea1[,1])
gsea1 = t(gsea1)
gsea1 = cbind.data.frame(expr1[,1:4],gsea1)

result_df_BP = data.frame()
for (i in 5:ncol(gsea1)) {
  t1 = gsea1[,c(1,2,i)]
  colnames(t1)[3]='Y'
  LM1 = lmer(Y~stage+(1|id),t1,REML = F)
  EM1 = emmeans(LM1,'stage')
  CT1 = contrast(EM1, 'tukey')
  AN1 = anova(LM1)
  EM1 = summary(EM1)
  CT1 = summary(CT1)
  t1 = data.frame(
    gene = colnames(gsea1)[i],
    stage_p = AN1$`Pr(>F)`,
    stage1_mean = EM1[1,'emmean'],
    stage2_mean = EM1[2,'emmean'],
    stage3_mean = EM1[3,'emmean'],
    stage1_se = EM1[1,'SE'],
    stage2_se = EM1[2,'SE'],
    stage3_se = EM1[3,'SE'],
    stage1v2_beta = CT1[1,'estimate'],
    stage1v3_beta = CT1[2,'estimate'],
    stage2v3_beta = CT1[3,'estimate'],
    stage1v2_se = CT1[1,'SE'],
    stage1v3_se = CT1[2,'SE'],
    stage2v3_se = CT1[3,'SE'],
    stage1v2_p = CT1[1,'p.value'],
    stage1v3_p = CT1[2,'p.value'],
    stage2v3_p = CT1[3,'p.value']
  )
  result_df_BP = rbind(result_df_BP,t1)
}
write.csv(result_df_BP,'Result/A02_DE_table_BP.csv',row.names = F)

result_df_BP = read.csv('Result/A02_DE_table_BP.csv')

plot_df = result_df_BP%>%
  arrange(stage_p)%>%
  dplyr::select(gene,stage1_mean,stage2_mean,stage3_mean)
plot_df = plot_df[1:30,]%>%as.data.frame()
plot_df = data.frame(plot_df[,-1],row.names = plot_df[,1])
colnames(plot_df)=gsub('_.*','',colnames(plot_df))

c1 = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
#pdf("Result/C07_lm_heatmap.pdf",width = 5, height = 10)
pheatmap(as.matrix(plot_df),  scale="row",color = c1)

##### direct correlate with antibody titer genes #####
df1 = ab1%>%inner_join(expr1,by=c('id','day'))
df_cor = data.frame()
for (i in 6:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,3])
  t1 = data.frame(gene = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_cor = rbind(df_cor,t1)
}
df_cor = df_cor%>%arrange(p)%>%
  mutate(p_adj = p.adjust(p,method = 'fdr'))
write.csv(df_cor,'Result/A02_cor_ab_gene.csv',row.names = F)


df1 = ab1%>%inner_join(gsea1,by=c('id','day'))
df_cor_BP = data.frame()
for (i in 6:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,3])
  t1 = data.frame(pathway = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_cor_BP = rbind(df_cor_BP,t1)
}
df_cor_BP = df_cor_BP%>%arrange(p)%>%
  mutate(p_adj = p.adjust(p,method = 'fdr'))
write.csv(df_cor_BP,'Result/A02_cor_ab_path.csv',row.names = F)


##### correlate with M4 ab #####
df1 = ab1%>%filter(day == 120)%>%na.omit()
t1 = expr1%>%filter(stage=='stage1')
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day.x, -day.y)

df_M4 = data.frame()
for (i in 5:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,2])
  t1 = data.frame(gene = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_M4 = rbind(df_M4,t1)
}
df_M4 = df_M4%>%arrange(p)
write.csv(df_M4,'Result/A02_cor_with_M4_Ab.csv',row.names = F)

#pathways
df1 = ab1%>%filter(day == 120)%>%na.omit()
t1 = gsea1%>%filter(stage=='stage1')
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day.x, -day.y)

df_M4_BP = data.frame()
for (i in 5:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,2])
  t1 = data.frame(gene = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_M4_BP = rbind(df_M4_BP,t1)
}
df_M4_BP = df_M4_BP%>%arrange(p)
write.csv(df_M4_BP,'Result/A02_cor_with_M4_Ab_BP.csv',row.names = F)

##### correlate with M4 ab all data #####
df1 = ab1%>%filter(day == 120)%>%na.omit()
t1 = expr1
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day.x, -day.y)

df_M4_all = data.frame()
for (i in 5:ncol(df1)) {
  t1 = df1[,c(2,3,i)]
  colnames(t1)[3]='x'
  LM1 = lm(antibody_value ~ stage*x,t1)
  LM2 = lm(antibody_value ~ stage,t1)
  a1 = anova(LM1,LM2)
  t1 = data.frame(gene = colnames(df1)[i],
                  p = a1$`Pr(>F)`[2])
  df_M4_all = rbind(df_M4_all,t1)
}
df_M4_all = df_M4_all%>%arrange(p)
write.csv(df_M4_all,'Result/A02_cor_with_M4_Ab_all.csv',row.names = F)

# plot 
t1 = df1[,c('id','stage','antibody_value','FOXP3')]
t1 = gather(t1,key='pathway',value = 'score',-stage,-antibody_value,-id)
t1 = t1 %>%mutate(M4_group =(antibody_value>median(antibody_value)))%>%
  mutate(M4_group=factor(M4_group,labels = c('low','high')))
p = ggplot(t1,aes(x=stage,y=score,color=M4_group))+
  geom_boxplot()+facet_wrap(~pathway, scales = 'free')
plot(p)

#pathways
df1 = ab1%>%filter(day == 120)%>%na.omit()
t1 = gsea1
df1 = df1%>%inner_join(t1,by=c('id'))%>%
  dplyr::select(-day.x, -day.y)

df_M4_all_BP = data.frame()
for (i in 5:ncol(df1)) {
  t1 = df1[,c(2,3,i)]
  colnames(t1)[3]='x'
  LM1 = lm(antibody_value ~ stage*x,t1)
  LM2 = lm(antibody_value ~ stage,t1)
  a1 = anova(LM1,LM2)
  t1 = data.frame(gene = colnames(df1)[i],
                  p = a1$`Pr(>F)`[2])
  df_M4_all_BP = rbind(df_M4_all_BP,t1)
}
df_M4_all_BP = df_M4_all_BP%>%arrange(p)
write.csv(df_M4_all_BP,'Result/A02_cor_with_M4_Ab_all_BP.csv',row.names = F)

# plot 
t1 = df1[,c('id','stage','antibody_value',df_M4_all_BP$gene[1:10])]
t1 = gather(t1,key='pathway',value = 'score',-stage,-antibody_value,-id)
t1 = t1 %>%mutate(M4_group =(antibody_value>median(antibody_value)))%>%
  mutate(M4_group=factor(M4_group,labels = c('low','high')))%>%
  mutate(pathway = factor(pathway, levels = df_M4_all_BP$gene[1:10]))
p = ggplot(t1,aes(x=stage,y=score,color=M4_group))+
  geom_boxplot()+facet_wrap(~pathway, scales = 'free')
plot(p)

##### survival analysis with virus clearance #####
t1 = expr1%>%filter(stage=='stage1')
df1 = timeTo1%>%inner_join(t1,by=c('id'))
t1 = Surv(time = df1$time2prime,event = df1$prime)
df1 = cbind.data.frame(surv = t1,df1)

df_prime = data.frame()
for (i in 13:ncol(df1)) {
  t1 = df1[,c(1,i)]
  colnames(t1)=c('surv','x')
  fit = coxph(surv~x ,data = t1)
  fit = summary(fit)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  z = fit[1,'z'],
                  p = fit[1,'Pr(>|z|)'])
  df_prime = rbind(df_prime,t1)
}
df_prime = df_prime%>%arrange(p)

# pathway
t1 = gsea1%>%filter(stage=='stage1')
df1 = timeTo1%>%inner_join(t1,by=c('id'))
t1 = Surv(time = df1$time2prime,event = df1$prime)
df1 = cbind.data.frame(surv = t1,df1)

df_prime_BP = data.frame()
for (i in 13:ncol(df1)) {
  t1 = df1[,c(1,i)]
  colnames(t1)=c('surv','x')
  fit = coxph(surv~x ,data = t1)
  fit = summary(fit)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  z = fit[1,'z'],
                  p = fit[1,'Pr(>|z|)'])
  df_prime_BP = rbind(df_prime_BP,t1)
}
df_prime_BP = df_prime_BP%>%arrange(p)

##### survival analysis with progress #####
t1 = expr1%>%filter(stage=='stage1')
df1 = timeTo1%>%inner_join(t1,by=c('id'))
t1 = Surv(time = df1$time2prog,event = df1$progression)
df1 = cbind.data.frame(surv = t1,df1)

df_prog = data.frame()
for (i in 13:ncol(df1)) {
  t1 = df1[,c(1,i)]
  colnames(t1)=c('surv','x')
  fit = coxph(surv~x ,data = t1)
  fit = summary(fit)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  z = fit[1,'z'],
                  p = fit[1,'Pr(>|z|)'])
  df_prog = rbind(df_prog,t1)
}
df_prog = df_prog%>%arrange(p)

# pathway
t1 = gsea1%>%filter(stage=='stage1')
df1 = timeTo1%>%inner_join(t1,by=c('id'))
t1 = Surv(time = df1$time2prog,event = df1$progression)
df1 = cbind.data.frame(surv = t1,df1)

df_prog_BP = data.frame()
for (i in 13:ncol(df1)) {
  t1 = df1[,c(1,i)]
  colnames(t1)=c('surv','x')
  fit = coxph(surv~x ,data = t1)
  fit = summary(fit)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  z = fit[1,'z'],
                  p = fit[1,'Pr(>|z|)'])
  df_prog_BP = rbind(df_prog_BP,t1)
}
df_prog_BP = df_prog_BP%>%arrange(p)

##### survival analysis with symptome resolve #####
t1 = expr1%>%filter(stage=='stage1')
df1 = timeTo1%>%inner_join(t1,by=c('id'))
t1 = Surv(time = df1$time2res,event = df1$resolution)
df1 = cbind.data.frame(surv = t1,df1)

df_res = data.frame()
for (i in 13:ncol(df1)) {
  t1 = df1[,c(1,i)]
  colnames(t1)=c('surv','x')
  fit = coxph(surv~x ,data = t1)
  fit = summary(fit)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  z = fit[1,'z'],
                  p = fit[1,'Pr(>|z|)'])
  df_res = rbind(df_res,t1)
}
df_res = df_res%>%arrange(p)

# pathway
t1 = gsea1%>%filter(stage=='stage1')
df1 = timeTo1%>%inner_join(t1,by=c('id'))
t1 = Surv(time = df1$time2res,event = df1$resolution)
df1 = cbind.data.frame(surv = t1,df1)

df_res_BP = data.frame()
for (i in 13:ncol(df1)) {
  t1 = df1[,c(1,i)]
  colnames(t1)=c('surv','x')
  fit = coxph(surv~x ,data = t1)
  fit = summary(fit)$coefficients
  t1 = data.frame(gene = colnames(df1)[i],
                  z = fit[1,'z'],
                  p = fit[1,'Pr(>|z|)'])
  df_res_BP = rbind(df_res_BP,t1)
}
df_res_BP = df_res_BP%>%arrange(p)


##### plot #####
t1 = df_cor[1:30,]%>%arrange(cor)
p = ggplot(t1,aes(x=reorder(gene,cor),y = cor,fill = (cor>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('genes')
plot(p)

t1 = df_cor_BP[1:30,]%>%arrange(cor)
p = ggplot(t1,aes(x=reorder(pathway,cor),y = cor,fill = (cor>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('pathways')
plot(p)

t1 = df_M4[1:30,]%>%arrange(cor)
p = ggplot(t1,aes(x=reorder(gene,cor),y = cor,fill = (cor>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('genes')
plot(p)

t1 = df_M4_BP[1:30,]%>%arrange(cor)
p = ggplot(t1,aes(x=reorder(gene,cor),y = cor,fill = (cor>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('pathway')
plot(p)

t1 = df_M4_all[1:30,]%>%arrange(p)
p = ggplot(t1,aes(x=reorder(gene,p),y = -log10(p)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('genes')
plot(p)

t1 = df_M4_all_BP[1:30,]%>%arrange(p)
p = ggplot(t1,aes(x=reorder(gene,p),y = -log10(p)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('genes')
plot(p)

t1 = df_prime[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('gene')
plot(p)

t1 = df_prime_BP[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('pathway')
plot(p)

t1 = df_prog[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('gene')
plot(p)

t1 = df_prog_BP[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('pathway')
plot(p)

t1 = df_res[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('gene')
plot(p)

t1 = df_res_BP[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('pathway')
plot(p)

##### correlate with antibody #####
ab1$day = as.numeric(ab1$day)
ab1 = ab1%>%inner_join(demo1,by= 'id')


# tx
p = ggplot(ab1,aes(x=factor(day),y = antibody_value,group=id,color=tx))+
  geom_line()
plot(p)

t1 = ab1 %>% filter(day==5)  
p = ggplot(ab1,aes(x=tx,y = antibody_value))+
  geom_violin()+geom_boxplot(width=0.5)
plot(p)

# age
p = ggplot(ab1,aes(x=day,y = antibody_value,group=id,color=(age<50)))+
  geom_line()
plot(p)

t1 = ab1 %>% filter(day==5)  
p = ggplot(ab1,aes(x=(age<50),y = antibody_value))+
  geom_violin()+geom_boxplot(width=0.5)
plot(p)

# sex
p = ggplot(ab1,aes(x=day,y = antibody_value,group=id,color=(sex==1)))+
  geom_line()
plot(p)

t1 = ab1 %>% filter(day==5)  
p = ggplot(ab1,aes(x=(sex==1),y = antibody_value))+
  geom_violin()+geom_boxplot(width=0.5)
plot(p)

# virus shedding 
t1 = timeTo1%>%dplyr::select(id, time2prime)%>%mutate(prime5 = time2prime<=5)
ab1 = ab1%>%inner_join(t1,by='id')


p = ggplot(ab1,aes(x=day,y = antibody_value,group=id,color=prime5))+
  geom_line()
plot(p)

t1 = ab1 %>% filter(day==5)  
p = ggplot(ab1,aes(x=prime5,y = antibody_value))+
  geom_violin()+geom_boxplot(width=0.5)
plot(p)

# progress 
t1 = timeTo1%>%dplyr::select(id, time2prog)%>%mutate(prog5 = time2prog<=5)
ab1 = ab1%>%inner_join(t1,by='id')


p = ggplot(ab1,aes(x=day,y = antibody_value,group=id,color=prog5))+
  geom_line()
plot(p)

t1 = ab1 %>% filter(day==5)  
p = ggplot(ab1,aes(x=prog5,y = antibody_value))+
  geom_violin()+geom_boxplot(width=0.5)
plot(p)

# resolve 
t1 = timeTo1%>%dplyr::select(id, time2res)%>%mutate(res5 = time2res<=5)
ab1 = ab1%>%inner_join(t1,by='id')

p = ggplot(ab1,aes(x=day,y = antibody_value,group=id,color=res5))+
  geom_line()
plot(p)

t1 = ab1 %>% filter(day==5)  
p = ggplot(ab1,aes(x=res5,y = antibody_value))+
  geom_violin()+geom_boxplot(width=0.5)
plot(p)
  

##### M4 vs time2event #####
load('Result/A01_compiled_data.rda')
df1 = ab1%>%filter(day == 120)%>%na.omit()
t1 = info1%>%
  dplyr::select(id,ab0,prime,time2prime,progression,
                time2prog,resolution,time2res)%>%
  unique()
df1 = df1%>%inner_join(t1,by = 'id')
df_event = df1%>%
  dplyr::select(id,ab0,antibody_value,
                prime,progression,resolution)

df_time = df1%>%
  dplyr::select(id,ab0,antibody_value,prime=time2prime,
                progression = time2prog,resolution = time2res)
df_event = df_event%>%
  gather(key='event_type',value = 'event',
         prime:resolution)

df_time = df_time%>%
  gather(key='event_type',value = 'time',
         prime:resolution)
df1 = df_event%>%
  inner_join(df_time,by=c("id", "ab0", "antibody_value",'event_type'))
df1$surv = Surv(time = df1$time,event = df1$event)

ab_surv = data.frame()
for (e1 in unique(df1$event_type)){
  df_sub = df1%>%filter(event_type==e1)
  fit = coxph(surv~ab0+antibody_value ,data =df_sub)
  fit = summary(fit)$coefficients
  t1 = data.frame(event_type = e1,
                  z = fit[2,'z'],p = fit[1,'Pr(>|z|)'])
  ab_surv = rbind(ab_surv,t1)
  df_sub = df_sub%>%
    mutate(ab120 = antibody_value>median(antibody_value))%>%
    mutate(ab120 = factor(ab120,labels = c('low','high')))
    
  km <- survfit(surv ~ ab120, 
                data = df_sub, conf.type = "log-log")
  p = ggsurvplot(km, mark.time = T)
  plot(p$plot)
}


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

##### get data #####
olink = cbind('samples'=1:nrow(olink),olink)
info1=olink[,1:3]
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

expr1 = info1%>%dplyr::select(samples)%>%inner_join(olink,by='samples')
expr1 =  t(as.matrix(expr1[,-(1:3)]))
colnames(expr1)=info1$samples

###### get gene symbol #####
prot1 = data.frame(id = rownames(expr1))
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
rownames(expr1)=prot1$symbol

##### PCA #####
expr1 = expr1[,colnames(expr1)!='180']
info1 = info1%>%filter(samples!=180)
t1 = expr1#%>%apply(2,scale)#%>%t()
pca1 = prcomp(t(t1))$x
pca1 = cbind.data.frame(info1,pca1[,1:2])
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

##### differential expression #####
expr1 = t(expr1)
expr1 = cbind.data.frame(id = info1$id,day=info1$day,
                         ab0 = info1$ab0,expr1)
expr1 = cbind.data.frame(stage =NA,expr1)
expr1$stage[expr1$day==0&expr1$ab0=='Ab0_neg']='stage1'
expr1$stage[expr1$day==5&expr1$ab0=='Ab0_neg']='stage2'
expr1$stage[expr1$day==0&expr1$ab0=='Ab0_pos']='stage2'
expr1$stage[expr1$day==5&expr1$ab0=='Ab0_pos']='stage3'
expr1 = expr1%>%na.omit()

result_df = data.frame()
for (i in 5:ncol(expr1)) {
  t1 = expr1[,c(1,2,i)]
  colnames(t1)[3]='Y'
  LM1 = lmer(Y~stage+(1|id),t1,REML = F)
  EM1 = emmeans(LM1,'stage')
  CT1 = contrast(EM1, 'tukey')
  AN1 = anova(LM1)
  EM1 = summary(EM1)
  CT1 = summary(CT1)
  t1 = data.frame(
    gene = colnames(expr1)[i],
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
  result_df = rbind(result_df,t1)
}
write.csv(result_df,'Result/A04_DE_table.csv',row.names = F)

result_df = read.csv('Result/A04_DE_table.csv')

##### heatmap #####
result_df = result_df%>%
  mutate(stage_p_adj = p.adjust(stage_p,method = 'fdr'))
plot_df = result_df%>%
  arrange(stage_p)%>%
  filter(stage_p_adj<0.05)%>%
  dplyr::select(gene,stage1_mean,stage2_mean,stage3_mean)
plot_df = plot_df[1:min(30,nrow(plot_df)),]%>%as.data.frame()
plot_df = data.frame(plot_df[,-1],row.names = plot_df[,1])
colnames(plot_df)=gsub('_.*','',colnames(plot_df))

c1 = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
#pdf("Result/C07_lm_heatmap.pdf",width = 5, height = 10)
pheatmap(as.matrix(plot_df),  scale="row",color = c1)

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
t1 = sapply(m_list, function(x){
  sum(!is.na(intersect(x,prot1$symbol)))
  })
m_list = m_list[t1>=5]
m_list = m_list[!grepl('^regulation_of|negative_regu|positive_regu',names(m_list))]


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
write.csv(result_df_BP,'Result/A04_DE_table_BP.csv',row.names = F)

result_df_BP = read.csv('Result/A04_DE_table_BP.csv')

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
df1 = ab1%>%mutate(day = factor(day))%>%
  inner_join(expr1,by=c('id','day'))
df_cor = data.frame()
for (i in 6:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,3])
  t1 = data.frame(gene = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_cor = rbind(df_cor,t1)
}
df_cor = df_cor%>%arrange(p)%>%
  mutate(p_adj = p.adjust(p,method = 'fdr'))
write.csv(df_cor,'Result/A04_cor_ab_gene.csv',row.names = F)


df1 = ab1%>%mutate(day = factor(day))%>%
  inner_join(gsea1,by=c('id','day'))
df_cor_BP = data.frame()
for (i in 6:ncol(df1)) {
  t1 = cor.test(df1[,i],df1[,3])
  t1 = data.frame(pathway = colnames(df1)[i],
                  cor = t1$estimate,p = t1$p.value)
  df_cor_BP = rbind(df_cor_BP,t1)
}
df_cor_BP = df_cor_BP%>%arrange(p)%>%
  mutate(p_adj = p.adjust(p,method = 'fdr'))
write.csv(df_cor_BP,'Result/A04_cor_ab_path.csv',row.names = F)


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
write.csv(df_M4,'Result/A04_cor_with_M4_Ab.csv',row.names = F)

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
write.csv(df_M4_BP,'Result/A04_cor_with_M4_Ab_BP.csv',row.names = F)


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
t1 = df_cor[1:30,]%>%filter(p_adj<0.05)%>%arrange(cor)
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

t1 = df_prime[1:30,]%>%arrange(z)
p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  xlab('gene')
plot(p)

# t1 = df_prime_BP[1:30,]%>%arrange(z)
# p = ggplot(t1,aes(x=reorder(gene,z),y = z,fill = (z>0)))+
#   geom_bar(stat = 'identity')+
#   coord_flip()+
#   xlab('pathway')
# plot(p)

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


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
library(lme4)
library(lmerTest)
library(emmeans)
library(gplots)
library(stringr)
library(corrplot)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))|colnames(RNAseq)%in%t1]
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


###### combined with T cell data #####
tcell = fread('Data/Tcell_Data/Lambda D28 T cell associated with antibody.csv',
              data.table=F)
tcell = tcell%>%gather(key = 'cell',value = 'cell_value',5:8)
tcell = tcell%>%mutate(tcell_var = paste(Stim,cell,sep = '_'))%>%
  filter(tcell_var%in%c('MN_MCD4/IFNg-IL21-TNF+',"SS1_MCD4/IFNg-IL21-TNF+",
                        "MEDIA_CD4+/Tfh/ICOS+",
                        "MN_CD4+/Tfh/CD137+OX40+","SS1_CD4+/Tfh/CD137+OX40+",
                        "MN_CD4+CD45RA-/CD137+OX40+","SS1_CD4+CD45RA-/CD137+OX40+"))
tcell = tcell %>%select(id,tcell_var,cell_value)
RNAseq  = RNAseq[,grep('GO_|^olink_|Info_day|Info_id',colnames(RNAseq))]
RNAseq = RNAseq%>%gather(key = 'gene',value = 'gene_value',-Info_day,-Info_id)
RNAseq = tcell%>%inner_join(RNAseq,by=c('id'='Info_id'))

f1 = function(Info_day,cell_value,gene_value){
  # t1 = RNAseq%>%
  #   filter(gene == "GO_activated_T_cell_proliferation"&
  #            tcell_var == "MEDIA_CD4+/Tfh/ICOS+")
  # Info_day = t1$Info_day
  # cell_value = t1$cell_value
  # gene_value = t1$gene_value
  
  t1 = data.frame(Info_day=Info_day,cell_value=cell_value,
                  gene_value=gene_value)%>%na.omit()
  LM1 = lm(cell_value~poly(Info_day,1)+gene_value,t1)
  LM1 = summary(LM1)
  
  t1 = data.frame(beta = LM1$coefficients['gene_value','Estimate'], 
                  t= LM1$coefficients['gene_value','t value'],
                  p = LM1$coefficients['gene_value','Pr(>|t|)'])
  return(t1)
}
dfResult = RNAseq%>%group_by(gene,tcell_var)%>%
  summarise(result = f1(Info_day,cell_value,gene_value))
dfResult  = dfResult %>%tidyr::unpack(cols=result)



dfResult = dfResult%>%
  filter(grepl('olink|GO_',gene))%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  mutate(nc = nchar(gene))%>%
  filter(nc<40)
write.csv(dfResult,'Result/A21_associations.csv',row.names = F)

sig1 = dfResult$gene[dfResult$pAdj<0.01]
dfResult = dfResult%>%filter(gene%in%sig1)
dfResult$gene = as.character(dfResult$gene)


##### plot heatmap #####
t1 = dfResult%>%
  as.data.frame()%>%
  select(gene,tcell_var,t)%>%
  spread(key = tcell_var,value = t)%>%
  as.data.frame()
t1 = as.data.frame(t1[,-1],row.names = t1[,1])%>%t()

p1 = dfResult%>%
  as.data.frame()%>%
  select(gene,tcell_var,pAdj)%>%
  spread(key = tcell_var,value = pAdj)%>%
  as.data.frame()
p1 = as.data.frame(p1[,-1],row.names = p1[,1])%>%t()

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
CLrow = hclust(dist(t1),method = 'ward.D')
CLcol = hclust(dist(t(t1)),method = 'ward.D')

t1 = t1[rev(CLrow$order),CLcol$order]
p1 = p1[rev(CLrow$order),CLcol$order]
pdf('Result/A21_heatmap_all_association.pdf',width =30,height = 30)
corrplot::corrplot(as.matrix(t1),is.corr = F,sig.level = 0.05,#pch.cex=2,
                   p.mat = as.matrix(p1),insig = "label_sig",tl.srt = 45,
                   cl.pos = "b",cl.ratio = 1, col = rev(col2(50)))
dev.off()


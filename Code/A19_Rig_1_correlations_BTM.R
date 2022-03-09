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
library(corrplot)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
#t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
#       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
#RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))|colnames(RNAseq)%in%t1]
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

##### rig1 protein vs gene #####
plotDf = RNAseq[,c("olink_DDX58_O95786_OID01018","Gene_DDX58",
                   'Gene_IFNA1','Gene_IFNB1','Gene_IFNG','EVE_IFN-a2',
                   'olink_IFN-gamma_P01579_OID05547',
                   'GO_antiviral IFN signature (M75)',
                   'GO_RIG-1 like receptor signaling (M68)',
                   'GO_type I interferon response (M127)')]

colnames(plotDf)=c("Protein DDX58","Gene DDX58",
                   'Gene IFNA1','Gene IFNB1','Gene IFNG','Protein IFN-a2',
                   'Protein IFN-gamma',
                   'BTM antiviral IFN signature',
                   'BTM RIG-1 like receptor signaling',
                   'BTM type I IFN response')
cor1 = cor(plotDf,use='pair',method = 'spearman')
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

pdf('Result/A19_DDX58_corrplot_BTM.pdf',5,5)
corrplot(cor1,order = 'hclust', col = rev(col2(50)))
dev.off()

##### rig1 protein vs gene time adjust #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
plotDf =df1[,c('Info_onset',"olink_DDX58_O95786_OID01018","Gene_DDX58",
               'Gene_IFNA1','Gene_IFNB1','Gene_IFNG','EVE_IFN-a2',
               'EVE_IFNg','olink_IFN-gamma_P01579_OID05547',
               'GO_response_to_interferon-gamma',
               'GO_response_to_type_I_interferon',
               'GO_RIG-I_signaling_pathway',
               'GO_MDA-5_signaling_pathway')]

colnames(plotDf)=c('Info_onset',"olink_DDX58","Gene_DDX58",
                   'Gene_IFNA1','Gene_IFNB1','Gene_IFNG','EVE_IFN-a2',
                   'EVE_IFNg','olink_IFN-gamma',
                   'GO_response_to_interferon-gamma',
                   'GO_response_to_type_I_interferon',
                   'GO_RIG-I_signaling_pathway',
                   'GO_MDA-5_signaling_pathway')
plotDf = plotDf%>%
  mutate(row=1:n())%>%
  gather(key='gene',value='value',-Info_onset,-row)%>%
  na.omit()%>%
  group_by(gene)%>%
  mutate(value = lm(value~poly(Info_onset,2))$residuals)%>%
  spread(key = gene,value = value)%>%
  dplyr::select(-Info_onset,-row)

cor1 = cor(plotDf,use='pair',method = 'spearman')
corrplot(cor1,order = 'hclust')

p = ggplot(plotDf,aes(olink_DDX58,Gene_DDX58 ))+
  geom_point()+geom_smooth(method = 'lm')+theme_bw()
plot(p)
cor.test(plotDf$olink_DDX58,plotDf$Gene_DDX58)



##### Top correlation with DDX58 ####
df1 = RNAseq[,grep('^GO_|^olink_|EVE_',colnames(RNAseq))]
df1 = sapply(colnames(df1),
             function(x){
               cor(df1[,x],df1[,'olink_DDX58_O95786_OID01018'],
                   use='pair',method = 'spearman')
             })
df1 = data.frame(gene = names(df1),
                 cor = df1)
df1 = df1%>%arrange(desc(abs(cor)))%>%
  mutate(gene = gsub('olink_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))

p = ggplot(df1[2:11,],aes(y=reorder(gene,cor),x=cor))+
  geom_bar(stat = 'identity')+theme_bw()
pdf('Result/A19_top_cor_Rig1.pdf',width = 2,height = 3)
plot(p)
dev.off()
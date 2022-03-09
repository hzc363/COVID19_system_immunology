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

##### load data #####
load('Result/A01_compiled_data.rda')
vaccineRNAseq = read.csv('Data/GSE169159_raw_counts_with_genenames.csv')
vaccineRNAseq = data.frame(vaccineRNAseq[,-ncol(vaccineRNAseq)],
                           row.names = vaccineRNAseq[,ncol(vaccineRNAseq)])

vaccineRNAseq = t(vaccineRNAseq)
vaccineRNAseq = cbind.data.frame(Sample_ID=rownames(vaccineRNAseq),vaccineRNAseq)
vaccineRNAseq = vaccineRNAseq%>%
  mutate(Participant=gsub('_.*|PID','',Sample_ID),.after=Sample_ID)%>%
  mutate(Day=gsub('\\..*|.*_Day|.*_BL','',Sample_ID),.after=Participant)
vaccineRNAseq$Day[vaccineRNAseq$Day=='']='0'
vaccineRNAseq$Day = as.numeric(vaccineRNAseq$Day)
vaccineRNAseq$Participant = as.numeric(vaccineRNAseq$Participant)

##### get go score ####
expr1 = vaccineRNAseq[,-c(1:3)]
rownames(expr1)=vaccineRNAseq$Sample_ID
expr1 = t(expr1)
#expr1 = apply(expr1, 2, function(x){x/sum(x)})
expr1 <- vst(expr1)
expr1 = cbind.data.frame(gene_symbol = rownames(expr1),expr1)

go_df = go_df%>%group_by(go_id)%>%mutate(N=n())%>%
  filter(N>=5)%>%as.data.frame()%>%select(-N)
go1 = inner_join(go_df,expr1,by="gene_symbol")
go1 = go1%>%dplyr::select(-go_id,-gene_id,-gene_symbol)%>%
  group_by(Term)%>%summarise_all(mean)
go1 = data.frame(go1[,-1],row.names =go1$Term )
go1 = t(go1)
go1 = cbind.data.frame(Sample_ID = rownames(go1),go1)
colnames(go1)[-1]=paste0('GO_',colnames(go1)[-1])

colnames(vaccineRNAseq)[2:3]=paste0('Info_',colnames(vaccineRNAseq)[2:3])
colnames(vaccineRNAseq)[4:ncol(vaccineRNAseq)]=paste0('Gene_',colnames(vaccineRNAseq)[4:ncol(vaccineRNAseq)])

vaccineRNAseq = vaccineRNAseq%>%
  inner_join(go1,by='Sample_ID')

#### join with Olink data #####
vaccine1 = read.csv('Data/Olink_Pfizer vaccine_Prabhu_BP Lab.csv')
vaccine1 = vaccine1%>%select(-QC_Warning)
colnames(vaccine1)= gsub('_|\\-| ','',colnames(vaccine1))

colnames(vaccine1)[1:5]=c("Sample_ID","Info_Participant","Info_Day","Info_Age","Info_Gender")
colnames(vaccine1)[6:ncol(vaccine1)]=paste0('olink_',toupper(colnames(vaccine1)[6:ncol(vaccine1)]))
vaccineRNAseq = vaccineRNAseq%>%
  inner_join(vaccine1,by=c('Info_Participant','Info_Day'))

# 
# vaccineRNAseq = vaccineRNAseq%>%
#   dplyr::select(-Sample_ID.x,-Sample_ID.y)%>%
#   arrange(Info_Participant,Info_Day)%>%
#   group_by(Info_Participant,Info_Age,Info_Gender)%>%
#   summarise_all(function(x){x-x[1]})
  

write.csv(vaccineRNAseq,'Result/A24_RNA_seq_go.csv',row.names = F)


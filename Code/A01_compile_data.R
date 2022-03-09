library(data.table)
library(dplyr)
library(tidyr)
require(GO.db)
library(org.Hs.eg.db)
require(pracma)
require(qusage)
require(ggplot2)
##### questions for pras #####
# In 'Data/ClinicalData/lambda_baseline_symptoms.csv', what does sympday mean?

##### read antibody data #####
ab1 = fread('Data/Antibodies/lambda_antibodies_followup.csv',data.table=F)
ab1 = ab1%>%filter(!is.na(id))
ab1 = gather(ab1,key = 'day',value = 'antibody_value',-id)
ab1$day[ab1$day=='M4']='D120'
ab1 = ab1 %>%mutate(day = gsub('D|M','',day))%>%mutate(day = as.integer(day))

##### get time 2 measure #####
timeTo1 = fread('Data/ClinicalData/lambda_time_to_event.csv',data.table=F)
timeTo1 = timeTo1%>%dplyr::rename(id = participantId)

##### get demographics data #####
demo1 = fread('Data/ClinicalData/lambda_demographics.csv',data.table = F)
demo1 = demo1%>%
  dplyr::select(id = participant_id, sex, 
         age = age_c,bmi,race,ethnic)

##### symptom cluster #####
symp1 = fread('Data/ClinicalData/SymptomsSeverityandClusters.csv',data.table=F)

##### symptom onset #####
sympOnset = fread('Data/ClinicalData/Lambda_symptom_onset.csv',
                  data.table = F)

sympOnset = sympOnset%>%
  dplyr::select(id=participantId,onset = symptom_onset)

##### get cytokine data #####
eve1 = fread('Data/EVE_Cytokines/SampleList_Discovery.csv',data.table = F)
eve1 = eve1 %>% dplyr::rename(day=date)%>%dplyr::select(-Well)

##### olink data #####
olink = fread('Data/OLINK Proteomics/olink.csv',data.table = F)
colnames(olink)=paste(olink[1,],olink[2,],olink[3,],sep = '_')
olink = olink[-c(1:3,184:187),-c(94,95)]
olink[,-1]= apply(olink[,-1], 2, as.numeric)
olink = cbind.data.frame(day = gsub('_L.*|_P.*','',olink$`Assay_Uniprot ID_OlinkID`),olink)
olink = cbind.data.frame(id = gsub('_.*','',olink$`Assay_Uniprot ID_OlinkID`),olink)
olink = olink%>%
  dplyr::select(-`Assay_Uniprot ID_OlinkID`)%>%
  mutate(day = gsub('.*_','',day))%>%
  mutate(day = as.numeric(day))
olink$day[is.na(olink$day)]=0

olink2 = fread('Data/OLINK Proteomics/olink2.csv',data.table = F)
colnames(olink2)=paste(olink2[1,],olink2[2,],olink2[3,],sep = '_')
olink2 = olink2[-c(1:3,175:178),-c(94,95)]
olink2[,-1]= apply(olink2[,-1], 2, as.numeric)
olink2 = cbind.data.frame(day = gsub('_L.*|_P.*','',olink2$`Assay_Uniprot ID_OlinkID`),olink2)
olink2 = cbind.data.frame(id = gsub('_.*','',olink2$`Assay_Uniprot ID_OlinkID`),olink2)
olink2 = olink2%>%
  dplyr::select(-`Assay_Uniprot ID_OlinkID`)%>%
  mutate(day = gsub('.*_','',day))%>%
  mutate(day = as.numeric(day))
olink2$day[is.na(olink2$day)]=0

olink = inner_join(olink,olink2,by=c('id','day'))

##### RNA-seq data #####
RNAseq = fread('Data/RNAseq/Result02_RNAseq_count_with_reSeq.tsv',data.table = F)
RNAseq = RNAseq%>%dplyr::select(-target_id)%>%
  filter( (gene!='') & (!is.na(gene)) )%>%
  group_by(gene)%>%summarise_all(sum)%>%
  as.data.frame()
RNAseq = data.frame(RNAseq[,-1],row.names = RNAseq[,1])
RNAseq = t(RNAseq)
RNAseq = cbind.data.frame(day = gsub('.*_','',rownames(RNAseq)),RNAseq)
RNAseq = cbind.data.frame(id = gsub('L_','',rownames(RNAseq)),RNAseq)
RNAseq = RNAseq[-nrow(RNAseq),]

RNAseq = RNAseq%>%
  mutate(id = gsub('_.*','',id))%>%
  mutate(day = as.numeric(day))

##### T cell data #####
tcell1 = fread('Data/Tcell_Data/Lambda AIM D28_complete_bgsub.csv')
tcell1 = tcell1%>%
  filter(!grepl('^E2',Sample_file))%>%
  mutate(day = 28)%>%
  dplyr::select(id = SampleID,day,Stim,
                bgsub_CD4p_CD137pOX40p:bgsub_CD8p_CD8pCD45RAn_CD137pCD69p)%>%
  filter(Stim %in%c('MN','SS1'))%>%
  filter(!grepl('LRS16',id))%>%
  gather(key = 'subset',value = 'value',-id,-day,-Stim)%>%
  mutate(subset= paste(Stim,subset,sep = '_'))%>%
  dplyr::select(-Stim)%>%
  spread(key = subset,value = value)

tcell2 = fread('Data/Tcell_Data/Lambda ICS D28_complete_bgsub_correct.csv')
tcell2 = tcell2%>%
  #filter(!grepl('^E2',Sample_file))%>%
  mutate(day = 28)%>%
  dplyr::select(id = SampleID,day,Stim,
                bgsub_CD4p_CD107a:`bgsub_Memory CD8p_CD107anIFNynIL10nIL21nTNFp`)%>%
  filter(Stim %in%c('MN','SS1'))%>%
  filter(!grepl('LRS16',id))%>%
  gather(key = 'subset',value = 'value',-id,-day,-Stim)%>%
  mutate(subset= paste(Stim,subset,sep = '_'))%>%
  dplyr::select(-Stim)%>%
  spread(key = subset,value = value)
tcell = full_join(tcell1,tcell2,by=c('id','day'))

##### get immune related GO terms #####
# filter1 <- as.list(GOBPANCESTOR)
# filter1 = sapply(filter1, function(x){'GO:0002376'%in% x})
# filter1 = names(filter1)[filter1]
# filter1=c(filter1,'GO:0070265','GO:0070266','GO:0006915')
# # get all GO mappings
# go_df=toTable(GOTERM)[,c(1,3,4)]%>%
#   filter(Ontology=='BP')%>%
#   unique()%>%
#   dplyr::select(-Ontology)%>%
#   mutate(Term=gsub(" ",'_',Term))%>%
#   filter(!grepl('^regulation_of|negative_regu|positive_regu',Term))%>%
#   filter(go_id %in% filter1)
# 
# go_list = mget(go_df$go_id, org.Hs.egGO2ALLEGS, ifnotfound = NA)
# names(go_list)= go_df$go_id
# 
# # transform mapping into data frame
# go_list = lapply(names(go_list),
#                  function(x){data.frame(go_id=x,gene_id = go_list[[x]])})
# go_list = do.call(rbind, go_list)
# go_list = na.omit(go_list)
# go_df = inner_join(go_df,go_list,by='go_id')
# 
# # map gene id to symbol
# id2symb = mget(unique(go_df$gene_id),org.Hs.egSYMBOL,ifnotfound = NA)
# id2symb = data.frame(gene_id = names(id2symb),
#                      gene_symbol = unlist(id2symb))
# go_df = inner_join(go_df,id2symb,by='gene_id')

##### BTM ######
BTM=read.gmt('Data/SupplementaryData_TutorialPackage/BTM_for_GSEA_20131008.gmt')
go_df= lapply(1:length(BTM), function(i){
  data.frame(Term=names(BTM)[i],gene_symbol=BTM[[i]])})
go_df = Reduce(rbind,go_df)
go_df = go_df%>%mutate(go_id=gsub('.* \\(|)$','',Term))
go_df = go_df %>%
  #mutate(Term = gsub('\\..\\)$',')',Term))%>%
  #mutate(Term = gsub('\\(.*\\) \\(','(',Term))%>%
  #mutate(go_id = gsub('\\..\\)$',')',go_id))%>%
  unique()%>%
  filter(!grepl('TBA',Term))

id2symb = mget(unique(go_df$gene_symbol),org.Hs.egSYMBOL2EG,ifnotfound = NA)
id2symb = data.frame(gene_symbol = names(id2symb),
                     gene_id = unlist(id2symb))%>%
  na.omit()
go_df = inner_join(go_df,id2symb,by='gene_symbol')
go_df = go_df[,c("go_id","Term","gene_id","gene_symbol")]

##### Lab data #####
lab1 = read.csv('Data/ClinicalData/lambda_labs.csv')
lab1 = lab1%>%filter(!grepl('Unsche',redcap_event_name))%>%
  mutate(redcap_event_name = gsub('.* ','',redcap_event_name))%>%
  mutate(redcap_event_name = gsub('Screening','0',redcap_event_name))%>%
  dplyr::rename(id = participant_id, day = redcap_event_name)%>%
  dplyr::select(-redcap_event_namef,-random)

#####op auc #####
op1 = read.csv('Data/ClinicalData/lambda_OP_ct.csv')
op1$CT.value = as.numeric(op1$CT.value)
op1$CT.value[op1$CT.value>=40]=40
op1$CT.value = 41-op1$CT.value
op1$CT.value[is.na(op1$CT.value)]=0

p = ggplot(op1,
           aes(x=Sample.Day,y=CT.value,group=participant_id))+
  geom_line()
plot(p)

opAUC = op1%>%group_by(participant_id)%>%
  summarise(auc=trapz(Sample.Day,CT.value))%>%
  dplyr::rename(id = participant_id)

##### save data #####
t1 = timeTo1%>%dplyr::select(id,tx)%>%unique()
demo1 = demo1%>%full_join(t1,by='id')%>%mutate(tx = gsub('Saline ','',tx))
save(ab1,demo1,eve1,symp1,olink,RNAseq,
     timeTo1,tcell,go_df,lab1,sympOnset,opAUC,
     file = 'Result/A01_compiled_data.rda')



##### count data #####
RNAseq = fread('Data/RNAseq/Result02_RNAseq_count_with_reSeq.tsv',data.table = F)
RNAseq = RNAseq%>%dplyr::select(-L_294_00_reseq)
write.csv(RNAseq,'Data_submission/RNAseq_counts.csv',row.names = F)

##### basic information #####
colNames = c('Sample name','title','source name',	'organism',
             'characteristics: treatment',	'characteristics: age',	
             'characteristics: subject',	
             'characteristics: time',	
             'characteristics: sex',	'characteristics:symptom severity',
             'characteristics: viral shedding',
             'description','processed data file',
             'raw file1',	'raw file2')

dfAll = data.frame('Sample name' = colnames(RNAseq)[-c(1:2)],
                   check.names = F)%>%
  mutate(title=gsub('L_','COVID19 patient ',`Sample name`))%>%
  mutate(title=gsub('_0',' day ',title))%>%
  mutate(`source name`='Whole blood')%>%
  mutate(organism='homo sapiens')%>%
  mutate(`characteristics: subject`=gsub('L_|_0.*','',`Sample name`))%>%
  mutate(`characteristics: subject`=as.integer(`characteristics: subject`))%>%
  mutate(`characteristics: time`=gsub('.*_0','day ',`Sample name`))

##### demographic information #####
demo1 = fread('Data/ClinicalData/lambda_demographics.csv',data.table = F)
demo1 = demo1%>%
  dplyr::select(`characteristics: subject` = participant_id, 
                `characteristics: sex`=sex, 
                `characteristics: age` = age_c)%>%
  mutate(`characteristics: sex` = factor(`characteristics: sex`,
                                         labels = c('Male','Female')))
dfAll = left_join(dfAll,demo1,by=c('characteristics: subject'))

##### severity #####
symp1 = fread('Data/ClinicalData/SymptomsSeverityandClusters.csv',data.table=F)
symp1 = symp1%>%select(`characteristics: subject` = id,
                       `characteristics:symptom severity`=severitycat)
dfAll = left_join(dfAll,symp1,by=c('characteristics: subject'))

##### viral shedding #####
load('Result/A01_compiled_data.rda')
opAUC = opAUC%>%dplyr::select(`characteristics: subject` = id, 
                      `characteristics: viral shedding`=auc)
dfAll = left_join(dfAll,opAUC,by=c('characteristics: subject'))

##### data files #####
dfAll = dfAll%>%
  mutate(`processed data file`='RNAseq_counts.csv')

file1 = list.files('Data/RNAseq/MD5_files/',full.names = T)
file1 = lapply(file1, function(x){fread(x,data.table = F,header = F)})
file1 = do.call(rbind,file1)

fq1 = file1%>%select(-V1)%>%
  mutate(`Sample name` = gsub('_CK.*','',V2))%>%
  mutate(end = gsub('.*_|\\.fq.*','',V2))%>%
  mutate(end = paste0('raw file',end))%>%
  arrange(`Sample name`,end,V2)%>%
  group_by(`Sample name`,end)%>%
  summarise(V2 = V2[1])%>%
  spread(key = end,value = V2)
dfAll = left_join(dfAll,fq1,by=c('Sample name'))
write.csv(dfAll,'Result/A22_sample_table.csv',row.names = F)

fq2 = data.frame('file name'=c(fq1$`raw file1`,fq1$`raw file2`),
                 'file type' = 'fastq',
                 check.names = F)
fq2 = left_join(fq2,file1,by=c('file name'='V2'))
colnames(fq2)[3]='file checksum'

write.csv(fq2 ,'Result/A22_raw_files.csv',row.names = F)

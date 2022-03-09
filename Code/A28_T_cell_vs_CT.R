library(dplyr)
require(pracma)


##### data prep #####
# t cell data 
tcell = read.csv('Data/D5 D28_booleans_AIM.csv',
                 check.names = F)

## auc
op1 = read.csv('Data/ClinicalData/lambda_OP_ct.csv')
op1 = op1%>%select(participant_id,Sample.Day,CT.value)
op1$CT.value = as.numeric(op1$CT.value)
op1$CT.value[op1$CT.value>=40]=40
op1$CT.value = 41-op1$CT.value
op1$CT.value[is.na(op1$CT.value)]=0
op1 = op1%>%select(participant_id,Sample.Day,CT.value)


##### single day #####
tcell = inner_join(op1,tcell,
                   by=c('participant_id'='id','Sample.Day'='days'))

f1 = function(x,y){cor(x,y,method = 'spearman',
         use = 'complete')}
f2 = function(x,y){cor.test(x,y,method = 'spearman',
                       use = 'complete')$p.value}
cor1 = tcell%>%
  select(-ExperimentID)%>%
  group_by(Sample.Day,Stim)%>%
  summarise(`MCD4/IFNg+IL21-TNF+`=f1(`MCD4/IFNg+IL21-TNF+`,CT.value),
            `MCD4/IFNg+IL21-TNF-`=f1(`MCD4/IFNg+IL21-TNF-`,CT.value),
            `MCD4/IFNg-IL21+TNF-`=f1(`MCD4/IFNg-IL21+TNF-`,CT.value),
            `MCD4/IFNg-IL21-TNF+`=f1(`MCD4/IFNg-IL21-TNF+`,CT.value),
            `Tfh/ICOS+`=f1(`Tfh/ICOS+`,CT.value),
            `Tfh/CD137+OX40+`=f1(`Tfh/CD137+OX40+`,CT.value),
            `CD4+CD45RA-/CD137+OX40+`=f1(`CD4+CD45RA-/CD137+OX40+`,CT.value))

p1 = tcell%>%
  select(-ExperimentID)%>%
  group_by(Sample.Day,Stim)%>%
  summarise(`MCD4/IFNg+IL21-TNF+`=f2(`MCD4/IFNg+IL21-TNF+`,CT.value),
            `MCD4/IFNg+IL21-TNF-`=f2(`MCD4/IFNg+IL21-TNF-`,CT.value),
            `MCD4/IFNg-IL21+TNF-`=f2(`MCD4/IFNg-IL21+TNF-`,CT.value),
            `MCD4/IFNg-IL21-TNF+`=f2(`MCD4/IFNg-IL21-TNF+`,CT.value),
            `Tfh/ICOS+`=f2(`Tfh/ICOS+`,CT.value),
            `Tfh/CD137+OX40+`=f2(`Tfh/CD137+OX40+`,CT.value),
            `CD4+CD45RA-/CD137+OX40+`=f2(`CD4+CD45RA-/CD137+OX40+`,CT.value))

plotdf = cor1[,-c(1:2)]
rownames(plotdf)=paste('day',cor1$Sample.Day,cor1$Stim)

p1 = as.matrix(p1[,-c(1:2)])
#p1 = matrix(p.adjust(p1,'fdr'),nrow = nrow(p1),ncol = ncol(p1))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
rownames(p1)=rownames(plotdf)
pdf('Result/CT_vs_t_cells.pdf',width = 10,height = 6)
corrplot::corrplot(as.matrix(plotdf),is.corr = F,
                   method = 'square',col = rev(col2(50)),
                   sig.level = 0.05,
                   #p.mat = (p1),insig = "label_sig",
                   cl.pos = "b",cl.ratio = 1)
dev.off()


##### auc before and after #####
t1 = op1%>%mutate(cutoffday=5)
t2 = op1%>%mutate(cutoffday=28)

op1 = rbind(t1,t2)%>%mutate(isAfter = Sample.Day>cutoffday)%>%
  group_by(participant_id,cutoffday,isAfter)%>%
  summarise(auc=trapz(Sample.Day,CT.value))%>%
  ungroup()%>%
  mutate(isAfter=factor(isAfter,labels = c('Before','After')))
t1 = op1%>%group_by(participant_id,cutoffday)%>%
  summarise(isAfter=isAfter[1],auc = sum(auc))%>%
  mutate(isAfter='all')
op1 = rbind(op1,t1)

##### auc cor before and after #####
tcell = read.csv('Data/Tcell_Data/D5 D28_MCD4 boolean_AIM_Tfh.csv',
                 check.names = F)
tcell = inner_join(op1,tcell,
                   by=c('participant_id'='id','cutoffday'='days'))


f1 = function(x,y){cor(x,y,method = 'spearman',
                       use = 'complete')}
f2 = function(x,y){cor.test(x,y,method = 'spearman',
                            use = 'complete')$p.value}
cor2 = tcell%>%
  select(-ExperimentID)%>%
  group_by(cutoffday,Stim,isAfter)%>%
  summarise(`MCD4/IFNg+IL21-TNF+`=f1(`MCD4/IFNg+IL21-TNF+`,auc),
            `MCD4/IFNg+IL21-TNF-`=f1(`MCD4/IFNg+IL21-TNF-`,auc),
            `MCD4/IFNg-IL21+TNF-`=f1(`MCD4/IFNg-IL21+TNF-`,auc),
            `MCD4/IFNg-IL21-TNF+`=f1(`MCD4/IFNg-IL21-TNF+`,auc),
            `Tfh/ICOS+`=f1(`Tfh/ICOS+`,auc),
            `Tfh/CD137+OX40+`=f1(`Tfh/CD137+OX40+`,auc),
            `CD4+CD45RA-/CD137+OX40+`=f1(`CD4+CD45RA-/CD137+OX40+`,auc))
p1 = tcell%>%
  select(-ExperimentID)%>%
  group_by(cutoffday,Stim,isAfter)%>%
  summarise(`MCD4/IFNg+IL21-TNF+`=f2(`MCD4/IFNg+IL21-TNF+`,auc),
            `MCD4/IFNg+IL21-TNF-`=f2(`MCD4/IFNg+IL21-TNF-`,auc),
            `MCD4/IFNg-IL21+TNF-`=f2(`MCD4/IFNg-IL21+TNF-`,auc),
            `MCD4/IFNg-IL21-TNF+`=f2(`MCD4/IFNg-IL21-TNF+`,auc),
            `Tfh/ICOS+`=f2(`Tfh/ICOS+`,auc),
            `Tfh/CD137+OX40+`=f2(`Tfh/CD137+OX40+`,auc),
            `CD4+CD45RA-/CD137+OX40+`=f2(`CD4+CD45RA-/CD137+OX40+`,auc))


col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

plotdf = cor2%>%filter(isAfter=='Before')%>%as.data.frame()
rownames(plotdf)=paste('day',plotdf$cutoffday,plotdf$Stim)
plotdf = plotdf[,-c(1:3)]

psub = p1%>%filter(isAfter=='Before')%>%as.data.frame()
psub = as.matrix(psub[,-c(1:3)])

corrplot::corrplot(as.matrix(plotdf),is.corr = F,
                   method = 'square',col = rev(col2(50)),
                   sig.level = 0.05,
                   p.mat = (psub),insig = "label_sig",
                   cl.pos = "b",cl.ratio = 1)


plotdf = cor2%>%filter(isAfter=='After')%>%as.data.frame()
rownames(plotdf)=paste('day',plotdf$cutoffday,plotdf$Stim)
plotdf = plotdf[,-c(1:3)]

psub = p1%>%filter(isAfter=='After')%>%as.data.frame()
psub = as.matrix(psub[,-c(1:3)])
psub[is.na(psub)]=1

corrplot::corrplot(as.matrix(plotdf),is.corr = F,
                   method = 'square',col = rev(col2(50)))

plotdf = cor2%>%filter(isAfter=='all')%>%as.data.frame()
rownames(plotdf)=paste('day',plotdf$cutoffday,plotdf$Stim)
plotdf = plotdf[,-c(1:3)]
psub = p1%>%filter(isAfter=='all')%>%as.data.frame()
psub = as.matrix(psub[,-c(1:3)])

corrplot::corrplot(as.matrix(plotdf),is.corr = F,
                   method = 'square',col = rev(col2(50)),
                   sig.level = 0.05,
                   p.mat = (psub),insig = "label_sig",
                   cl.pos = "b",cl.ratio = 1)

#### time vs t cell ####
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
RNAseq = RNAseq%>%filter(Info_day==0)%>%
  select(Info_onset,Info_id)%>%unique()
combined = tcell%>%
  inner_join(RNAseq,by=c('participant_id'='Info_id'))%>%
  mutate(time_since_onset = 0-Info_onset)%>%
  filter(isAfter=='all')%>%
  select(time_since_onset,auc,`CD4+CD45RA-/CD137+OX40+`)
plot(combined)




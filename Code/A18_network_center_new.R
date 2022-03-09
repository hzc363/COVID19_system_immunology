##### get all results #####
severity1 = read.csv('Result/A12_severity_results.csv')
clear1 = read.csv('Result/A13_clearance_results_auc.csv')
tcell1 = read.csv('Result/A14_t_cell_results.csv')
Ab1 = read.csv('Result/A15_Ab_M4_results_new.csv')
Ab2 = read.csv('Result/A15_Ab_M7_results_new.csv')


tAll = list(severity1,clear1,tcell1,Ab1,Ab2)
f1 = function(x,y){inner_join(x,y,by='gene')}
tAll = Reduce(f1,tAll)
tAll = tAll%>%
  dplyr::select(gene,
                severity_t=t.x,severity_p=p.x,
                clear_t=t.y,clear_p=p.y,
                tcell_t = t.x.x,tcell_p=p.x.x,
                Ab_t = t.y.y,Ab_p=p.y.y,
                Ab2_t = t, Ab2_p = p,
                rig_p = p_rig.x,rig_t = t_rig.x)%>%
  filter(grepl('olink_|EVE_|GO_',gene))%>%
  mutate(type = grepl('GO_',gene))%>%
  filter(nchar(gene)<40)%>%
  group_by(type)%>%
  mutate(severity_p=p.adjust(severity_p,'fdr'),
         clear_p=p.adjust(clear_p,'fdr'),
         tcell_p=p.adjust(tcell_p,'fdr'),
         Ab_p=p.adjust(Ab_p,'fdr'),
         Ab2_p=p.adjust(Ab2_p,'fdr'),
         rig_p = p.adjust(rig_p,'fdr'))%>%
  mutate(sum = (severity_p<0.05) +(clear_p<0.05)+(tcell_p<0.05)+
           (Ab_p<0.05)+(Ab2_p<0.05))%>%
  #filter(rig_p<0.05)%>%
  filter(sum>=3)


##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
#RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))]
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

##### join with clinical data #####
neut_ab = read.csv('Data/Antibodies/D28_M7 Neut IgG Binding.csv')
df_combine = neut_ab%>%
  dplyr::select(Info_id=ID,D28_IgG=d28_IgG_AUC,D28_Neut,
                M7_IgG = M7_IgG_AUC,M7_Neut)%>%
  mutate(D28_IgG = log2(D28_IgG),
         D28_Neut=log2(D28_Neut),
         M7_IgG=log2(M7_IgG) ,
         M7_Neut=log2(M7_Neut))%>%
  mutate(Info_id = as.character(Info_id))
RNAseq = inner_join(df_combine,RNAseq,by='Info_id')

RNAseq = tcell%>%
  dplyr::select(Info_id=id,
                mCD8_CD107a=`SS1_bgsub_Memory CD4p_IFNy`)%>%
  inner_join(RNAseq,by='Info_id')

RNAseq = opAUC%>%mutate(id = as.character(id))%>%
  right_join(RNAseq,by=c('id'='Info_id'))

##### olink data #####
RNAseq = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  as.data.frame()

olink = tAll$gene[grepl('olink_',tAll$gene)]
olink = expand.grid(olink ,olink )
olink$Var1 = as.character(olink$Var1)
olink$Var2 = as.character(olink$Var2)
olink$p = NA; olink$t = NA
for (i in 1:nrow(olink)) {
  t1 = data.frame(y=RNAseq[,olink$Var1[i]],
                  x=RNAseq[,olink$Var2[i]],
                  time = RNAseq[,'Info_onset'])
  t1 = lm(y~poly(time,2)+x,t1)
  t1 = summary(t1)
  olink$p[i]=t1$coefficients['x','Pr(>|t|)']
  olink$t[i]=t1$coefficients['x','t value']
  
}
olink$p[olink$Var1==olink$Var2] = 1
olink$p = p.adjust(olink$p,'fdr')
connect1 = olink%>%
  mutate(sig = p<0.05)%>%
  group_by(Var1)%>%
  summarise(N = sum(sig))

t1 = olink%>%
  filter(p<0.0000000005)
g <- graph_from_data_frame(t1, directed=F)
plot(g)

##### BN #####
olink = tAll$gene[grepl('olink|EVE|GO',tAll$gene)]
# 'auc','M7_Neut','D28_Neut',
#'Info_severitycat','mCD8_CD107a','Info_onset'
dfBN = RNAseq[,c(olink)]%>%
  na.omit()#%>%
  #mutate(Info_severitycat = grepl('^Se',Info_severitycat)*1)

#t1 = colnames(dfBN)[!grepl('severitycat|onset',colnames(dfBN))]
#dfBN[,t1] = apply(dfBN[,t1], 2, 
#                  function(x){lm(x~poly(dfBN$Info_onset,2))$residuals})
#dfBN = dfBN%>%dplyr::select(-Info_onset)

library(bnlearn)
bn = hc(dfBN)
g1 = graphviz.plot(bn, shape = "rectangle", layout = "fdp")

##### timeline ######
olink = tAll$gene[grepl('olink|GO',tAll$gene)&
                    !grepl('peripheral',tAll$gene)]

dfPlot = RNAseq[,c(olink)]
time1 = RNAseq[,'Info_onset']
timeDf = data.frame()
for (i in 1:(ncol(dfPlot)*1000)) {
  g1 = sample(1:ncol(dfPlot),1)
  t1 = data.frame(time1 = time1,y = dfPlot[,g1])
  t1 = na.omit(t1)
  t1 = t1%>%sample_n(nrow(t1),replace = T)
  t1 = lm(y~poly(time1,2),t1)
  t1 = which.max(predict(t1,data.frame(time1=seq(0,15,0.1))))
  t1 = seq(0,15,0.1)[t1]
  t1 = data.frame(gene = colnames(dfPlot)[g1],max = t1)
  timeDf = rbind(timeDf,t1)
}

protCL=read.csv('Result/A11_gene_CL.csv')
timeRange = timeDf%>%
  mutate(isSig = gene %in% protCL$prot)%>%
  mutate(gene = gsub('olink_|GO_|EVE_','',gene))%>%
  mutate(gene = gsub('_Q.*|_P.*|_O.*','',gene))%>%
  mutate(gene = gsub('_',' ',gene))%>%
  arrange(max)%>%
  group_by(gene,isSig)%>%
  summarise(peakMean = mean(max),
            peak5 = quantile(max,0.05),
            peak95 = quantile(max,0.95),
            peakSD = sd(max))


p = ggplot(data=timeRange, aes(x=peakMean, y=reorder(gene,peakMean), 
                               xmin=peakMean-peakSD, xmax=peakMean+peakSD))+
  geom_pointrange() +theme_bw()+
  geom_text(aes(x = 18,label=factor(isSig,labels = c('','*'))),size = 5)
pdf('Result/A18_time_range.pdf',width = 6,height = 7)
plot(p)
dev.off()


t1 = data.frame(time1 = time1,y = dfPlot$GO_B_cell_activation)
t1 = na.omit(t1)
t1 = t1%>%sample_n(nrow(t1),replace = T)
p = ggplot(t1,aes(x=time1,y=y))+
  geom_point()+
  geom_smooth(method = 'lm',formula = y~poly(x,2),se = F)+
  theme_bw()
pdf('Result/A18_example_fit.pdf',width = 2,height = 2)
plot(p)
dev.off()




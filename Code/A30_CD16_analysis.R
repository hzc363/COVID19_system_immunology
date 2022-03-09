library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

##### prepare data#####
CD16Df = read.csv('Data/Result03_CD16_count.csv')

CD16Df$count[CD16Df$count<0]=0
CD16Df = CD16Df%>%
  mutate(id = gsub('.*L_','',fn))%>%
  mutate(id = gsub('_.*','',id))%>%
  mutate(id = as.integer(id))%>%
  mutate(sample = gsub('_CKD.*','',fn))%>%
  mutate(sample = gsub('.*\\/','',sample))%>%
  mutate(day = gsub('.*_','',sample))%>%
  mutate(CD16 = gsub('_RC','',seqID))%>%
  mutate(isRC = grepl('RC',seqID))%>%
  mutate(fqPair = gsub('.*_','',fn))%>%
  mutate(fqPair = gsub('\\..*','',fqPair))
  
#### grep vs kalisto #####
RNAseq = fread('Data/RNAseq/Result02_RNAseq_count_with_reSeq.tsv',data.table = F)
RNAseq = RNAseq[,-1]%>%filter(gene=="FCGR3A")%>%
  group_by(gene)%>%summarise_all(sum)

RNAseq = data.frame(sample = colnames(RNAseq)[-1],
                    map_count =unlist(RNAseq[1,-1]))


CD16Df_sum = CD16Df%>%
  group_by(sample)%>%summarise(count=sum(count))
compare1 = inner_join(CD16Df_sum,RNAseq,by='sample')

cor(compare1$count,compare1$map_count)
plot(compare1$count,compare1$map_count)

##### compare severity, t cell and ab #####
CD16Df_sum = CD16Df%>%
  filter(isRC==F)%>%
  group_by(id,CD16)%>%
  summarise(count = sum(count))%>%
  group_by(id)%>%
  mutate(total = sum(count))%>%
  mutate(percent = count/total)%>%
  select(CD16, percent, id,total)%>%
  spread(key = CD16,value = percent)

p = ggplot(CD16Df_sum,aes(x=reorder(id,CD16_T),y=CD16_T))+
  geom_bar(stat = 'identity')
plot(p)

p = ggplot(CD16Df_sum,aes(x=reorder(id,CD16_G),y=CD16_G))+
  geom_bar(stat = 'identity')
plot(p)

CD16Df_sum$genotype = 'GT'
CD16Df_sum$genotype[CD16Df_sum$CD16_G>0.75]='GG'
CD16Df_sum$genotype[CD16Df_sum$CD16_G<0.25]='TT'


load('Result/A01_compiled_data.rda')
neut_ab = read.csv('Data/Antibodies/D28_M7 Neut IgG Binding.csv')
neut_ab = neut_ab%>%select(id=ID,D28_Neut)%>%
  mutate(id = as.integer(id))

tcell = tcell%>%
  dplyr::select(id,mCD8_IFNg=`SS1_bgsub_Memory CD4p_IFNy`)%>%
  mutate(id = as.integer(id))

CD16Df_sum = CD16Df_sum%>%
  left_join(symp1[,c('id','severitycat')],by='id')%>%
  left_join(neut_ab,by='id')%>%
  left_join(tcell,by='id')

cor(CD16Df_sum$CD16_T,CD16Df_sum$D28_Neut,method = 'spearman',use = 'complete')  
p = ggplot(CD16Df_sum,aes(x=D28_Neut,y=CD16_T))+
  geom_point()
plot(p)

cor(CD16Df_sum$CD16_T,CD16Df_sum$mCD8_IFNg,method = 'spearman',use = 'complete')  
p = ggplot(CD16Df_sum,aes(x=mCD8_IFNg,y=CD16_T))+
  geom_point()
plot(p)

t1 = lm(mCD8_IFNg~genotype,CD16Df_sum)
summary.aov(t1)
p = ggplot(CD16Df_sum,aes(x=genotype,y=mCD8_IFNg,color=genotype))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)
plot(p)

t1 = lm(D28_Neut~genotype,CD16Df_sum)
summary.aov(t1)
p = ggplot(CD16Df_sum,aes(x=genotype,y=D28_Neut,color=genotype))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)
plot(p)

plot_df = CD16Df_sum%>%group_by(genotype,severitycat)%>%summarise(N = n())
p = ggplot(plot_df,aes(x=genotype,y=N,fill=severitycat))+
  geom_bar(stat = 'identity')
plot(p)

plot_df = CD16Df_sum%>%
  group_by(genotype,severitycat)%>%
  summarise(N = n())%>%
  group_by(genotype)%>%
  mutate(total = sum(N))%>%
  mutate(percent = N/total*100)

p = ggplot(plot_df,aes(x=genotype,y=percent,fill=severitycat))+
  geom_bar(stat = 'identity')
plot(p)

p = ggplot(plot_df,aes(x=genotype,y=percent,fill=severitycat=='Severe'))+
  geom_bar(stat = 'identity')
plot(p)


CD16Df_sum = CD16Df_sum%>%mutate(is_severe = severitycat=='Severe')
L1 = glm(is_severe~genotype,
          data=CD16Df_sum,family = 'binomial')
L2 = glm(is_severe~1,
         data=CD16Df_sum,family = 'binomial')
anova(L1,L2,test = 'LRT')



p = ggplot(CD16Df_sum,aes(x=severitycat,y=CD16_T,color=severitycat))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)
plot(p)

p = ggplot(CD16Df_sum,aes(x=severitycat,y=CD16_G,color=severitycat))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)
plot(p)

LM = lm(CD16_T~severitycat,CD16Df_sum)
summary.aov(LM)

##### day vs CD16 #####
CD16Df_sum = CD16Df%>%
  group_by(day,id,CD16)%>%
  summarise(count = sum(count))%>%
  group_by(id,day)%>%
  mutate(total = sum(count))%>%
  mutate(percent = count/total)%>%
  select(CD16, percent, id,total)%>%
  spread(key = CD16,value = percent)

p = ggplot(CD16Df_sum,aes(x=day,y=CD16_T,group=id))+
  geom_point()+geom_line()
plot(p)

fit = lm(CD16_T~day+id,CD16Df_sum)
summary.aov(fit)

##### fqtype vs RC #####
CD16Df_sum = CD16Df%>%
  group_by(sample,fqPair,isRC)%>%
  summarise(count = sum(count))%>%ungroup()%>%
  mutate(seqencing_direction = factor(fqPair,labels = c('Reverse','Forward')))%>%
  mutate(matched_sequence = factor(isRC,labels = c('Foward','Reverse Complement')))

p = ggplot(CD16Df_sum,aes(x=seqencing_direction,y=count,color=matched_sequence))+
  geom_boxplot()
plot(p)



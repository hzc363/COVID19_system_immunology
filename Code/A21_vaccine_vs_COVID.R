library(dplyr)
library(tidyr)
library(ggplot2)

##### load data #####
IC1 = read.csv('Data/COVID_vaccine_IC50.csv')
EC1 = read.csv('Data/COVID_vaccine_EC50.csv')
CD4_1 = read.csv('Data/COVID_vaccine_CD4.csv')
CD8_1 = read.csv('Data/COVID_vaccine_CD8.csv')
olink1 = read.csv('Data/Olink_Pfizer vaccine_Prabhu_BP Lab.csv')

##### join data #####
# take the max of olink measure
olink1  =olink1 %>%
  rename(PID = Participant)%>%
  filter(Day==21)%>%
  select(-Day,-Sample_ID,-QC_Warning)%>%
  group_by(PID,Age,Gender)%>%
  summarise_all(max)

# gather olink data
olink1 = olink1%>%gather(key = 'olink_protein',value = 'olink_value',IL8:CSF_1)

# join with IC50
IC1 = IC1%>%select(PID, IC50_D42=D_42)
olink1 = olink1%>%full_join(IC1,by='PID')

# join with EC50
EC1 = EC1%>%select(PID, EC50_D42=D42)#%>%
  #mutate(EC50_D42=log2(EC50_D42))
olink1 = olink1%>%full_join(EC1,by='PID')

# join with CD4
CD4_1  = CD4_1 %>%select(PID, CD4_D28=D28)
olink1 = olink1%>%full_join(CD4_1, by='PID')

# join with CD8
CD8_1  = CD8_1 %>%select(PID, CD8_D28=D28)
olink1 = olink1%>%full_join(CD8_1, by='PID')

# gather t and b response
olink1 = olink1%>%
  gather(key='outcome',value = 'outcome_value',IC50_D42:CD8_D28)
olink1 = olink1%>%na.omit()

##### calculate the t value #####
f1 = function(outcome_value,olink_value){
  fit = lm(outcome_value~olink_value)
  fit = summary(fit)
  t1 = fit$coefficients[2,,drop=F]%>% as.data.frame()
  return(t1)
}

result = olink1%>%group_by(olink_protein, outcome)%>%
  summarise(LM = f1(outcome_value,olink_value))%>%
  unpack(LM)%>%
  group_by(outcome)%>%
  mutate(pAdj = p.adjust(`Pr(>|t|)`,'fdr'))%>%
  mutate(olink_protein = gsub(' |-|_','',olink_protein))%>%
  mutate(olink_protein = toupper(olink_protein))

##### compare with COVID-19 data #####
tcell1 = read.csv('Result/A14_t_cell_results.csv')
Ab1 = read.csv('Result/A15_Ab_M4_results_new.csv')

tcell1 = tcell1%>%filter(grepl('olink',gene))%>%
  mutate(gene = gsub('olink_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  mutate(gene = gsub(' |-|_','',gene))%>%
  mutate(gene = toupper(gene))%>%
  mutate(covid_pAdj= p.adjust(p,'fdr'))%>%
  select(gene,t_covid = t,covid_pAdj)%>%
  mutate(covid_outcome = 'CD4_D28')

Ab1 = Ab1%>%filter(grepl('olink',gene))%>%
  mutate(gene = gsub('olink_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  mutate(gene = gsub(' |-|_','',gene))%>%
  mutate(gene = toupper(gene))%>%
  mutate(covid_pAdj= p.adjust(p,'fdr'))%>%
  select(gene,t_covid = t,covid_pAdj)%>%
  mutate(covid_outcome = 'EC50_D42')

covid_result = rbind(tcell1,Ab1)

result = result%>%
  inner_join(covid_result ,by=c('olink_protein'='gene',
                                'outcome'='covid_outcome'))

label1 = result%>%
  filter(covid_pAdj<0.05)%>%
  group_by(outcome)%>%
  filter(`Pr(>|t|)`<0.05)
  #mutate(pAdj = p.adjust(`Pr(>|t|)`,'fdr'))#%>%
  #filter(pAdj<0.05)

p = ggplot(result,aes(y = `t value`,x = t_covid))+
  geom_point()+ #aes(color= covid_pAdj<0.05)
  theme_bw()+geom_smooth(method = 'lm')+
  ggrepel::geom_label_repel(data = label1,aes(label=olink_protein),size =3)+
  #geom_text(data = label1,aes(y = `t value`,x = t_covid,label=olink_protein))+
  facet_wrap(~outcome,scales = 'free')+
  xlab('Association (t stat) in COVID-19 cohort')+
  ylab('Association (t stat) in  vaccine cohort')
pdf('Result/A21_vaccine_vs_covid.pdf',height = 5,width = 8)
plot(p)
dev.off()

cor1 = result%>%group_by(outcome)%>%
  summarise(cor = cor(`t value`,t_covid,method = 'spearman'),
            p = cor.test(`t value`,t_covid,method = 'spearman')$p.value)

##### plot individual protein #####


IC1 = read.csv('Data/COVID_vaccine_IC50.csv')
EC1 = read.csv('Data/COVID_vaccine_EC50.csv')
CD4_1 = read.csv('Data/COVID_vaccine_CD4.csv')
CD8_1 = read.csv('Data/COVID_vaccine_CD8.csv')
olink1 = read.csv('Data/Olink_Pfizer vaccine_Prabhu_BP Lab.csv')

EC1 = EC1%>%select(PID, D42)
olink1 = olink1%>%
  filter(Day==21)%>%
  select(PID=Participant,TWEAK)

dfPlot = inner_join(EC1,olink1,by='PID')
p = ggplot(dfPlot,aes(TWEAK, log2(D42)))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)  

cor.test(log2(dfPlot$D42),dfPlot$TWEAK,method = 'spearman')  


CD4_1 = CD4_1%>%select(PID, D28)

dfPlot = inner_join(CD4_1,olink1,by='PID')
p = ggplot(dfPlot,aes(TWEAK, D28))+
  geom_point()+theme_bw()+geom_smooth(method = 'lm')
plot(p)  

cor.test(dfPlot$D28,dfPlot$TWEAK,method = 'spearman')  



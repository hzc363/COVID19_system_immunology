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
##### get gene list #####
load('Result/A01_compiled_data.rda')

pathways = c('T-helper_1_type_immune_response','complement_activation',
             'immunoglobulin_biosynthetic_process',
             'B_cell_activation','lymphocyte_mediated_immunity')

geneDf = go_df%>%filter(Term%in% pathways)%>%unique()
write.csv(geneDf,'A27_top_go_genes.csv',row.names = T)

geneDf$Term%>%table()

##### load data #####
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
  filter(grepl('Gene_',gene))

tAll = tAll%>%
  mutate(gene = gsub('Gene_','',gene))%>%
  inner_join(geneDf,by=c('gene'='gene_symbol'))%>%
  group_by(Term)%>%
  mutate(severity_p=p.adjust(severity_p,'fdr'),
         clear_p=p.adjust(clear_p,'fdr'),
         tcell_p=p.adjust(tcell_p,'fdr'),
         Ab_p=p.adjust(Ab_p,'fdr'),
         Ab2_p=p.adjust(Ab2_p,'fdr'),
         rig_p = p.adjust(rig_p,'fdr'))

t1 = tAll%>%
  select(gene,severity_t,clear_t,tcell_t,Ab_t,Ab2_t,Term)%>%
  gather(key = 'clinical',value = 't_stat',-gene,-Term)%>%
  mutate(clinical=gsub('_.*','',clinical))

p1 = tAll%>%
  select(gene,severity_p,clear_p,tcell_p,Ab_p,Ab2_p,Term)%>%
  gather(key = 'clinical',value = 'p_stat',-gene,-Term)%>%
  mutate(clinical=gsub('_.*','',clinical))

dfPlot = t1%>%inner_join(p1, by=c('gene','Term','clinical'))

p = ggplot(dfPlot,aes(x=t_stat))+
  geom_histogram()+
  geom_point(y=-1,aes(color=p_stat<0.05))+
  facet_wrap(clinical~Term,scales = 'free')
plot(p)

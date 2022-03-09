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
           (Ab_p<0.05))%>%
  #filter(rig_p<0.05)%>%
  filter(sum>=3)

##### olink heapmap #####
olink = tAll%>%filter(grepl('olink_',gene))%>%
  mutate(gene = gsub('olink_|EVE_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})
  
t1 = as.matrix(olink[,c('severity_t','clear_t',
                           'tcell_t','Ab_t','Ab2_t')])
rownames(t1)=olink$gene
p1 = as.matrix(olink[,c('severity_p','clear_p',
                       'tcell_p','Ab_p','Ab2_p')])
rownames(p1)=olink$gene
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
CL = hclust(dist(t1),method = 'ward.D')
t1 = t1[rev(CL$order),]
p1 = p1[rev(CL$order),]
pdf('Result/A17_heatmap_all_association.pdf',width =15,height = 15)
corrplot::corrplot(t1,is.corr = F,sig.level = 0.05,pch.cex=2,
                   p.mat = p1,insig = "label_sig",tl.srt = 45,
                   cl.pos = "b",cl.ratio = 1, col = rev(col2(50)))
dev.off()

##### plot go #####
olink = tAll%>%
  filter(grepl('GO_',gene))%>%
  filter(!grepl('peripheral',gene))%>%
  mutate(gene = gsub('GO_','',gene))%>%
  mutate(gene = gsub('_',' ',gene))%>%
  group_by(gene)%>%summarise_all(function(x){x[1]})

t1 = as.matrix(olink[,c('severity_t','clear_t',
                        'tcell_t','Ab_t','Ab2_t')])
rownames(t1)=olink$gene
p1 = as.matrix(olink[,c('severity_p','clear_p',
                        'tcell_p','Ab_p','Ab2_p')])
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
CL = hclust(dist(t1),method = 'ward.D')
t1 = t1[rev(CL$order),]
p1 = p1[rev(CL$order),]
pdf('Result/A17_heatmap_all_association_GO.pdf',width =10,height = 10)
corrplot::corrplot(t1,is.corr = F,sig.level = 0.05,pch.cex=2,
                   p.mat = p1,insig = "label_sig",tl.srt = 45,
                   cl.pos = "b",cl.ratio = 1, col = rev(col2(50)))
dev.off()

##### validation #####
olinkValid = read.csv('Data/Olink_validation.csv',check.names = F)
olinkValid = olinkValid%>%
  filter(Group!=0)%>%
  mutate(y = Group>1)%>%
  mutate(y = factor(y,labels = c('mild','severe')))%>%
  mutate(Info_age = Age)%>%
  mutate(Info_sex = (Gender=='F')*1)
olinkValid = olinkValid%>%na.omit()
t1 = tAll%>%filter(severity_p<0.05)#%>%
  #filter(clear_p>0.05&tcell_p>0.05&Ab_p>0.05)
colnames(olinkValid)=gsub('olink_','',colnames(olinkValid))
colnames(olinkValid)=gsub('_.*','',colnames(olinkValid))
olinkValid = olinkValid[,c('y',t1$gene)]
olinkValid = olinkValid%>%
  gather(key='prot',value = 'NPX',-y)
p = ggplot(olinkValid,aes(x=y,y = NPX,color=y))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  theme_bw()+facet_wrap(~prot,scales = 'free')
pdf('Result/A17_validation.pdf',width = 8,height = 5)
plot(p)
dev.off()




##### plot receptors #####
t1 = c('Gene_CCR2','Gene_CSF1R',#'Gene_LIFR','Gene_IL6ST',
       'Gene_OSMR','Gene_PDCD1','Gene_CXCR3')

tAll = list(severity1,clear1,tcell1,Ab1)
f1 = function(x,y){inner_join(x,y,by='gene')}
tAll = Reduce(f1,tAll)
tAll = tAll%>%
  dplyr::select(gene,
                severity_t=t.x,severity_p=p.x,
                clear_t=t.y,clear_p=p.y,
                tcell_t = t.x.x,tcell_p=p.x.x,
                Ab_t = t.y.y,Ab_p=p.y.y)%>%
  filter(gene%in%t1)

t1 = as.matrix(tAll[,c('severity_t','clear_t',
                       'tcell_t','Ab_t')])
rownames(t1)=tAll$gene
p1 = as.matrix(tAll[,c('severity_p','clear_p',
                       'tcell_p','Ab_p')])
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
CL = hclust(dist(p1<0.05),method = 'ward.D')
t1 = t1[rev(CL$order),]
p1 = p1[rev(CL$order),]
pdf('Result/A17_heatmap_receptor_association.pdf',width =10,height = 5)
corrplot::corrplot(t(t1),is.corr = F,sig.level = 0.05,
                   p.mat = t(p1),insig = "label_sig",
                   cl.pos = "b",cl.ratio = 1, col = rev(col2(50)))
dev.off()

##### Associations with receptors #####
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
RNAseq = RNAseq%>%mutate(Info_onset = Info_day-Info_onset)%>%
  mutate(Info_time2prime = Info_time2prime-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')
  
t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
RNAseq = RNAseq[,c('Info_severitycat','Info_onset','Info_time2prime',t1)]%>%
  gather(key = 'gene',value = 'expression',
         -Info_severitycat,-Info_onset,-Info_time2prime)
p = ggplot(RNAseq,aes(x=Info_severitycat,
                      y=expression,color = Info_severitycat))+
  geom_boxplot()+facet_wrap(~gene,scales = 'free')+
  geom_jitter(width = 0.1)
plot(p)

p = ggplot(RNAseq,aes(color=Info_severitycat,
                      y=expression,x = Info_onset))+
  geom_point()+facet_wrap(~gene,scales = 'free')+
  geom_smooth(method='lm')
plot(p)

f1 = function(y,x,xc){
  a = lm(y~x+poly(xc,2))
  a = summary.aov(a)
  a[[1]][1,'Pr(>F)']}

RNAseq%>%group_by(gene)%>%
  summarise(p = f1(expression,Info_severitycat,Info_onset))





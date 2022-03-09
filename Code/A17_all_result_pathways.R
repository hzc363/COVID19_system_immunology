##### get all results #####
severity1 = read.csv('Result/A12_severity_results.csv')
clear1 = read.csv('Result/A13_clearance_results_auc.csv')
tcell1 = read.csv('Result/A14_t_cell_results.csv')
Ab1 = read.csv('Result/A15_Ab_M4_results_new.csv')
Ab2 = read.csv('Result/A15_Ab_M7_results_new.csv')

GO_filter = c('antiviral IFN signature (M75)',
               'type I interferon response (M127)',
               'RIG-1 like receptor signaling (M68)')
#GO_filter = c('RIG-I signaling pathway','response to type I interferon',
#              'type I interferon signaling pathway','apoptotic process',
#              'necrotic cell death')


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
  filter(grepl('GO_',gene))%>%
  mutate(gene = gsub('GO_','',gene))%>%
  mutate(gene = gsub('_',' ',gene))%>%
  filter(gene%in%GO_filter)%>%
  filter(nchar(gene)<40)%>%
  # mutate(severity_p=p.adjust(severity_p,'fdr'),
  #        clear_p=p.adjust(clear_p,'fdr'),
  #        tcell_p=p.adjust(tcell_p,'fdr'),
  #        Ab_p=p.adjust(Ab_p,'fdr'),
  #        Ab2_p=p.adjust(Ab2_p,'fdr'),
  #        rig_p = p.adjust(rig_p,'fdr'))%>%
  mutate(sum = (severity_p<0.05) +(clear_p<0.05)+(tcell_p<0.05)+
           (Ab_p<0.05))#%>%
  #filter(rig_p<0.05)%>%
  #filter(sum>=3)


##### plot go #####
olink = tAll%>%
  filter(!grepl('peripheral',gene))%>%
  mutate(gene = gsub('GO_','',gene))%>%
  mutate(gene = gsub('_',' ',gene))%>%
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
#col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                           "#FDDBC7", "#FFFFFF"))
#CL = hclust(dist(t1),method = 'ward.D')
#t1 = t1[rev(CL$order),]
#p1 = p1[rev(CL$order),]
pdf('Result/A17_heatmap_all_association_GO.pdf',width =5,height = 5)
corrplot::corrplot(t1,is.corr = F,sig.level = 0.05,pch.cex=2,
                   p.mat = p1,insig = "label_sig",tl.srt = 45,
                   cl.pos = "b",cl.ratio = 1, 
                   col = rev(col2(50)),col.lim = c(-2,4))
dev.off()

##### severity vs others #####
dfPlot = RNAseq%>%
  dplyr::select(D28_IgG,M7_IgG,Info_severitycat,)







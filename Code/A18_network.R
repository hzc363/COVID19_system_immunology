##### get all results #####
severity1 = read.csv('Result/A12_severity_results.csv')%>%
  mutate(clin = 'severity')
clear1 = read.csv('Result/A13_clearance_results_auc.csv')%>%
  mutate(clin = 'virus')
tcell1 = read.csv('Result/A14_t_cell_results.csv')%>%
  mutate(clin = 'tcell')
Ab1 = read.csv('Result/A15_Ab_M4_results_new.csv')%>%
  mutate(clin = 'Ab')

tAll = list(severity1,clear1,tcell1,Ab1)
tAll = Reduce(rbind,tAll)
tAll = tAll%>%
  filter(pAdj<0.05)%>%
  filter(grepl('olink',gene))%>%
  mutate(gene = gsub('olink_','',gene))%>%
  mutate(gene = gsub('_.*','',gene))%>%
  filter(nchar(gene)<50)%>%
  dplyr::select(gene,clin)

g <- graph_from_data_frame(tAll, directed=F)

#colrs <- c("gray50", "tomato", "gold")
V(g)$color[names(V(g))%in%unique(tAll$clin) ] <- 'tomato'
V(g)$color[!names(V(g))%in%unique(tAll$clin) ] <- 'gray50'
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
l <- do.call("layout_with_dh", list(g)) 
plot(g,  layout=l)
library(data.table)

mbaa = fread('Data/41586_2020_2588_MOESM3_ESM.csv')

mbaa = mbaa%>%filter(!is.na(DFSO))%>%
  na.omit()%>%mutate(ID=gsub('\\..*','',ID))%>%
  filter(DFSO<20)%>%
  select(ID, DFSO,CXCL1=CXCL1orGROa,IFNy,TRAIL,MCP1=CCL2orMCP1,
         CXCL10=CXCL10orIP10,MCP2=CCL8OrMCP2)%>%
  group_by(ID)%>%
  mutate(N = length(unique(DFSO)))%>%
  filter(N>2)%>%
  gather(key='gene',value = 'value',CXCL1:MCP2)

p = ggplot(mbaa,aes(x=DFSO, y=value))+
  geom_point()+facet_wrap(~gene,scales = 'free')+
  geom_line(aes(group=ID),alpha=0.2)+
  geom_smooth(method = 'lm',formula = 'y~poly(x,2)')+
  theme_bw()
pdf('Result/A33_misfiring_trajectory.pdf',height = 5,width = 7)
plot(p)
dev.off()

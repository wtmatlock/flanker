library(tidyverse)

library(igraph)

z=1
# load some metadata (in this case the source of isolate eg human/cow etc)
source<-read_csv('../REHAB_isolate_source.csv')
library(RColorBrewer)

for(i in seq(100,3000,100)){
  print(i)
  q<-read.delim(paste('~/gn/flanker/mash_',i,'.tsv',sep=''),sep='\t',header=F,stringsAsFactors = F)
  names(q)<-c('X1','X2','mash_dist','X4','X5')
  x<-unique(unlist(q[c(1,2)]))
  x<-data.frame(str_replace_all(x,'.fasta',''))
  names(x)<-c('X1')
  x$name<-str_replace_all(x$X1,'_.*','')
  x<-left_join(x,source,by=c("name"="Alternative isolate name"))
  x$Niche<-ifelse(is.na(x$Niche),'human',x$Niche)
  x$Niche<-ifelse(grepl('WTP',x$Niche),'WTP',x$Niche)
  q$X1<-str_replace_all(q$X1,'.fasta.*','')
  q$X2<-str_replace_all(q$X2,'.fasta.*','')
  q<-filter(q, mash_dist == 0)
  q<-select(q,X1,X2)
  names(q)<-c('X1','X2')
 
  q<-filter(q, X1 != X2)

  q<- q %>%  rowwise() %>%
    mutate(HH=paste0(sort(c(X1,X2)), collapse =',')) %>%
    distinct(HH, .keep_all = TRUE) %>% select(-HH)
  q<-graph_from_data_frame(q,vertices = x,directed =F )
  V(q)$color<-ifelse(V(q)$Niche =='human','blue',ifelse(V(q)$Niche=='WTP','lightblue',ifelse(V(q)$Niche=='Cattle','brown',ifelse(V(q)$Niche=='Pig','pink',ifelse(V(q)$Niche=='Sheep','yellow','black')))))
  png(filename = paste(sprintf("out_%03d",z),'.png',sep=''),width = 800,height = 800)
  
  plot(q,vertex.label=NA,vertex.size=2.0,width=0.1,main=i)

  legend("bottomleft",legend=c("human","WTP","Catlle","Pig","Sheep","Other"),fill=c("blue","lightblue","brown","pink","yellow","black"))
  dev.off()
  z=z+1
  
}



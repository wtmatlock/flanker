##R script to take output from mash distance matrix and turn into png files for each window 

library(tidyverse)

library(igraph)
out=NULL
z=1
graphs=list()
for(i in seq(100,5000,100)){
  print(i)
  q<-read.delim(paste('~/gn/flanker/mash_',i,'.tsv',sep=''),sep='\t',header=F,stringsAsFactors = F)
  names(q)<-c('X1','X2','mash_dist','X4','X5')
  x<-unique(unlist(q[c(1,2)]))
  x<-str_replace_all(x,'.fasta','')
  q$X1<-str_replace_all(q$X1,'.fasta.*','')
  q$X2<-str_replace_all(q$X2,'.fasta.*','')
  q<-select(q,X1,X2,mash_dist)
  names(q)<-c('X1','X2','weights')
  q<-filter(q, weights == 0)
  q<-filter(q, X1 != X2)
  q<-graph_from_data_frame(q,vertices = x)
  png(filename = paste(i,'.png'),width = .1,height = .1)
  plot(q)
  dev.off()
  
}


# these png files can the be passed to ffmpeg eg:
# ffmpeg -framerate 1 -i out_%03d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

library(tidyverse)
library(ape)
library(igraph)
library(phytools)
library(ggtree)
library(gggenes)
library(optparse)

option_list= list(
  make_option(c("-i", "--input_directory"), type="character", default=NULL, help='lol',metavar='character'),
  make_option(c("-o", "--out_prefix"), type="character", default=NULL, help='list of taxa',metavar='character'),
  make_option(c("-t", "--tree"), type="character", default=NULL, help='file with dates for big tree, tab separated',metavar='character'),
  make_option(c("-w", "--windows"), type="character", default=NULL, help='file with dates for big tree, tab separated',metavar='character')
  
);

opt_parser= OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


flanker_output<-read_csv(paste(opt$i,'/flanker_matrix', sep=''))

flanker_output<-flanker_output %>% 
  mutate(window = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[3]))

flanker_output$window<-as.numeric(flanker_output$window)
x<-unique(flanker_output$window) %>% sort()
flanker_output<-flanker_output %>% 
  mutate(guuid = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[6]))

flanker_output<-flanker_output %>% 
  mutate(gene = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[5]))

gene<-select(flanker_output,guuid,gene)
flanker_output<-select(flanker_output,guuid,window,group)
flanker_output<-flanker_output %>% pivot_wider(id_cols = guuid,names_from=window,values_from=group,names_sort=TRUE)

flanker_output$ID<-flanker_output %>% group_indices(flanker_output[,2:opt$windows])

test<-select(flanker_output,guuid,ID) %>% distinct()

flanker_output<-select(flanker_output,-ID)
flanker_output<-flanker_output %>% pivot_longer(-guuid,names_to = "window",values_to="group")

flanker_output<-left_join(flanker_output,test,by=c("guuid"="guuid"))


flanker_output$window<-as.numeric(flanker_output$window)
flanker_output<-left_join(flanker_output,gene,by=c("guuid"="guuid"))
p0<-ggplot(flanker_output) + aes(x=window,y=interaction(guuid,ID),fill=group) +geom_tile() + theme_minimal()  + facet_wrap(~gene)

jpeg(paste(opt$out_prefix, '_1.jpg',sep = ''))
p0
dev.off()


##########################

tree<-read.tree(paste(opt$tree,'.tree',sep=''))
tree<-midpoint.root(tree)
x<-data.frame(tree$tip.label)
print(x)

x<-x %>% 
  mutate(id = map_chr(tree.tip.label, function(s) rev(strsplit(s, "_")[[1]])[1]))

tree$tip.label<-x$id
x<-left_join(x,gene,by=c("id"="guuid"))
x<-select(x,id,gene) %>% distinct()

g<-ggtree(tree)
p1<- g %<+% x + geom_tippoint(aes(color=gene)) + geom_tiplab()

px<-p1 + geom_facet(panel="flankergram", data=flanker_output,aes(x=window,fill=group), geom=geom_tile) 
names(flanker_output)<-c('id','window','group','cluster','gene')


jpeg(paste(opt$out_prefix, '_2.jpg',sep = ''))
px
dev.off()

##################################

filenames<-list.files(paste(opt$input_directory,sep=''),pattern="*.gbk",full.names = TRUE)
print(filenames)
filenames<-str_replace_all(filenames,'//','/')
files<-as.list(filenames)

library(genbankr)
files2<-lapply(files, readGenBank)


for (i in 1:length(files2)){
  assign(paste(paste("df", i, sep="_"), "summary", sep="."), data.frame(genes(files2[[i]])))
  t<-paste(paste("df", i, sep="_"), "summary", sep=".")
  df<-get(t)
  df$seqnames<-paste('g',i,sep='')
  assign(paste(paste("df", i, sep="_"), "summary", sep="."), df)
}

files2<-lapply(ls(pattern='df_*'), get)

gggene_configure<-function(seqname){
  seqname<-select(seqname,seqnames,gene,start,end,strand)
  seqname$gene<-ifelse(is.na(seqname$gene),'hypothetical',seqname$gene)
  seqname$strand<-ifelse(seqname$strand=='+','forward','reverse')
  seqname$direction<-ifelse(seqname$strand =='forward',1,-1)
  
  return(seqname)
  
}

dfs<-lapply(ls(pattern = 'df_'),get)


all<-bind_rows(dfs)
all2<-gggene_configure(all)
all<-all2
all<-select(all,seqnames,gene,start,end,strand,direction)
names(all)<-c('molecule','gene','start','end','strand','direction')
filenames<-data.frame(filenames)
filenames<-filenames %>% 
  mutate(isolate = map_chr(filenames, function(s) rev(strsplit(s, "_")[[1]])[6])) 
filenames$isolate<-str_replace_all(filenames$isolate,'.gbk','')
filenames<-select(filenames,-filenames)
filenames$molecule<-paste('g',1:nrow(filenames),sep = '')

filenames$isolate<-str_replace_all(filenames$isolate,'^.*/','')
print(filenames)

all<-left_join(all,filenames,by=c("molecule"="molecule"))
all<-select(all,isolate,gene,start,end,strand,direction)
names(all)<-c('id','gene','start','end','strand','direction')



g<-ggplot(all,aes(xmin=start,xmax=end,y=id,fill=gene,forward=direction)) +
  geom_gene_arrow() + theme_genes() 

jpeg(paste(opt$out_prefix, '_3.jpg',sep = ''))
g
dev.off()


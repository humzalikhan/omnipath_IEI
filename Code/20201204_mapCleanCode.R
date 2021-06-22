# loading ----
library(tidyverse)
library(readxl)
library(rmutil)
library(OmnipathR)
#library(biomaRt)
library(HGNChelper)
library(RColorBrewer)
library(tidyr)
library(dnet)
library(gprofiler2)

theme_manish <- function(legend="none") {
  theme_bw(base_size=10, base_line_size = 0.5, base_rect_size = 0.5) +
    theme(axis.text = element_text(size = rel(0.8), colour = "black", angle = 0)) +
    theme(axis.title = element_text(size = rel(1), colour = "black", angle = 0)) +
    theme(axis.ticks = element_line(colour = "black", size = rel(0.25))) +
    theme(panel.border = element_rect(colour = "black", size = rel(0.25))) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(panel.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(legend.position=legend) +
    theme(strip.background = element_rect(size=rel(0.25)) )
}

theme_manish_legend<-function(legend){
  theme_manish(legend=legend)
}

theme_humza <- function() {
  theme_bw(base_size=10, base_line_size = 0.5, base_rect_size = 0.5) +
    theme(axis.text = element_text(size = rel(0.8), colour = "black", angle = 0)) +
    theme(axis.title = element_text(size = rel(1), colour = "black", angle = 0)) +
    theme(axis.ticks = element_line(colour = "black", size = rel(0.25))) +
    theme(panel.border = element_rect(colour = "black", size = rel(0.25))) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(size=rel(.5))) +
    theme(panel.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(legend.position="none") +
    theme(strip.background = element_rect(size=rel(0.25)) )
}

theme_humza_legend <- function() {
  theme_bw(base_size=10, base_line_size = 0.5, base_rect_size = 0.5) +
    theme(axis.text = element_text(size = rel(0.8), colour = "black", angle = 0)) +
    theme(axis.title = element_text(size = rel(1), colour = "black", angle = 0)) +
    theme(axis.ticks = element_line(colour = "black", size = rel(0.25))) +
    theme(panel.border = element_rect(colour = "black", size = rel(0.25))) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_line()) +
    theme(panel.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(strip.background = element_rect(size=rel(0.25)) )
}

gNOMADpLIs<-read_tsv("/Users/humzakhan/Desktop/Khan_Butte_Supplement/pLI_GDI/gNOMAD_pLI.tsv")

# import data -------
pidList<-read_xlsx("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Lists/ConfirmedPIDList/PID_Genes_FilteredForGenes_ChrRemoved.xlsx")

genes<-pidList %>% as_tibble() %>% dplyr::select('Genetic defect')
genes<-unique(genes) %>% na.omit()

colnames(genes)<-"ConfirmedPIDGene"

# PTM -----
#phosphoInteractions <- import_omnipath_enzsub(resources=c("SIGNOR")) #use database freeze

phosphoInteractions<-read_csv("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Omnipath_DatabaseFreeze/20201208_PTMData.csv")

pidGeneAsKinase<-as_tibble()
for(i in c(1:nrow(genes))) {
  gene<-as.character(genes[i,1])
  phosphointeraction<-phosphoInteractions %>% dplyr::filter(enzyme_genesymbol==gene)
  pidGeneAsKinase<-rbind(phosphointeraction,pidGeneAsKinase)
}

pidGeneAsSubstrate<-as_tibble()
for(i in c(1:nrow(genes))) {
  gene<-as.character(genes[i,1])
  phosphointeraction<-phosphoInteractions %>% filter(substrate_genesymbol==gene)
  pidGeneAsSubstrate<-rbind(phosphointeraction,pidGeneAsSubstrate)
}

unique(pidGeneAsKinase %>% as_tibble() %>% dplyr::select(substrate_genesymbol))
#above is novel substrates of PID genes

unique(pidGeneAsSubstrate %>% as_tibble() %>% dplyr::select(enzyme_genesymbol))
#above is novel kinases of PID genes

pidGeneAsSubstrate %>% as_tibble() %>% dplyr::count(modification)
pidGeneAsKinase %>% as_tibble() %>% dplyr::count(modification)

pidGeneAsSubstrate2<-pidGeneAsSubstrate %>% as_tibble() %>% dplyr::select(substrate_genesymbol,enzyme_genesymbol)
pidGeneAsSubstrate2$InteractionType<-"PIDGeneIsPhosphoSubstrate"
colnames(pidGeneAsSubstrate2)<-c("PIDGene","Target","InteractionType")

pidGeneAsKinase2<-pidGeneAsKinase %>% as_tibble() %>% dplyr::select(enzyme_genesymbol,substrate_genesymbol)
pidGeneAsKinase2$InteractionType<-"PIDGeneIsPhosphoKinase"
colnames(pidGeneAsKinase2)<-c("PIDGene","Target","InteractionType")

#transcription factor --------
#tfInteractions<-import_dorothea_interactions() #use database freeze

tfInteractions<-read_csv("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Omnipath_DatabaseFreeze/20201208_TFData.csv")

pidGeneAsTF<-as_tibble()
for(i in 1:nrow(genes)){
  gene<-as.character(genes[i,1])
  interactions_AB <- dplyr::filter(tfInteractions, source_genesymbol == gene,dorothea_level=="A"|dorothea_level=="B")
  pidGeneAsTF<-rbind(interactions_AB,pidGeneAsTF)
}

pidGeneAsTF<-pidGeneAsTF %>% as_tibble() 

unique(pidGeneAsTF$target_genesymbol)

pidGeneAsTF2<-pidGeneAsTF %>% dplyr::select(source_genesymbol,target_genesymbol)
pidGeneAsTF2$InteractionType<-"PIDGeneIsTF"

colnames(pidGeneAsTF2)<-c("PIDGene","Target","InteractionType")

# targets of PID gene TFs 

pidGeneAsTranscribed<-as_tibble()
for(i in 1:nrow(genes)){
  gene<-as.character(genes[i,1])
  interactions_AB <- dplyr::filter(tfInteractions, target_genesymbol == gene,dorothea_level=="A"|dorothea_level=="B")
  pidGeneAsTranscribed<-rbind(interactions_AB,pidGeneAsTranscribed)
}

pidGeneAsTranscribed<-pidGeneAsTranscribed %>% as_tibble() 

pidGeneAsTranscribed2<-pidGeneAsTranscribed %>% dplyr::select(source_genesymbol,target_genesymbol)
pidGeneAsTranscribed2$InteractionType<-"PIDGeneIsTranscribed"

colnames(pidGeneAsTranscribed2)<-c("Target","PIDGene","InteractionType")

# binding path 1 before failed ------
TotalInteractions<-rbind(pidGeneAsKinase2,pidGeneAsSubstrate2,pidGeneAsTF2,pidGeneAsTranscribed2) #CleanedComplexComponentsList adds 4000 unique genes
TotalInteractions<-unique(TotalInteractions)

length(unique(TotalInteractions$Target))

TotalInteractions$Group<-""
for(i in c(1:nrow(TotalInteractions))){
  #print(TotalInteractions[i,]$PIDGene)
  x<-((pidList %>% ungroup() %>% dplyr::filter(`Genetic defect`==TotalInteractions[i,]$PIDGene))$`Major category`)[1]
  #print(x)
  TotalInteractions[i,]$Group<-x
}

genesAnalyzed<-unique(TotalInteractions$PIDGene)

names(TotalInteractions)[names(TotalInteractions) == 'Target'] <- 'gene'

TotalInteractions <- merge(TotalInteractions,gNOMADpLIs, by  = "gene",all.x=T) %>% as_tibble() %>% 
  dplyr::select(gene,PIDGene,InteractionType,Group,pLI)

# failed ----
genes<-genes %>% mutate(Analyzed=ConfirmedPIDGene %in% TotalInteractions$PIDGene)
failedgenes<-(genes %>% filter(Analyzed=="FALSE"))$ConfirmedPIDGene %>% as_tibble()

correctedgenes<-as_tibble()

for(i in c(1:nrow(failedgenes))){
  #print(failedgenes[i,1])
  x<-hgnc.table %>% as_tibble() %>% filter(Symbol==as.character(failedgenes[i,1]))
  #print(x)
  correctedgenes<-rbind(x,correctedgenes)
}

# correctedgenes2<-as_tibble()
# 
# for(i in c(1:nrow(failedgenes))){
#   #print(failedgenes[i,1])
#   x<-hgnc.table %>% as_tibble() %>% filter(Approved.Symbol==as.character(failedgenes[i,1]))
#   #print(x)
#   correctedgenes2<-rbind(x,correctedgenes)
# }
# 
# intersect(correctedgenes,correctedgenes2)

corrected<-correctedgenes %>% dplyr::select(Approved.Symbol)

pidGeneAsKinase3<-as_tibble()
for(i in c(1:nrow(corrected))) {
  gene<-as.character(corrected[i,1])
  phosphointeraction<-phosphoInteractions %>% dplyr::filter(enzyme_genesymbol==gene)
  pidGeneAsKinase3<-rbind(phosphointeraction,pidGeneAsKinase3)
}

pidGeneAsSubstrate3<-as_tibble()
for(i in c(1:nrow(corrected))) {
  gene<-as.character(corrected[i,1])
  phosphointeraction<-phosphoInteractions %>% filter(substrate_genesymbol==gene)
  pidGeneAsSubstrate3<-rbind(phosphointeraction,pidGeneAsSubstrate3)
}


pidGeneAsSubstrate4<-pidGeneAsSubstrate3 %>% as_tibble() %>% dplyr::select(substrate_genesymbol,enzyme_genesymbol)
pidGeneAsSubstrate4$InteractionType<-"PIDGeneIsPhosphoSubstrate"
colnames(pidGeneAsSubstrate4)<-c("PIDGene","Target","InteractionType")

pidGeneAsKinase4<-pidGeneAsKinase3 %>% as_tibble() %>% dplyr::select(enzyme_genesymbol,substrate_genesymbol)
pidGeneAsKinase4$InteractionType<-"PIDGeneIsPhosphoKinase"
colnames(pidGeneAsKinase4)<-c("PIDGene","Target","InteractionType")

# tf

pidGeneAsTF3<-as_tibble()
for(i in 1:nrow(corrected)){
  gene<-as.character(corrected[i,1])
  interactions_AB <- dplyr::filter(tfInteractions, source_genesymbol == gene,dorothea_level=="A"|dorothea_level=="B")
  pidGeneAsTF3<-rbind(interactions_AB,pidGeneAsTF3)
}

pidGeneAsTF3<-pidGeneAsTF3 %>% as_tibble() 

pidGeneAsTF4<-pidGeneAsTF3 %>% dplyr::select(source_genesymbol,target_genesymbol)
pidGeneAsTF4$InteractionType<-"PIDGeneIsTF"

colnames(pidGeneAsTF4)<-c("PIDGene","Target","InteractionType")

# targets of PID gene TFs 

pidGeneAsTranscribed3<-as_tibble()
for(i in 1:nrow(corrected)){
  gene<-as.character(corrected[i,1])
  interactions_AB <- dplyr::filter(tfInteractions, target_genesymbol == gene,dorothea_level=="A"|dorothea_level=="B")
  pidGeneAsTranscribed3<-rbind(interactions_AB,pidGeneAsTranscribed3)
}

pidGeneAsTranscribed3<-pidGeneAsTranscribed3 %>% as_tibble() 

pidGeneAsTranscribed4<-pidGeneAsTranscribed3 %>% dplyr::select(source_genesymbol,target_genesymbol)
pidGeneAsTranscribed4$InteractionType<-"PIDGeneIsTranscribed"

colnames(pidGeneAsTranscribed4)<-c("Target","PIDGene","InteractionType")

View(pidGeneAsTranscribed4)

TotalInteractions2<-TotalInteractions %>%  dplyr::select(-c(Group)) #InPID
names(TotalInteractions2)[names(TotalInteractions2)=='gene']<-"Target"
TotalInteractions<-rbind(TotalInteractions2[,1:3],pidGeneAsSubstrate4,pidGeneAsTranscribed4) 
#pidGeneAsKinase4,pidGeneAsTF4

TotalInteractions$Group<-""
for(i in c(1:nrow(TotalInteractions))){
  #print(TotalInteractions[i,]$PIDGene)
  x<-((pidList %>% ungroup() %>% dplyr::filter(`Genetic defect`==TotalInteractions[i,]$PIDGene))$`Major category`)[1]
  #print(x)
  TotalInteractions[i,]$Group<-x
}

TotalInteractions

TotalInteractions 

# pathways -----
a<-Sys.time()
# Omnipath signaling
#ppInteractionsOmniimport<- import_omnipath_interactions() #use database freeze
ppInteractions<-read_csv("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Omnipath_DatabaseFreeze/20201208_OmniPathData.csv") %>% dplyr::select(-X1)

# all together

OPI_g = interaction_graph(interactions = ppInteractions)

makePathsEasier <- function(fromGene,toGene) {
  b<-Sys.time()
  x<-all_shortest_paths(OPI_g,from = fromGene,to=toGene)
  placeholder<-NULL
  if(length(x$res)!=0){
    for(i in c(1:length(x$res))){
     
      y<-x$res[[i]] %>% as.integer() %>% as_tibble() %>% dplyr::mutate(Gene=vertex_attr(OPI_g,index=V(OPI_g))$name[value])
      
      y$Pathway<-i
      y<-y %>% dplyr::mutate(distance=dplyr::row_number())
      y<-y %>% dplyr::mutate(DistFromStart=max(y$distance))
      placeholder<-rbind(placeholder,y)
      
    }
    placeholder$Start<-fromGene
    placeholder$End<-toGene
    a<-Sys.time()
    print(a-b)
    return(placeholder)
  }
}

makePathsEvenEasier <- function(fromGene,toGene) {
  
 # a<-Sys.time()
  x<-all_shortest_paths(OPI_g,from = fromGene,to=toGene)
 # c<-Sys.time()
  #print(paste("allshortpath",c-a))

  if(length(x$res)!=0){
    df <- as.data.frame(t(sapply(x$res, as_ids)))
    cols<-colnames(df)
    last<-cols[length(cols)]
    df<-df %>% mutate(DistFromStart=ncol(df)) %>%  mutate(Pathway=row_number(),Start=fromGene,End=toGene) 
    df<-df %>% pivot_longer(cols=!c(DistFromStart,Pathway,Start,End)) %>% rename(distance=name,Gene=value) %>% mutate(distance=substr(distance,2,nchar(distance)))
    
    return(df)
    }
  return(NULL)
}

genesInDatabase<-as_tibble()
failedGenes<-as_tibble()
for(i in c(1:nrow(genes))){
  skip_to_next <- FALSE
  fromGene<-as.character(genes[i,1])
  tryCatch(all_shortest_paths(OPI_g,from = fromGene,to='STAT3'), error = function(e) { skip_to_next <<- TRUE}) #STAT3 is a filler
  if(!skip_to_next) {
    genesInDatabase<-rbind(fromGene,genesInDatabase)
  }
  else{
    failedGenes<-rbind(failedGenes,fromGene)
  }
}

genesInDatabase<-genesInDatabase %>% dplyr::rename(Gene= X.CD3D.)
geneVector<-as.vector(genesInDatabase$Gene)

SignalPathways<-as_tibble()
y<-Sys.time()
for (i in geneVector){
  z<-map_dfr(geneVector,~makePathsEvenEasier(.,i))
  SignalPathways<-rbind(SignalPathways,z)
}
b<-Sys.time()
b-y

names(pidList)[names(pidList) == 'Genetic defect'] <- 'Start'
#names(pidList)[names(pidList) == 'End'] <- 'Start'

newtable <- merge(SignalPathways,pidList, by  = "Start",all.x=T) 

names(newtable)[names(newtable) == 'Major category'] <- 'Start Gene Group'

names(pidList)[names(pidList) == 'Start'] <- 'End'
newtable2 <- merge(newtable,pidList, by  = "End",all.x=T) 

names(newtable2)[names(newtable2) == 'Major category'] <- 'End Gene Group'

names(newtable2)[names(newtable2) == 'Gene'] <- 'gene'

Pathways<-newtable2 %>% dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`)

Pathways<-merge(Pathways,gNOMADpLIs,by='gene',all.x=T) %>% 
  as_tibble()

i<-Pathways %>% filter(is.na(pLI))

table(i$gene)

Pathways <- Pathways %>% 
  as_tibble() %>% 
  dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`,pLI)

genes<-genes %>% mutate(AnalyzedinPathways=(ConfirmedPIDGene %in% Pathways$Start|ConfirmedPIDGene %in% Pathways$End))
failedgenes<-(genes %>% filter(AnalyzedinPathways=="FALSE"))$ConfirmedPIDGene %>% as_tibble()

correctedgenes2<-as_tibble()

for(i in c(1:nrow(failedgenes))){
  #print(failedgenes[i,1])
  x<-hgnc.table %>% as_tibble() %>% filter(Symbol==as.character(failedgenes[i,1]))
  #print(x)
  correctedgenes2<-rbind(x,correctedgenes2)
}

fixedGenes<-as_tibble()
failedgenes2<-as_tibble()
for(i in c(1:nrow(correctedgenes2))){
  skip_to_next <- FALSE
  fromGene<-as.character(correctedgenes2[i,2])
  tryCatch(all_shortest_paths(OPI_g,from = fromGene,to='STAT3'), error = function(e) { skip_to_next <<- TRUE}) #STAT3 is a filler
  if(skip_to_next) {
    failedgenes2<-rbind(fromGene,failedgenes2)
  }
  else{
    fixedGenes<-rbind(fixedGenes,fromGene)
  }
}

fixedGenes<-as.vector(fixedGenes$X.FANCG.)

SignalPathways2<-as_tibble()
for (i in geneVector){
  z<-map_dfr(fixedGenes,~makePathsEvenEasier(.,i))
  SignalPathways2<-rbind(SignalPathways2,z)
}

SignalPathways3<-as_tibble()
for (i in fixedGenes){
  z<-map_dfr(geneVector,~makePathsEvenEasier(.,i))
  SignalPathways3<-rbind(SignalPathways3,z)
}

SignalPathways4<-as_tibble()
for (i in fixedGenes){
  z<-map_dfr(fixedGenes,~makePathsEvenEasier(.,i))
  SignalPathways4<-rbind(SignalPathways4,z)
}

#colnames(SignalPathways2)[5]<-"DistFromStart"

names(pidList)[names(pidList) == 'End'] <- 'Start'

names(SignalPathways2)[names(SignalPathways2) == 'Start'] <- 'Approved.Symbol'

x<-merge(SignalPathways2,correctedgenes2, by  = "Approved.Symbol",all.x=T) 

names(x)[names(x) == 'Approved.Symbol'] <- 'Start'
names(x)[names(x) == 'Symbol'] <- 'StartSymbolChanged'

SignalPathways2<-x %>% as_tibble()

names(pidList)[names(pidList) == 'Start'] <- 'StartSymbolChanged'

newtable <- merge(SignalPathways2,pidList, by  = "StartSymbolChanged",all.x=T) 

names(newtable)[names(newtable) == 'Major category'] <- 'Start Gene Group'

names(pidList)[names(pidList) == 'StartSymbolChanged'] <- 'End'


newtable2 <- merge(newtable,pidList, by  = "End",all.x=T)

names(newtable2)[names(newtable2) == 'Major category'] <- 'End Gene Group'

names(newtable2)[names(newtable2) == 'Gene'] <- 'gene'

Pathways2<-newtable2 %>% dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`)

Pathways2<-merge(Pathways2,gNOMADpLIs,by='gene',all.x=T) %>% 
  as_tibble()

Pathways2 <- Pathways2 %>% 
  as_tibble() %>% 
  dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`,pLI)

# signal pathways 3

#colnames(SignalPathways3)[5]<-"DistFromStart"

#names(pidList)[names(pidList) == 'Genetic defect'] <- 'Start'

names(SignalPathways3)[names(SignalPathways3) == 'End'] <- 'Approved.Symbol'

x<-merge(SignalPathways3,correctedgenes2, by  = "Approved.Symbol",all.x=T) 

names(x)[names(x) == 'Approved.Symbol'] <- 'End'
names(x)[names(x) == 'Symbol'] <- 'EndSymbolChanged'
#names(x)[names(x) == 'End'] <- 'Approved.Symbol'

# x<-merge(x,correctedgenes, by  = "Approved.Symbol") 
# names(x)[names(x) == 'Approved.Symbol'] <- 'End'
# names(x)[names(x) == 'Symbol'] <- 'EndSymbolChanged'

SignalPathways3<-x %>% as_tibble()

names(pidList)[names(pidList) == 'End'] <- 'EndSymbolChanged'

newtable <- merge(SignalPathways3,pidList, by  = "EndSymbolChanged",all.x=T) 

names(newtable)[names(newtable) == 'Major category'] <- 'End Gene Group'

names(pidList)[names(pidList) == 'EndSymbolChanged'] <- 'Start'
newtable2 <- merge(newtable,pidList, by  = "Start",all.x=T)

names(newtable2)[names(newtable2) == 'Major category'] <- 'Start Gene Group'

names(newtable2)[names(newtable2) == 'Gene'] <- 'gene'

Pathways3<-newtable2 %>% dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`)

Pathways3<-merge(Pathways3,gNOMADpLIs,by='gene',all.x=T) %>% 
  as_tibble()

Pathways3 <- Pathways3 %>% 
  as_tibble() %>% 
  dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`,pLI)

# pathways 4 --------

#colnames(SignalPathways4)[5]<-"DistFromStart"
#names(pidList)[names(pidList) == 'Genetic defect'] <- 'Start'

names(SignalPathways4)[names(SignalPathways4) == 'End'] <- 'Approved.Symbol'

x<-merge(SignalPathways4,correctedgenes2, by  = "Approved.Symbol",all.x=T) 

names(x)[names(x) == 'Approved.Symbol'] <- 'End'
names(x)[names(x) == 'Symbol'] <- 'EndSymbolChanged'

names(x)[names(x) == 'Start'] <- 'Approved.Symbol'

y<-merge(x,correctedgenes2, by  = "Approved.Symbol",all.x=T) 
names(y)[names(y) == 'Symbol'] <- 'StartSymbolChanged'

SignalPathways4<-y %>% as_tibble()

names(pidList)[names(pidList) == 'Start'] <- 'EndSymbolChanged'

newtable <- merge(SignalPathways4,pidList, by  = "EndSymbolChanged",all.x=T) 

names(newtable)[names(newtable) == 'Major category'] <- 'End Gene Group'

names(pidList)[names(pidList) == 'EndSymbolChanged'] <- 'StartSymbolChanged'
newtable2 <- merge(newtable,pidList, by  = "StartSymbolChanged",all.x=T)

names(newtable2)[names(newtable2) == 'Major category'] <- 'Start Gene Group'
names(newtable2)[names(newtable2) == 'Gene'] <- 'gene'
names(newtable2)[names(newtable2) == 'Start'] <- 'StartOld'
names(newtable2)[names(newtable2) == 'End'] <- 'EndOld'

names(newtable2)[names(newtable2) == 'StartSymbolChanged'] <- 'Start'
names(newtable2)[names(newtable2) == 'EndSymbolChanged'] <- 'End'

names(newtable2)[names(newtable2) == 'Gene'] <- 'gene'

Pathways4<-newtable2 %>% dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`)

Pathways4<-merge(Pathways4,gNOMADpLIs,by='gene',all.x=T) %>% 
  as_tibble()

Pathways4 <- Pathways4 %>% 
  as_tibble() %>% 
  dplyr::select(gene,Pathway,distance,DistFromStart,Start,End,`Start Gene Group`,`End Gene Group`,pLI)

# ------

DistFiltPathways<-Pathways %>% filter(DistFromStart<6)
DistFiltPathways2<-Pathways2 %>% filter(DistFromStart<6)
DistFiltPathways3<-Pathways3 %>% filter(DistFromStart<6)
DistFiltPathways4<-Pathways4 %>% filter(DistFromStart<6)

PathwaysTotal<-bind_rows(DistFiltPathways,DistFiltPathways2,DistFiltPathways3,DistFiltPathways4)
unique(PathwaysTotal$gene)
unique(TotalInteractions$Target)

length(unique(c(unique(PathwaysTotal$gene),unique(TotalInteractions$Target))))
#  --------
pidList<-pidList %>% dplyr::rename(PIDGene=StartSymbolChanged)
TotalInteractions<-TotalInteractions %>% dplyr::mutate(InPID=PIDGene %in% pidList$PIDGene)
NameChange<-TotalInteractions %>% filter(InPID=="FALSE")

for (i in c(1:nrow(NameChange))) {
  name<-correctedgenes %>% filter(Approved.Symbol==as.character(NameChange[i,]$PIDGene))
  NameChange[i,]$PIDGene<-name$Symbol
}

TotalInteractions<-TotalInteractions %>% filter(InPID=="TRUE")
TotalInteractions<-rbind(TotalInteractions,NameChange) %>% as_tibble()

TotalInteractions<-TotalInteractions %>% dplyr::select(-Group) %>% merge(pidList,by="PIDGene",all.x=T) %>% dplyr::select(Target,PIDGene,InteractionType,
                                                                                                                 "Major category",InPID) %>% 
  as_tibble() %>% unique()

TotalInteractions<-dplyr::rename(TotalInteractions,Group=`Major category`)
pidList<-dplyr::rename(pidList,Start=PIDGene)

PathwaysTotal<-PathwaysTotal %>% mutate(InPIDStart=Start %in% pidList$Start)
PathwaysTotal<-PathwaysTotal %>% mutate(InPIDEnd=End %in% pidList$Start)

NameChangePathwaysStart<-PathwaysTotal %>% filter(InPIDStart=="FALSE")
NameChangePathwaysEnd<-PathwaysTotal %>% filter(InPIDEnd=="FALSE")

correctedgenes<-dplyr::rename(correctedgenes,Start=Approved.Symbol)

NameChangePathwaysStart<-merge(NameChangePathwaysStart,correctedgenes,by="Start",all.x=T) %>% as_tibble()
NameChangePathwaysStart<-dplyr::rename(NameChangePathwaysStart,StartBadName=Start)
NameChangePathwaysStart<-dplyr::rename(NameChangePathwaysStart,Start=Symbol)

correctedgenes<-dplyr::rename(correctedgenes,End=Start)

NameChangePathwaysEnd<-merge(NameChangePathwaysEnd,correctedgenes,by="End",all.x=T) %>% as_tibble()
NameChangePathwaysEnd<-dplyr::rename(NameChangePathwaysEnd,EndBadName=End)
NameChangePathwaysEnd<-dplyr::rename(NameChangePathwaysEnd,End=Symbol)

NameChangePathwaysStart<-dplyr::rename(NameChangePathwaysStart,BadName=StartBadName)
NameChangePathwaysEnd<-dplyr::rename(NameChangePathwaysEnd,BadName=EndBadName)

NameChangePathways<-rbind(NameChangePathwaysStart,NameChangePathwaysEnd)
NotNameChangePathways<-PathwaysTotal %>% filter(InPIDEnd!="FALSE"&InPIDStart!="FALSE")

PathwaysTotal<-rbind(NotNameChangePathways,NameChangePathways %>% dplyr::select(-BadName))

PathwaysTotal %>% mutate(InPIDStart=Start %in% pidList$Start) %>% filter(InPIDStart=="FALSE")
PathwaysTotal %>% mutate(InPIDEnd=End %in% pidList$Start) %>% filter(InPIDEnd=="FALSE")

#------
toBind1<-TotalInteractions %>% dplyr::select(Target)

names(toBind1)[names(toBind1)=="Target"]<-'gene'

toBind2<-PathwaysTotal %>% dplyr::select(gene)

boundAssociatedGenes<-unique(bind_rows(toBind1,toBind2))

boundAssociatedGenes<-boundAssociatedGenes %>% mutate(Pathway1=gene %in% TotalInteractions$Target,Pathway2=(gene %in% PathwaysTotal$gene|gene%in% PathwaysTotal$gene))

boundAssociatedGenes<-boundAssociatedGenes %>% dplyr::mutate(PIDGene=gene %in% pidList$Start)

boundAssociatedGenes %>% filter(PIDGene=="FALSE")  # 2531
boundAssociatedGenes %>% filter(PIDGene=="TRUE") 

length(unique(c(unique(TotalInteractions$PIDGene),
                unique(PathwaysTotal$Start),
                unique(PathwaysTotal$End))))
#339 PID genes analyzed

itanList<-read_xlsx("~/Downloads/ITAN_PIDList.xlsx")

boundAssociatedGenes <-boundAssociatedGenes %>% filter(PIDGene=="FALSE") %>% mutate(ItanList=gene %in% itanList$
                                                                                      `Predicted PID gene candidate`)

table(boundAssociatedGenes$ItanList) # ITAN LIST OVERLAP: 1360 NEW, 1171 OLD

boundAssociatedGenes2<-merge(boundAssociatedGenes,gNOMADpLIs,by='gene',all.x=TRUE)

novelGenes<-boundAssociatedGenes2 %>% filter(PIDGene==FALSE) %>% dplyr::select(pLI,gene,PIDGene,ItanList,pNull,pRec)

novelGenes<-novelGenes %>% as_tibble()

colnames(novelGenes)<-c("pLI","AssociatedGene","PIDGene","ItanList","pNull","pRec")

novelGenes<-novelGenes %>% dplyr::select(-c(PIDGene,ItanList))
length(unique(novelGenes$AssociatedGene))
novelGenes<-unique(novelGenes) %>% as_tibble()

#final list include PID genes

novelGenes %>% filter(AssociatedGene=="TLR9")
novelGenes %>% filter(AssociatedGene=="CKS1B")
novelGenes %>% filter(AssociatedGene=="RBL1")

novelGenes<-novelGenes[(novelGenes$pLI<.00001|novelGenes$AssociatedGene!="TLR9" ),]
novelGenes<-novelGenes[(!is.na(novelGenes$pLI)|novelGenes$AssociatedGene!="CKS1B" ),]
novelGenes<-novelGenes[(!is.na(novelGenes$pLI)|novelGenes$AssociatedGene!="LSP1" ),]
novelGenes<-novelGenes[(novelGenes$pNull<.1|novelGenes$AssociatedGene!="RBL1" ),]

novelGenes<-novelGenes %>% mutate(InItan= AssociatedGene %in% itanList$`Predicted PID gene candidate`) %>% 
  mutate(InItan=case_when(InItan=="TRUE"~"KNOWN",InItan=="FALSE"~"NEW"))

# setwd("/Users/humzakhan/Desktop/OmniPath_PID/20201208_FinalGenes")
# write.csv(novelGenes,"20201208_NovelPIDCandidateGenes_nopLI_Filter.csv")

setwd("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Lists/202105_Lists")

novelGenes<-novelGenes %>% mutate(Pathway1=AssociatedGene %in% TotalInteractions$Target,
                      Pathway2=(AssociatedGene %in% PathwaysTotal$gene|AssociatedGene%in% PathwaysTotal$gene))

novelGenes<-novelGenes %>% mutate(BothPathways=(Pathway1==TRUE&Pathway2==TRUE))

write.csv(novelGenes,"202105_NovelIEICandidateGenes_nopLI_Filter.csv")

novelGenespLI<-novelGenes %>% filter(pLI>.9)
write.csv(novelGenespLI,"202105_NovelIEICandidateGenes_HighpLI_Filter.csv")

#write.csv(novelGenespLI,"20201208_NovelPIDCandidateGenes_pLI_Filter.csv")

i<-unique(c(unique(PathwaysTotal$End),unique(PathwaysTotal$Start),unique(TotalInteractions$PIDGene)))

i<-enframe(i)

names(i)<-c("num","PIDGene") #338, good sanity check

TotalInteractions<-TotalInteractions[order(TotalInteractions$Group),]

#TotalInteractions$InteractionType <-str_replace(TotalInteractions$InteractionType," \n", "")
TotalInteractions$InteractionType <-str_replace(TotalInteractions$InteractionType,"PIDGeneIsPhosphoSubstrate", "PIDGeneIsPTMSubstrate")
TotalInteractions$InteractionType <-str_replace(TotalInteractions$InteractionType,"PIDGeneIsPhosphoKinase", "PIDGeneIsPTMEffector")
# TotalInteractions$`End Gene Group`<-str_replace(PathwaysTotal$`End Gene Group`,"\n", ": ")

write.csv(TotalInteractions,"20201208_OrderedTotalInteractions.csv")

# PathwaysTotal$`Start Gene Group`<-str_replace(PathwaysTotal$`Start Gene Group`,"\n", ": ")
# PathwaysTotal$`End Gene Group`<-str_replace(PathwaysTotal$`End Gene Group`,"\n", ": ")

PathwaysTotal<-PathwaysTotal %>% dplyr::filter(Start!=End)
PathwaysTotal<-unique(PathwaysTotal)

PathwaysTotal<-PathwaysTotal[with(PathwaysTotal, order(Start,End,Pathway,distance)),]

colnames(PathwaysTotal)<-c("Associated Gene","Pathway Number","Position","Pathway Length","Start PID Gene","End PID Gene",
                           "Start Gene Group","End Gene Group","pLI")

write.csv(PathwaysTotal[,1:8],"20201208_Pathways.csv")

# database freeze -------
phosphoInteractions
tfInteractions
ppInteractionsTib

write.csv(phosphoInteractions,"20201208_PTMData")
write.csv(tfInteractions,"20201208_TFData")
write.csv(ppInteractionsTib,"20201208_OmniPathData")


### HLA-DPB1 and HLA-DRB3 manually merged https://www.ncbi.nlm.nih.gov/gene/3125

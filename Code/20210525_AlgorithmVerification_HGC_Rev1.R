library(tidyverse)
library(resample)
library(HGNChelper)
library(RColorBrewer)

set.seed(1932)

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

cibox<-function(x){
  b<-resample::bootstrap(x,R=30000,statistic=mean);
  m = b$observed;
  data.frame(lower=CI.percentile(b)[1],
             upper=CI.percentile(b)[2],
             ymin=CI.percentile(b)[1],
             ymax=CI.percentile(b)[2],
             middle=m)
}


setwd("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Omnipath_DatabaseFreeze")

pidList<-readxl::read_xlsx("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Lists/ConfirmedPIDList/PID_Genes_FilteredForGenes_ChrRemoved.xlsx") %>% 
  dplyr::select('Genetic defect','Major category')

x<-checkGeneSymbols(pidList$`Genetic defect`) %>% as_tibble() %>% filter(!is.na(x)&Approved==F)

genes<-readxl::read_xlsx("/Users/humzakhan/Desktop/Khan_Butte_Supplement/Lists/ConfirmedPIDList/PID_Genes_FilteredForGenes_ChrRemoved.xlsx") %>% 
  dplyr::select('Genetic defect') %>% unique() %>% na.omit()

Pathways<-read_csv("20201208_Pathways.csv")
TotalInteractions<-read_csv("20201208_OrderedTotalInteractions.csv")

fractionFound<-as_tibble()
fractionFound$PercentFound<-0

totalFoundGenes<-as_tibble()
fractionFound$Gene<-""


for (i in c(1:1000)) {
  ConfirmedPIDGeneSample<-genes %>% sample_frac(size=.8) 
  ToFind<-setdiff(genes,ConfirmedPIDGeneSample)
  
  subsetPath<-Pathways %>% filter(`Start PID Gene` %in% ConfirmedPIDGeneSample$`Genetic defect`& `End PID Gene`%in% ConfirmedPIDGeneSample$`Genetic defect`)
  subsetPTMTF<-TotalInteractions %>% filter(PIDGene %in% ConfirmedPIDGeneSample$`Genetic defect`)
  
  foundgenes<-unique(c(subsetPath$`Associated Gene`,subsetPTMTF$Target)) %>% as_tibble() %>% rename(Gene=value)
  totalFoundGenes<-rbind(foundgenes,totalFoundGenes)
  
  foundgenes<-foundgenes %>% mutate(FoundPIDGene=Gene %in% ToFind$'Genetic defect')
  numFound<-foundgenes %>% count(FoundPIDGene) %>% filter(FoundPIDGene==TRUE) %>% dplyr::select(n) %>% pull()
                              
  found<-tibble(PercentFound=numFound/nrow(ToFind)*100)
  
 
  
  fractionFound<-rbind(fractionFound,found)
}

totalFoundGenes<-unique(totalFoundGenes) %>% mutate(FoundPIDGene=Gene %in% pidList$'Genetic defect')
totalFoundGenes %>% mutate()

totalFoundPIDGenes<-totalFoundGenes %>% filter(FoundPIDGene==T)

setdiff(genes %>% rename(Gene='Genetic defect'),totalFoundPIDGenes %>% select(Gene)) %>% view()

mean<-mean(fractionFound$PercentFound)

fractionFound$b<-"T"
# 
# ggplot(fractionFound,aes(x=PercentFound,fill=b))+
#   geom_density(color="#62D7FF",alpha=.1)+
#   geom_vline(xintercept=mean,linetype='dashed')+
#   scale_fill_manual(values="#62D7FF")+
#   theme_manish()
# 
# ggplot(fractionFound,aes(x=PercentFound,fill=b))+
#   geom_histogram(color="#62D7FF",alpha=.1,bins=20)+
#   geom_vline(xintercept=mean,linetype='dashed')+
#   scale_fill_manual(values="#62D7FF")+
#   theme_manish()+
#   ggtitle("")
  
ggplot(fractionFound,aes(y=PercentFound,fill=b,x=b))+
  geom_jitter(color="#62D7FF",alpha=.2)+
  geom_vline(xintercept=mean,linetype='dashed')+
  scale_fill_manual(values="#62D7FF")+
  theme_manish()+
  #stat_summary(alpha=0.1,outlier.shape = NA,fun.data=cibox, position=position_dodge(0.95),geom='boxplot') +
  geom_hline(yintercept = mean,linetype="dashed")+
  xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab("Percent Found")+
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),
        plot.title = element_text(size=5))+
  ylim(35,75)+
  ggtitle("Percent of Left-Out Known IEI Genes Rediscovered by Algorithm")

setwd("/Users/humzakhan/Desktop/OmniPath_PID/20210525_Revision1Figs")
ggsave("FigS1A_LeftOutRediscovered.png",height=50,width=60,units = "mm",dpi=3000)


# min(fractionFound$PercentFound)
# max(fractionFound$PercentFound)

# leave one out approach

fractionFound2<-as_tibble()
fractionFound2$PercentFound<-0
fractionFound2$Category<-0

for (i in c(1:9)) {
  notIncluded<-pidList %>% filter(grepl(i,`Major category`)) %>% dplyr::select('Genetic defect')
  
  included<-pidList %>% filter(!grepl(i,`Major category`)) %>% dplyr::select('Genetic defect')
  
  subsetPath<-Pathways %>% filter(`Start PID Gene` %in% included$`Genetic defect`& `End PID Gene`%in% included$`Genetic defect`)
  subsetPTMTF<-TotalInteractions %>% filter(PIDGene %in% included$`Genetic defect`)
  
  foundgenes<-unique(c(subsetPath$`Associated Gene`,subsetPTMTF$Target)) %>% as_tibble() %>% rename(Gene=value)
  
  foundgenes<-foundgenes %>% mutate(FoundPIDGene=Gene %in% notIncluded$'Genetic defect')
  numFound<-foundgenes %>% count(FoundPIDGene) %>% filter(FoundPIDGene==TRUE) %>% dplyr::select(n) %>% pull()
  
  found<-tibble(PercentFound=numFound/nrow(notIncluded)*100)  
  found$Category<-i
  fractionFound2<-rbind(fractionFound2,found)
}

fractionFound2$Category<-paste(fractionFound2$Category)
mean2<-mean(fractionFound2$PercentFound)

brewer<-c("#8DD3C7","#F6CF4B" , "#BEBADA", "#FB8072" , "#80B1D3" ,"#B3DE69","#FDB462" , "#FCCDE5" ,"#D9D9D9")

ggplot(fractionFound2,aes(x=Category,y=PercentFound,color=Category))+
  geom_point()+
  xlab("")+
  geom_hline(yintercept=mean2,linetype="dashed")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_color_manual(values=brewer)+
  ylab("Percent of Category Found")+
  theme_manish()+
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),
        plot.title = element_text(size=5))+
  xlab("Table Left Out")+
  ggtitle("Percent of Left-Out Known IEI Genes Rediscovered by Algorithm")

ggsave("FigS2A_LeftOutCategory_Rediscovered.png",height=50,width=60,units = "mm",dpi=3000)


hgcDist<-read_tsv("/Users/humzakhan/Downloads/candidates_in_core_genes.txt")

median<-median(hgcDist$Distance)

ggplot(hgcDist,aes(x=Distance))+
  geom_histogram()+
  geom_vline(xintercept = median,linetype='dashed')+
  theme_bw()

ggsave("hgcPlot.png",height=50,width=60,units = "mm",dpi=3000)


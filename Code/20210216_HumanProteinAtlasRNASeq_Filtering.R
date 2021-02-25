library(biomaRt)
library(reshape2)
library(HGNChelper)
library(tidyverse)
library(readxl)
# ------

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

List<-read_csv("20201208_NovelPIDCandidateGenes_nopLI_Filter.csv") %>% 
  dplyr::select(-X1) %>% 
  rename(Gene=AssociatedGene)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listEnsGeneName <- getBM(filters = "external_gene_name", 
                         attributes = c("ensembl_gene_id","external_gene_name"),
                         values = List$Gene, mart = mart) %>% as_tibble()

colnames(listEnsGeneName)<-c("ENSGene","Gene")

length(unique(List$Gene))
List2<-List %>% merge(listEnsGeneName,all.x=T) %>% as_tibble()
length(unique(List2$Gene))

RNAseqBloodCell<-read_tsv("~/Downloads/rna_blood_cell.tsv") %>% rename(ENSGene=Gene,Gene=`Gene name`,Cell=`Blood cell`)

#Merged<-List2 %>% merge(RNAseqBloodCell,all.x=T,by='ENSGene')

Merged<-List2 %>% merge(RNAseqBloodCell,all.x=T,by='Gene') %>% as_tibble()

MergedMinusNonCoding<-Merged %>% as_tibble() %>% filter(!is.na(TPM))

List<-List %>% mutate(RNASeq=Gene %in%MergedMinusNonCoding$Gene)

MergedMinusNonCoding<-MergedMinusNonCoding %>% mutate(TPM2=TPM+1)

# ggplot(MergedMinusNonCoding, aes(x=TPM2,fill=Cell))+
#   geom_histogram()+
#   scale_x_log10()+
#   #facet_grid(.~Cell)+
#   NULL

#length(unique(MergedMinusNonCoding$Gene))

nrow(MergedMinusNonCoding %>% filter(TPM2>2))/nrow(MergedMinusNonCoding)

RNASeqList<-MergedMinusNonCoding %>% split(list(MergedMinusNonCoding$Cell)) %>% map(function(x){x %>% dplyr::select(-c(ENSGene.x,ENSGene.y)) %>% unique()})

setwd('/Users/humzakhan/Desktop/OmniPath_PID/20210216_HumanProteomeRNASeq')

RNASeqList %>% map(function(x) x %>% dplyr::select(Gene, TPM) %>% unique())

RNASeqList %>% map(function(x){x %>% filter(TPM>1) %>% dplyr::select(Gene, TPM,Cell) %>% unique()}) %>% 
  map(function(x){write.csv(x,paste("20210216_",unique(x$Cell),"FilteredGenes.csv",sep=""))})

#RNASeqListFiltered$`naive B-cell` %>% filter(Gene=="DLEU1")

pidList<-read_xlsx("20210216_IEIListRNASeq.xlsx") %>% rename(Gene=`Genetic defect`)

genes<-read_xlsx("20210216_IEIListRNASeq.xlsx") %>% as_tibble() %>% dplyr::select('Genetic defect')
genes<-unique(genes)  
  
RNAseqBloodCell %>% filter(Gene=="TRA")
RNAseqBloodCell %>% filter(ENSGene=="ENSG00000277734") 


unique(RNAseqBloodCell$Cell)

myeloid<-c("monocyte","eosino","basophil","neutrophil")
tcell<-c("T-cell","T-reg")
bcell<-c("B-cell")

RNAseqBloodCell<-RNAseqBloodCell %>% mutate(CellType=case_when(
  grepl(paste(bcell,collapse = "|"),Cell) ~ "B cell",
  grepl(" DC",Cell) ~ "Dendritic Cell",
  grepl(paste(myeloid,collapse = "|"),Cell) ~ "Myeloid",
  grepl(paste(tcell,collapse="|"),Cell) ~ "T cell",
  grepl("NK-cell",Cell) ~ "NK cell",
  grepl("total",Cell) ~ "Total PBMC"
), TPM2=TPM+1) 


RNAseqBloodCellMeans<-RNAseqBloodCell %>% group_by(CellType, Gene) %>% summarise(mean=mean(TPM2)) 

RNAseqBloodCellMeans %>% split(list(RNAseqBloodCellMeans$CellType)) %>%
  map(function(x) {x %>% filter(mean>2)}) %>% 
  map(function(x) nrow(x)/(nrow(RNAseqBloodCellMeans)/6)*100)
  
#FIG S3
RNAseqBloodCellMeans %>% filter(CellType!="Total PBMC") %>% 
  ggplot(aes(x=mean))+
  geom_histogram(bins=100)+
  scale_x_log10(breaks=c(1,100,1000,10000),labels=c(1,100,1000,10000))+
  scale_y_log10()+ #rows removed with scale log y -- bins with 0 genes in them, doesn't matter
  facet_wrap(CellType~.,strip.position = "top",nrow =5)+
  theme_manish()+
  theme(axis.text.x = element_text(size=2.3),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),
        plot.title = element_text(size=5),strip.text.x = (element_blank()))+
  xlab("Transcripts per Million (TPM) + 1")+
  ylab("Number of Genes")+
  ggtitle("RNAseq Expression of HPA-Recorded Genes")+
  #geom_vline(xintercept = 1)+
  NULL

ggsave("AllGenes.png",height=50,width=60,units = "mm",dpi=3000)


#FIG S4
!is.na(mean)

RNAseqPID<-merge(RNAseqBloodCell,pidList,all.y=T)
(RNAseqPID %>% filter(is.na(TPM)) %>% as_tibble())$Gene
RNAseqPID %>% as_tibble()
# IGHM, TRAC, and IGKC not included... weird, ask manish about this 

#RNU4ATAC, TERC is noncoding


View(RNAseqPID %>% filter(is.na(CellType)))

RNAseqBloodCellMeans %>% split(list(RNAseqBloodCellMeans$CellType)) %>%
  map(function(x) {x %>% filter(mean>2)}) %>% 
  map(function(x) nrow(x)/(nrow(RNAseqBloodCellMeans)/6)*100)

RNAseqPIDSum<-RNAseqPID %>% as_tibble() %>% group_by(Gene,CellType) %>% summarise(mean=mean(TPM2))

RNAseqPIDSum<-left_join(RNAseqPIDSum,pidList %>% dplyr::select(Gene,`Major category`),by="Gene") %>% unique()

RNAseqPIDSum %>% split(list(RNAseqPIDSum$CellType)) %>%
  map(function(x) {x %>% filter(mean>2)}) %>% 
  map(function(x) nrow(x)/(nrow(RNAseqPIDSum)/6)*100)

unique((RNAseqPIDSum %>% filter(mean==1))$Gene)
# complement, AIRE, BRCA2, foxp3 in some non T cell types

unique((RNAseqPIDSum %>% filter(mean==1))$Gene)

RNAseqPIDSum %>% summarise(mean=mean(mean)) %>% filter(mean==1)
# 0 TPM in all immune cell types: AIRE, complement (CFHR too), ALPI: brush border enzyme for intestinal host defense, 
# CFTR, ficolin (complement activator thru lectin pathway), 
# IL17F-->didnt look at activated T cells
x


brewer<-c("#8DD3C7","#FFFFB3" , "#BEBADA", "#FB8072" , "#80B1D3" ,"#B3DE69","#FDB462" , "#FCCDE5" ,"#D9D9D9")
ggplot(RNAseqPIDSum %>% filter(CellType!='Total PBMC',!is.na(`Major category`)),aes(x=mean,fill=`Major category`))+
  geom_histogram(bins=400)+
  scale_x_log10(breaks=c(1,100,1000,10000),labels=c(1,100,1000,10000))+
  facet_wrap(CellType~.,strip.position = "top",nrow =5)+
  geom_vline(xintercept = 2)+
  theme_manish()+
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),legend.text = element_text(size=2),legend.key.size = unit(0, "mm"),
        legend.title = element_text(size=2),
        plot.title = element_text(size=5), strip.text.x = element_blank())+
  xlab("Transcripts per Million (TPM) + 1")+
  ylab("Number of Genes")+
  ggtitle("RNAseq Expression of all IEI Genes")+
  scale_fill_manual(values=brewer)+
  NULL

ggsave("AllIEI.png",height=50,width=60,units = "mm",dpi=3000)


RNAseqPIDSum %>% split(list(RNAseqPIDSum$CellType)) %>%
  map(function(x) {x %>% filter(mean>2)}) %>% 
  map(function(x) nrow(x)/(nrow(RNAseqPIDSum)/6)*100)

b<-RNAseqPIDSum %>% split(list(RNAseqPIDSum$CellType)) %>%
  map(function(x) {x %>% filter(mean<2)})

RNAseqPIDSum %>% filter(CellType=="B cell",grepl("3",`Major category`))

b$`B cell` %>% filter(grepl("3",`Major category` )) # SLC39A7, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6561116/ recently described for b cell dev

unique(b$`B cell`$`Major category`)
  
RNAseqPIDSum %>% split(list(RNAseqPIDSum$CellType)) %>%
  map(function(x) {x %>% filter(mean>2)}) %>% 
  map(function(x) nrow(x)/(nrow(RNAseqPIDSum)/6)*100)

brewer<-c("#8DD3C7")

RNAseqPIDSum %>% split(list(RNAseqPIDSum$CellType)) %>% 
  map(function(x) {x %>% filter(grepl("1",`Major category`))}) %>% 
  map(function(x) {nrow(x %>% filter(mean>2))/nrow(x) })

# percentiles

RNAseqPIDSum %>% filter(grepl("1",`Major category`), CellType!='Total PBMC') %>% 
  ggplot(aes(x=mean,fill=`Major category`))+
  geom_histogram(bins=20)+
  scale_x_log10()+
  facet_grid(.~CellType)+
  geom_vline(xintercept = 2)+
  theme_manish()+
  scale_y_continuous(breaks=c(0,2,4,6,8))+
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),legend.text = element_text(size=2),legend.key.size = unit(.33, "cm"),
        legend.title = element_text(size=2),
        plot.title = element_text(size=5), strip.text.x = element_blank())+
  xlab("Transcripts per Million (TPM) + 1")+
  ylab("Number of Genes")+
  ggtitle("RNAseq Expression of Table 1 (SCID+CID) Genes")+
  scale_fill_manual(values=brewer)+
  NULL

ggsave("SCIDCIDGenes.png",height=50,width=60,units = "mm",dpi=3000)


brewer<-c("#FFFFB3")
RNAseqPIDSum %>% filter(grepl("3",`Major category`), CellType=='B cell') %>% 
  ggplot(aes(x=mean,fill=`Major category`))+
  geom_histogram(bins=34)+
  scale_x_log10()+
  facet_grid(.~CellType)+
  geom_vline(xintercept = 2)+
  theme_manish()+
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),legend.text = element_text(size=2),legend.key.size = unit(.33, "cm"),
        legend.title = element_text(size=2),
        plot.title = element_text(size=5), strip.text.x = element_blank())+
  xlab("Transcripts per Million (TPM) + 1")+
  ylab("Number of Genes")+
  ggtitle("RNAseq Expression of Table 1 (SCID+CID) Genes")+
  scale_fill_manual(values=brewer)+
  NULL

ggsave("PADGenes.png",height=50,width=60,units = "mm",dpi=3000)

  
# length(unique(List$Gene))
#  length(List$Gene)
#  
# data.frame(table(List$Gene)) %>% filter(Freq!=1)
# List %>% filter(Gene=="HLA-DPB1")
# List %>% filter(Gene=="RBL1")
 
 

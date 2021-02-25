library(biomaRt)
library(reshape2)
library(HGNChelper)

# ------

List<-read_csv("20201208_NovelPIDCandidateGenes_nopLI_Filter.csv") %>% select(-X1) %>% 
  rename(Gene=AssociatedGene)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listEnsGeneName <- getBM(filters = "external_gene_name", 
                         attributes = c("ensembl_gene_id","external_gene_name"),
                         values = List$Gene, mart = mart) %>% as_tibble()

colnames(listEnsGeneName)<-c("ENSGene","Gene")

length(unique(List$Gene))
List2<-List %>% merge(listEnsGeneName,all.x=T) %>% as_tibble()
length(unique(List2$Gene))

RNAseqBloodCell<-read_tsv("~/Downloads/rna_blood_cell.tsv") %>% rename(ENSGene=Gene,Gene=`Gene name`)

#Merged<-List2 %>% merge(RNAseqBloodCell,all.x=T,by='ENSGene')

Merged<-List2 %>% merge(RNAseqBloodCell,all.x=T,by='Gene')

RNAseqBloodCell %>% filter(ENSGene=="ENSG00000278677")

NAs<-Merged %>% mutate(TPMThere=case_when(
  is.na(TPM) ~ "NotAnalyzed"
))

not<-NAs %>% filter(TPMThere=="NotAnalyzed")
length(unique((NAs %>% filter(TPMThere=="NotAnalyzed"))$Gene))
# 
# failedgenes<-unique(not$Gene)
# correctedgenes<-as_tibble()
# # 
# # for(i in c(1:length(failedgenes))){
# #   #print(failedgenes[i,1])
# #   x<-hgnc.table %>% as_tibble() %>% filter(Symbol==as.character(failedgenes[i]))
# #   #print(x)
# #   correctedgenes<-rbind(x,correctedgenes)
# # }

not %>% filter(Gene=="HLA-DQA1")
failedgenes %>% enframe() %>% filter(value=="HLA-DQA1")
RNAseqBloodCell %>% filter(Gene=="HLA-DQA1")

# many failed genes just have ENS duplicates (diff haplotypes, expressed differently) ie:

not %>% filter(Gene=="NXN")
failedgenes %>% enframe() %>% filter(value=="NXN")
RNAseqBloodCell %>% filter(Gene=="NXN")

# can just drop na, most will be there no problem

MergedPresent<-Merged %>% as_tibble() %>% filter(!is.na(TPM))

MergedPresent %>% filter(Gene=="ACPP")

Merged %>% filter(Gene=="ACPP")

List<-List %>% mutate(RNASeq=Gene %in%MergedPresent$Gene)

length(unique(Merged$Gene))
length(unique(MergedPresent$Gene)) #lost 13

List %>% filter(RNASeq==FALSE)

#hand code these last ones--based on gene cards 

RNAseqBloodCell %>% filter(Gene=="MARCHF1") # noncoding, probably did polyA pulldown for RNAseq experiments

unique((RNAseqBloodCell %>% filter(grepl("GPX",Gene)))$Gene)  # no GPX1 -- ok 

unique((RNAseqBloodCell %>% filter(grepl("Opo",Gene)))$Gene)  

RNAseqBloodCell %>% filter(Gene=="DRB4") # no HLA-DRB4

# DLEU1, LINC00299, LINC00474, MDS2: noncoding RNAs. Keep in interactions, but not in RNAseq

# MIR17HG: miRNA

#CEMIP2 -> TMEM2
# DIPK2B -> CXorf36
# HLA-DRB3 -> HLA-DPB1 
# TENT4A -> PAPD7
# GPX1, OFCC1, and HLA-DRB4 protein coding but not found





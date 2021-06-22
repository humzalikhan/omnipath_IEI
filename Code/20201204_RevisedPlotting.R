library(tidyverse)
library(readxl)
library(rmutil)
library(biomaRt)
library(HGNChelper)
library(OmnipathR)
library(RColorBrewer)


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

gNOMADpLIs<-read_tsv("/pLI_GDI_Data/gNOMAD_pLI.tsv")

# import data -------
pidList<-read_xlsx("PID_Genes_FilteredForGenes_ChrRemoved.xlsx")

genes<-pidList %>% dplyr::select('Genetic defect')
genes<-unique(genes)

colnames(genes)<-"ConfirmedPIDGene"

# DICOVERY PATHWAY 1 PLOTTING ------
TotalInteractions$InteractionType<-as.character(TotalInteractions$InteractionType)
#TotalInteractions[TotalInteractions$InteractionType=="Complex",]$InteractionType<-"Complex"
TotalInteractions[TotalInteractions$InteractionType=="PID GeneisPTMEffector",]$InteractionType<-"PID Gene is PTM Effector"
TotalInteractions[TotalInteractions$InteractionType=="PID GeneisPTMSubstrate",]$InteractionType<-"PID Gene is PTM \n Substrate"
TotalInteractions[TotalInteractions$InteractionType=="PIDGeneIsTF",]$InteractionType<- "PID Gene is TF"
TotalInteractions[TotalInteractions$InteractionType=="PIDGeneIsTranscribed",]$InteractionType<- "PID Gene is Transcribed"
#
# HighpLI[HighpLI$InteractionType=="Complex",]$InteractionType<-"Gene Makes Complex with PID Gene"
# HighpLI[HighpLI$InteractionType=="PID Gene Modifies Gene",]$InteractionType<-"PID Gene PT Modifies Gene"
# HighpLI[HighpLI$InteractionType=="PID Gene is PTM Substrate",]$InteractionType<-"Gene PT Modifies PID Gene"
# HighpLI[HighpLI$InteractionType=="PID Gene isTF",]$InteractionType<- "PID Gene is TF"
# HighpLI[HighpLI$InteractionType=="PID Gene is Transcribed",]$InteractionType<- "Gene Is TF to PID Gene"

x<-table(TotalInteractions$InteractionType) %>% as.data.frame() %>% as_tibble()

TotalInteractions <- within(TotalInteractions, 
                            InteractionType <- factor(InteractionType, 
                                                      levels=names(sort(table(InteractionType), 
                                                                        decreasing=TRUE))))

brewer <- brewer.pal(4, "Set3")
brewer[2]<-"#F6CF4B"

setwd("/Users/humzakhan/Desktop/OmniPath_PID/20210525_Revision1Figs")

ggplot(TotalInteractions,aes(x=InteractionType,fill=InteractionType))+
  geom_bar()+
  ggtitle("Types of Interactions with IEI Gene")+
  theme_manish()+
  ylab("Associated Genes (repeats included)")+
  xlab("Interaction Type")+
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),
        plot.title = element_text(size=5))+
  scale_fill_manual(values=brewer)+
  NULL

ggsave("Fig1a_InteractionTypeCounts_TotalInteractions.png",height=50,width=60,units = "mm",dpi=3000)

TotalInteractions[grepl("Table 1",TotalInteractions$Group),]$Group<-"T1\nCellular and Humoral"
TotalInteractions[grepl("Table 3",TotalInteractions$Group),]$Group<-"T3\nAb Deficiencies"
TotalInteractions[grepl("Table 7",TotalInteractions$Group),]$Group<-"T7\nAutoinflammatory"
TotalInteractions[grepl("Table 4",TotalInteractions$Group),]$Group<-"T4\nImmune Dysregulation"
TotalInteractions[grepl("Table 2",TotalInteractions$Group),]$Group<-"T2\nCombined immunodeficiencies"
TotalInteractions[grepl("Table 9",TotalInteractions$Group),]$Group<-"T9\nBM Failure"
TotalInteractions[grepl("Table 5",TotalInteractions$Group),]$Group<-"T5\nPhagocyte Defects"
TotalInteractions[grepl("Table 6",TotalInteractions$Group),]$Group<-"T6\nInnate Immunity Defects"
TotalInteractions[grepl("Table 8",TotalInteractions$Group),]$Group<-"T8\nComplement Defects"

brewer<-c("#8DD3C7","#F6CF4B" , "#BEBADA", "#FB8072" , "#80B1D3" ,"#B3DE69","#FDB462" , "#FCCDE5" ,"#D9D9D9")

TotalInteractions <- within(TotalInteractions, 
                            InteractionType <- factor(InteractionType, 
                                                      levels=names(sort(table(InteractionType), 
                                                                        decreasing=TRUE))))

tableCountsog<-table(TotalInteractions$Group) %>% as.data.frame() %>% as_tibble()

#$Var1<-reorder(tableCountsog$Var1,tableCountsog$Freq) %>% fct_rev()

ggplot(tableCountsog,aes(x=Var1,y=Freq,fill=Var1))+
  geom_col()+
  ggtitle("IEI-related Gene Group")+
  theme_manish_legend("right")+
  ylab("Number of Associated Genes (repeats included)")+
  xlab("")+
  theme(axis.text.x = element_text(size=1),axis.title.y = element_text(size=5),plot.title = element_text(size=5),axis.ticks.x = element_blank(),
        legend.title = element_blank(),legend.text = element_text(size=3),legend.key.size = unit(3, "mm"),legend.key.width = unit(3, "mm"))+
  scale_fill_manual(values=brewer)

ggsave("Discovery1_MajorGroupCounts_TotalInteractions.png",height=50,width=60,units = "mm",dpi=3000)

brewer <- c("#F6CF4B","#BEBADA","#8DD3C7","#FB8072")	
TotalInt2<-TotalInteractions %>% dplyr::select(PIDGene,InteractionType) %>% unique()	
TotalInt2 <- within(TotalInt2, 	
                    InteractionType <- factor(InteractionType, 	
                                              levels=names(sort(table(InteractionType), 	
                                                                decreasing=TRUE))))	
tableCountsPIDTxn<-table(TotalInt2$InteractionType) %>% as.data.frame() %>% as_tibble()	
ggplot(TotalInt2,aes(x=InteractionType,fill=InteractionType))+	
  geom_bar()+	
  ggtitle("IEI-related Gene Group")+	
  theme_manish()+	
  ylab("Number of IEI Genes")+	
  xlab("")+	
  theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=5),
        axis.title.x = element_text(size=5),
        plot.title = element_text(size=5))+
  scale_fill_manual(values=brewer)	

ggsave("IEIGene_Txn.png",height=50,width=60,units = "mm",dpi=3000)

pidList[grepl("Table 1",pidList$`Major category`),]$`Major category`<-"T1\nCellular and Humoral"
pidList[grepl("Table 3",pidList$`Major category`),]$`Major category`<-"T3\nAb Deficiencies"
pidList[grepl("Table 7",pidList$`Major category`),]$`Major category`<-"T7\nAutoinflammatory"
pidList[grepl("Table 4",pidList$`Major category`),]$`Major category`<-"T4\nImmune Dysregulation"
pidList[grepl("Table 2",pidList$`Major category`),]$`Major category`<-"T2\nCombined immunodeficiencies"
pidList[grepl("Table 9",pidList$`Major category`),]$`Major category`<-"T9\nBM Failure"
pidList[grepl("Table 5",pidList$`Major category`),]$`Major category`<-"T5\nPhagocyte Defects"
pidList[grepl("Table 6",pidList$`Major category`),]$`Major category`<-"T6\nInnate Immunity Defects"
pidList[grepl("Table 8",pidList$`Major category`),]$`Major category`<-"T8\nComplement Defects"

# pidList1 <- within(pidList, 
#                    X <- factor('Major category', 
#                                levels=names(sort(table(pidList$`Major category`), 
#                                                  decreasing=TRUE))))

tableCounts<-table(pidList$`Major category`) %>% as.data.frame() %>% as_tibble()

#tableCounts$Var1<-fct_relevel(tableCounts$Var1, levels(tableCountsog$Var1))

#tableCounts$Var1<-reorder(tableCounts$Var1,tableCounts$Freq) %>% fct_rev()

#brewer<-c("#80B1D3","#8DD3C7", "#FB8072", "#B3DE69" ,"#FDB462" ,"#FCCDE5", "#FFFFB3" ,"#BEBADA" ,"#D9D9D9")
brewer<-c("#8DD3C7","#F6CF4B" , "#BEBADA", "#FB8072" , "#80B1D3" ,"#B3DE69","#FDB462" , "#FCCDE5" ,"#D9D9D9")
ggplot(tableCounts,aes(x=Var1,y=Freq,fill=Var1))+
  geom_col()+
  ggtitle("Confirmed IEI Gene Group")+
  theme_manish_legend("right")+
  ylab("Number of Genes in Group")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size=5),plot.title = element_text(size=5),axis.ticks.x = element_blank(),
        legend.title = element_blank(),legend.text = element_text(size=3),legend.key.size = unit(3, "mm"),legend.key.width = unit(3, "mm"))+
  scale_fill_manual(values=brewer)

ggsave("Fig2aConfirmedPID_GroupCounts.png",height=50,width=60,units = "mm",dpi=3000)

PathwaysTotal<-PathwaysTotal %>% filter(`Start PID Gene`!=`End PID Gene`)

PathwaysTotal[grepl("Table 1",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T1: Cellular and Humoral"
PathwaysTotal[grepl("Table 3",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T3: Ab Deficiencies"
PathwaysTotal[grepl("Table 7",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T7: Autoinflammatory"
PathwaysTotal[grepl("Table 4",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T4: Immune Dysregulation"
PathwaysTotal[grepl("Table 2",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T2: Combined immunodeficiencies"
PathwaysTotal[grepl("Table 9",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T9: BM Failure"
PathwaysTotal[grepl("Table 5",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T5: Phagocyte Defects"
PathwaysTotal[grepl("Table 6",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T6: Innate Immunity Defects"
PathwaysTotal[grepl("Table 8",PathwaysTotal$`Start Gene Group`),]$`Start Gene Group`<-"T8: Complement Defects"

PathwaysTotal[grepl("Table 1",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T1: Cellular and Humoral"
PathwaysTotal[grepl("Table 3",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T3: Ab Deficiencies"
PathwaysTotal[grepl("Table 7",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T7: Autoinflammatory"
PathwaysTotal[grepl("Table 4",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T4: Immune Dysregulation"
PathwaysTotal[grepl("Table 2",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T2: Combined immunodeficiencies"
PathwaysTotal[grepl("Table 9",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T9: BM Failure"
PathwaysTotal[grepl("Table 5",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T5: Phagocyte Defects"
PathwaysTotal[grepl("Table 6",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T6: Innate Immunity Defects"
PathwaysTotal[grepl("Table 8",PathwaysTotal$`End Gene Group`),]$`End Gene Group`<-"T8: Complement Defects"

PathwaysTotal<-PathwaysTotal %>% drop_na()

PathwaysTotal$SumGroup<-paste(PathwaysTotal$`Start Gene Group`,PathwaysTotal$`End Gene Group`)

groupCounts<-table(PathwaysTotal$SumGroup) %>% as.data.frame() %>% as_tibble()

groupCounts<-separate(groupCounts,Var1,into=c("StartGroup","EndGroup"),sep=" T",remove=F)

SplitPathways<-PathwaysTotal %>% dplyr::select(-c(`Pathway Number` ,Position, `Pathway Length`)) %>% 
  split(list(PathwaysTotal$`Start PID Gene` ,PathwaysTotal$`End PID Gene` ))

#SplitPathways<-PathwaysTotal %>% dplyr::select(-c(`Pathway Length`,`Pathway Number`,Position)) %>% split(list(PathwaysTotal$`Start PID Gene`,PathwaysTotal$`End PID Gene`))

x<-map_dfr(SplitPathways,~unique(.))

x<-x %>% filter(`Associated Gene`!=`Start PID Gene`&`Associated Gene`!=`End PID Gene`)

CountofGenes<-table(x$`Associated Gene`) %>% as.data.frame() %>% as_tibble() # somehow this seems quite important 

CountofGenes$Var1<-as.factor(CountofGenes$Var1)

CountofGenes<-CountofGenes %>% dplyr::mutate(ConfirmedPID=Var1 %in% genes$ConfirmedPIDGene)

CountofGenes$Var1 <- reorder(CountofGenes$Var1,CountofGenes$Freq)

colors <- c("#fabab5","#9dd1d2")

i<-CountofGenes %>% filter(ConfirmedPID=="TRUE")

i[order(i$Freq),]
i[order(i$Freq,decreasing = T),]

png("test.png", height=50,width=60,units="mm",res=3000)

ggplot(CountofGenes,aes(x=Var1,y=Freq,color=ConfirmedPID))+
  geom_point()+
  ggtitle("Gene Frequency in IEI Gene-IEI Gene Pathways")+
  theme_manish_legend(legend="right")+
  scale_color_manual(values = colors,name = "Validated IEI Gene?", labels = c("No","Yes"))+
  theme(axis.text.x=element_text(),axis.text.y=element_text(size=0))+
  theme(axis.text.x = element_text(size=3),axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=5),
        axis.title.y = element_text(size=5),
        plot.title = element_text(size=5),
        legend.title = element_text(size=3),
        legend.text = element_text(size=2.5))+
  ylab("Count")+
  xlab("Gene")+
  coord_flip()
dev.off()

CountofGenes[order(CountofGenes$Freq,decreasing = T),]

CountofGenes %>% dplyr::filter(Var1=="NFKB1")

table(CountofGenes$ConfirmedPID)

ggsave("Fig3cDistanceFiltered_UNIQUEPERPATHWAY_FrequencyGeneAppearance_Unfiltered_Path2.png",height=50,width=60,units="mm")

# circle plots -----

PIDGenesInPathways<-x %>% mutate(PIDGene=`Associated Gene` %in% genes$ConfirmedPIDGene) %>% filter(PIDGene==TRUE)

NonPIDGenesInPathways<-x %>% mutate(PIDGene=`Associated Gene` %in% genes$ConfirmedPIDGene) %>% filter(PIDGene==FALSE)

groupCountsunique<-table(NonPIDGenesInPathways$SumGroup) %>% as.data.frame() %>% as_tibble()
groupCountsunique<-separate(groupCountsunique,Var1,into=c("StartGroup","EndGroup"),sep="(?<= T)",remove=F)

groupCountsunique[grepl("1",groupCountsunique$StartGroup),]$StartGroup<-"T1"#: Cellular and Humoral"
groupCountsunique[grepl("3",groupCountsunique$StartGroup),]$StartGroup<-"T3"#: Ab Deficiencies"
groupCountsunique[grepl("7",groupCountsunique$StartGroup),]$StartGroup<-"T7"#: Autoinflammatory"
groupCountsunique[grepl("4",groupCountsunique$StartGroup),]$StartGroup<-"T4"#: Immune Dysregulation"
groupCountsunique[grepl("2",groupCountsunique$StartGroup),]$StartGroup<-"T2"#: Combined immunodeficiencies"
groupCountsunique[grepl("9",groupCountsunique$StartGroup),]$StartGroup<-"T9"#: BM Failure"
groupCountsunique[grepl("5",groupCountsunique$StartGroup),]$StartGroup<-"T5"#: Phagocyte Defects"
groupCountsunique[grepl("6",groupCountsunique$StartGroup),]$StartGroup<-"T6"#: Innate Immunity Defects"
groupCountsunique[grepl("8",groupCountsunique$StartGroup),]$StartGroup<-"T8"#: Complement Defects"

groupCountsunique[grepl("1",groupCountsunique$EndGroup),]$EndGroup<-"T1: Cellular and Humoral"
groupCountsunique[grepl("3",groupCountsunique$EndGroup),]$EndGroup<-"T3: Ab Deficiencies"
groupCountsunique[grepl("7",groupCountsunique$EndGroup),]$EndGroup<-"T7: Autoinflammatory"
groupCountsunique[grepl("4",groupCountsunique$EndGroup),]$EndGroup<-"T4: Immune Dysregulation"
groupCountsunique[grepl("2",groupCountsunique$EndGroup),]$EndGroup<-"T2: Combined immunodeficiencies"
groupCountsunique[grepl("9",groupCountsunique$EndGroup),]$EndGroup<-"T9: BM Failure"
groupCountsunique[grepl("5",groupCountsunique$EndGroup),]$EndGroup<-"T5: Phagocyte Defects"
groupCountsunique[grepl("6",groupCountsunique$EndGroup),]$EndGroup<-"T6: Innate Immunity Defects"
groupCountsunique[grepl("8",groupCountsunique$EndGroup),]$EndGroup<-"T8: Complement Defects"

groupCountsunique2<-table(PIDGenesInPathways$SumGroup) %>% as.data.frame() %>% as_tibble()
groupCountsunique2<-separate(groupCountsunique2,Var1,into=c("StartGroup","EndGroup"),sep="(?<= T)",remove=F)

groupCountsunique2[grepl("1",groupCountsunique2$StartGroup),]$StartGroup<-"T1"#: Cellular and Humoral"
groupCountsunique2[grepl("3",groupCountsunique2$StartGroup),]$StartGroup<-"T3"#: Ab Deficiencies"
groupCountsunique2[grepl("7",groupCountsunique2$StartGroup),]$StartGroup<-"T7"#: Autoinflammatory"
groupCountsunique2[grepl("4",groupCountsunique2$StartGroup),]$StartGroup<-"T4"#: Immune Dysregulation"
groupCountsunique2[grepl("2",groupCountsunique2$StartGroup),]$StartGroup<-"T2"#: Combined immunodeficiencies"
groupCountsunique2[grepl("9",groupCountsunique2$StartGroup),]$StartGroup<-"T9"#: BM Failure"
groupCountsunique2[grepl("5",groupCountsunique2$StartGroup),]$StartGroup<-"T5"#: Phagocyte Defects"
groupCountsunique2[grepl("6",groupCountsunique2$StartGroup),]$StartGroup<-"T6"#: Innate Immunity Defects"
groupCountsunique2[grepl("8",groupCountsunique2$StartGroup),]$StartGroup<-"T8"#: Complement Defects"

groupCountsunique2[grepl("1",groupCountsunique2$EndGroup),]$EndGroup<-"T1: Cellular and Humoral"
groupCountsunique2[grepl("3",groupCountsunique2$EndGroup),]$EndGroup<-"T3: Ab Deficiencies"
groupCountsunique2[grepl("7",groupCountsunique2$EndGroup),]$EndGroup<-"T7: Autoinflammatory"
groupCountsunique2[grepl("4",groupCountsunique2$EndGroup),]$EndGroup<-"T4: Immune Dysregulation"
groupCountsunique2[grepl("2",groupCountsunique2$EndGroup),]$EndGroup<-"T2: Combined immunodeficiencies"
groupCountsunique2[grepl("9",groupCountsunique2$EndGroup),]$EndGroup<-"T9: BM Failure"
groupCountsunique2[grepl("5",groupCountsunique2$EndGroup),]$EndGroup<-"T5: Phagocyte Defects"
groupCountsunique2[grepl("6",groupCountsunique2$EndGroup),]$EndGroup<-"T6: Innate Immunity Defects"
groupCountsunique2[grepl("8",groupCountsunique2$EndGroup),]$EndGroup<-"T8: Complement Defects"

# Without unique per-pathway genes. A gene will be doubly counted if it occcurs in distinct start-end pathways 
ggplot(groupCountsunique,aes(x=StartGroup,y=EndGroup,color=StartGroup,size=Freq))+
  geom_point(alpha=.6)+
  scale_color_manual(values=brewer)+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  theme_manish_legend('right')+
  xlab("Starting Gene Group")+
  ylab("Ending Gene Group")+
  scale_size_continuous(range=c(.5,3))+
  theme(axis.title.x = element_text(size=5),axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5),axis.title.y = element_text(size=5),
        plot.title = element_text(size=5,hjust=1))+
  ggtitle("Candidate IEI Genes on Inter-IEI-Gene Pathways")

ggsave("CircleGraph_Pathway2_NonPIDGenes.png",height=50,width=60,units="mm",dpi=3000)

#making circles ~ 1/3 as small bc only 1/3 as large 

ggplot(groupCountsunique2,aes(x=StartGroup,y=EndGroup,color=StartGroup,size=Freq))+
  geom_point(alpha=.6)+
  scale_color_manual(values=brewer)+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  theme_manish_legend('right')+
  xlab("Starting Gene Group")+
  ylab("Ending Gene Group")+
  scale_size_continuous(range=c(.2,1))+
  theme(axis.title.x = element_text(size=5),axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5),axis.title.y = element_text(size=5),
        plot.title = element_text(size=5,hjust=1))+
  ggtitle("Known IEI Genes on Inter-IEI-Gene Pathways")

ggsave("CircleGraph_Pathway2_PIDGenes.png",height=50,width=60,units="mm",dpi=3000)

# 
# #write.csv(TotalInteractions %>% dplyr::select(-Group),"20200819_Type1DiscoveryNovelPIDGenes.csv")
# 
# PathwaysTotalToWrite<-PathwaysTotal %>% dplyr::select( `Associated Gene`,`Pathway Number`,Position, `Pathway Length`,
#                                                 `Start PID Gene`,`End PID Gene`)
# PathwaysTotalToWrite[order(PathwaysTotalToWrite$`Pathway Number`),]
# 
# #write.csv(,"20200819_Type1DiscoveryNovelPIDGenes.csv")

# pLI Supp -----


# -----
CombinedPIDGeneAsKinase<-rbind(pidGeneAsKinase,pidGeneAsKinase3)
x<-CombinedPIDGeneAsKinase %>% as_tibble() %>% dplyr::count(modification)

x$modification<-x$modification %>% fct_relevel(c("phosphorylation","dephosphorylation","cleavage","ubiquitination"))

brewer<-c("#F8766D","#CD9600","#7CAE00","#00C19F")

ggplot(x,aes(x=modification,y=n,fill=modification))+
  geom_col()+
  theme_manish_legend("right")+
 # scale_y_log10()+
  xlab("Modification")+
  ylab("Number of genes (repeats included)")+
  scale_fill_manual(values=brewer)+
  theme(axis.title.x = element_text(size=5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5),axis.title.y = element_text(size=5),
        plot.title = element_text(size=5),
        legend.title = element_blank(),legend.text = element_text(size=3),legend.key.size = unit(3, "mm"),legend.key.width = unit(3, "mm"))+
  ggtitle("Type of PTM by PID Gene")

ggsave("Supp1aPTMbyPIDGene.png",height=50,width=60,units="mm",dpi=3000)

CombinedPIDGeneAsSubstrate<-rbind(pidGeneAsSubstrate,pidGeneAsSubstrate3)

x<-CombinedPIDGeneAsSubstrate %>% as_tibble() %>% dplyr::count(modification)

x$modification<-reorder(x$modification,x$n) %>% fct_rev()	
x <- within(x,modification <- factor(modification, 	
                                     levels=names(sort(table(modification), 	
                                                       decreasing=TRUE))))	


ggplot(x,aes(x=modification,y=n,fill=modification))+
  geom_col()+
  theme_manish_legend("right")+
 # scale_y_log10()+
  xlab("Modification")+
  ylab("Number of genes (repeats included)")+
 # scale_x_discrete(limits = rev(levels(x$modification)))+
  theme(axis.title.x = element_text(size=5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5),axis.title.y = element_text(size=5),
        plot.title = element_text(size=5),
        legend.title = element_blank(),legend.text = element_text(size=3),legend.key.size = unit(3, "mm"),legend.key.width = unit(3, "mm"))+
  ggtitle("Type of PTM Received by PID Gene")

ggsave("Supp1bPIDGenePTM.png",height=50,width=60,units="mm",dpi=3000)


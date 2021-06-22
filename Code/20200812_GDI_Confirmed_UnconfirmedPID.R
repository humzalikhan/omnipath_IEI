library(tidyverse)
library(readxl)
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
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_line()) +
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

gNOMADpLIs<-read_tsv("~/Desktop/OmniPath_PID/pLI_GDI_Data/gNOMAD_pLI.tsv")

genes<-read_xlsx("~/Downloads/PID_Genes_FilteredForGenes_ChrRemoved.xlsx") %>% as_tibble() %>% dplyr::select(`Genetic defect`,Inheritance)
genes<-unique(genes)

genes<-genes %>% filter(!is.na(`Genetic defect`))

names(genes)[names(genes)=="Genetic defect"]<-"gene"

joined<-inner_join(genes,gNOMADpLIs,by="gene",all.x=T) 

notFound<-setdiff(genes$gene,joined$gene)

joined<-joined %>% mutate(InheritanceSimplified=case_when(
  grepl("AR",Inheritance) ~ "AR",
  grepl("AD",Inheritance) ~ "AD",
  grepl("X",Inheritance) ~ "X-Linked"
))

#setwd("~/Desktop/OmniPath_PID/20201206_FinalFigs/")

setwd("~/Desktop/OmniPath_PID/20210525_Revision1Figs/")

ggplot(joined %>% filter(InheritanceSimplified!="NA"),aes(x=pLI))+
  geom_histogram()+
  theme_manish()+
  ggtitle("Confirmed IEI Gene pLIs")+
  xlab("pLI")+
  ylab("Number of genes")+
  facet_grid(InheritanceSimplified~.,scales = "free_y")+
  theme(strip.text.y = element_text(size=4,angle = 0),axis.text.x = element_text(size=4),axis.text.y = element_text(size=4),
        axis.title.y = element_text(size=6),
        axis.title.x = element_text(size=6),title = element_text(size=6))
# 
# theme(axis.text.x = element_text(size=2.5),axis.text.y = element_text(size=3),axis.title.y = element_text(size=4),
#  axis.title.x = element_text(size=4),
#       plot.title = element_text(size=5))+
#   scale_fill_manual(values=brewer)+
#   NULL

ggsave("ConfirmedPIDGene_pLI.png",height=50,width=60,units = "mm",dpi=3000)


GDIs<-read_tsv("/Users/humzakhan/Desktop/OmniPath_PID/pLI_GDI_Data/20200812_GeneDamageIndex_ConfirmedPIDGenes.tsv")
names(GDIs)[names(GDIs)=="Gene"]<-"gene"
GDIs<-inner_join(GDIs,genes,by="gene") 

GDIs<-GDIs %>% mutate(InheritanceSimplified=case_when(
  grepl("AR",Inheritance) ~ "AR",
  grepl("AD",Inheritance) ~ "AD",
  grepl("X",Inheritance) ~ "X-Linked"
))

GDIs<-unique(GDIs)

brewer<-c("#F27465","#F2D863","#70F263")

GDIs$GDI_Prediction<-GDIs$GDI_Prediction %>% fct_relevel(c("Low","Medium","High"))

ggplot(GDIs  %>% filter(InheritanceSimplified!="NA",GDI_Score!="NA"),aes(x=GDI_Score,fill=GDI_Prediction))+
  geom_histogram(bins=50)+
  theme_manish_legend("right")+
  ggtitle("Confirmed IEI Gene GDI Scores")+
  xlab("GDI")+
  ylab("Number of genes")+
  facet_grid(InheritanceSimplified~.,scales = "free_y")+
  theme(strip.text.y = element_text(angle = 0))+
  theme(strip.text.y = element_text(size=3,angle = 0),axis.text.x = element_text(size=4),axis.text.y = element_text(size=4),
        axis.title.y = element_text(size=6),
        axis.title.x = element_text(size=5),legend.text = element_text(size=3),legend.key.size = unit(.3, "cm"),
        legend.title = element_text(size=3),
        plot.title = element_text(size=5),strip.text.x = (element_text(size=1,margin = margin(.05, 0, .05, 0, "cm"))))+
  scale_fill_manual(values = brewer)

ggsave("ConfirmedPIDGene_GDI.png",height=50,width=60,units = "mm",dpi=3000)

GDIs %>% filter(GDI_Score>5000)


ggplot(GDIs  %>% filter(GDI_Prediction!="NA"),aes(x=GDI_Prediction,fill=GDI_Prediction))+
  geom_bar()+
  theme_manish()+
  ggtitle("Confirmed IEI Gene GDI Scores")+
  xlab("Predicted Mutation Effect")+
  ylab("Number of genes")+
  theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4),
        axis.title.y = element_text(size=6),
        axis.title.x = element_text(size=6),title = element_text(size=6))+
  scale_fill_manual(values = brewer)

ggsave("Supp2c_ConfirmedPIDGene_GDI_Prediction.png",height=50,width=60,units = "mm",dpi=3000)

GDIs_Unconf<-read_tsv("~/Desktop/OmniPath_PID/pLI_GDI_Data/20201208_GDI_UnconfirmedIEIgenes.tsv") %>% rename(gene=Gene)

# ggplot(GDIs_Unconf,aes(x=GDI_Score))+
#   geom_histogram()+
#   theme_manish()

GDIs %>% filter(GDI_Score>5000)

names(GDIs)[names(GDIs)=="End"]<-"gene"

GDIs<-inner_join(GDIs,gNOMADpLIs,by="gene") 

GDIs$GDI_Score
GDIs$pLI

ggplot(GDIs %>% filter(InheritanceSimplified!="NA"),aes(x=GDI_Score,y=pLI))+
  geom_point()+
  theme_manish()+
  ggtitle("Confirmed PID Gene GDI Scores")+
  xlab("GDI Score")+
  facet_grid(InheritanceSimplified~.,scales = "free_y")+
  theme(strip.text.y = element_text(angle = 0))+
 # xlim(0,5000)+
  theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4),
        axis.title.y = element_text(size=6),
        axis.title.x = element_text(size=6),title = element_text(size=6))

ggsave("Supp2dConfirmedPIDGene_GDIbypLI.png",height=50,width=60,units = "mm",dpi=3000)


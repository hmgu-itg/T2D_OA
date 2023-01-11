library(data.table)
library(ggplot2)
library(RColorBrewer)

# 1=OA, 2=T2D, 3=both
hc.genes.dt <- data.table(genes=as.list(fread("/project_data/overlap_T2D_OA/hc_genes.txt", header=FALSE))$V1,
                          mQTL=c(NA,1,NA,NA,NA,NA,NA,2,2,1,NA,NA,NA,NA,NA,3,NA,3,1),
                          DEG=c(NA,NA,1,NA,NA,NA,1,NA,NA,NA,1,NA,1,1,2,NA,3,3,2),
                          `KO mice`=c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,NA,NA), 
                          OMIM=c(NA,NA,NA,1,1,2,NA,NA,NA,NA,NA,2,NA,NA,NA,NA,NA,NA,NA),
                          HC=c(3,NA,1,NA,1,2,NA,2,2,NA,NA,NA,NA,NA,NA,2,NA,NA,2),
                          Missense=c(3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,3),
                          OA.phen=c("both", "a_knee", "a_knee", "both", "both", "both", "z_hip", "both", 
                                    "a_knee", "both", "both", "both", "a_knee", "both","a_knee", "a_knee", 
                                    "a_knee", "a_knee", "a_knee"))

# order.genes <- c("HOXA9","APOE","ENHO","EPRS","FTO","GLIS3","IRX3","MSRA","MYO7A","OPA3","RTN2","TMEM119","WSCD2",
#                  "DIABLO","JADE2","RARRES2","SMARCD3","TCF7L2","TMEM176A")
hc.genes.dt <- hc.genes.dt[order(OA.phen)]

long.dt <- melt(hc.genes.dt, 
                id.vars=c("genes", "OA.phen"), 
                measure.vars=c("mQTL", "DEG", "KO mice", "OMIM", "HC", "Missense"))

#reorder factors
mylevels <- hc.genes.dt$genes
long.dt$genes <- factor(long.dt$gene,levels=mylevels)

my.colors <- RColorBrewer::brewer.pal(8, "Set2")[c(1,2,6)]

ggplot(long.dt, aes(variable,genes, fill=cut(value, breaks=0:3, labels=c("OA", "T2D", "both")))) +
  geom_tile(colour="black") +
  theme_classic() +
  scale_fill_manual(values=my.colors, name="", na.translate = F) +
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_text(angle=-50,hjust=1), legend.position="top",axis.text.y=element_text(face="italic")) +
  scale_x_discrete(position = "top") 
  # scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave("/project_data/overlap_T2D_OA/paper/figures/hc_genes_figure.png", width=2.8, height=6.5, dpi=320)

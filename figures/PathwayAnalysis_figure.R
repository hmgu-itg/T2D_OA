library(data.table)
library(ggplot2)
library(RColorBrewer)

my_colors <- RColorBrewer::brewer.pal(9, "RdPu")[c(3,5,8)]

dt <- fread("/project_data/overlap_T2D_OA/pathway_analysis/pathway_table.csv")
dt[, gene_set:=ifelse(grep("hip", GeneSet), "hip", ifelse(grep("knee", GeneSet), "knee", "all")), by=seq_len(nrow(dt))]

likely.dt <- dt[grep("likely", GeneSet)]
hc.dt <- dt[grep("hc", GeneSet)] %>% .[, gene.set:=c(rep("all",2), "knee", rep("hip", 18))]
# hc.dt[, gene_set:=ifelse(grep("hc_hip", GeneSet), "hip", ifelse(grep("hc_knee", GeneSet), "knee", "all")), by=seq_len(nrow(hc.dt))]

ggplot(hc.dt[c(seq(1,4),7)], aes(x=-log(p.adjust), y=Pathway, fill=gene.set)) +
  geom_col() +
  facet_grid(rows = vars(gene.set),scales = "free_y",space = "free_y") +
  theme_bw() + 
  scale_fill_manual(values=RColorBrewer::brewer.pal(8, "Paired")[c(3,2,4)]) + 
  #scale_fill_manual(values = my_colors) +
  ylab("") +
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.position = "none")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60))

ggsave(paste0("/project_data/overlap_T2D_OA/paper/figures/pathway_analysis.png"), width=9, height=3, dpi=320)


################### ConsensusPathDb #####################
my_colors <- RColorBrewer::brewer.pal(9, "RdPu")[c(3,5,8)]

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

dt <- fread("/project_data/overlap_T2D_OA/pathway_analysis/ConsensusPathDb/high_conf_genes.csv")
dt[, gene.set:=ifelse(grepl("hip", GeneSet), "hip", ifelse(grepl("knee", GeneSet), "knee", "all")), by=seq_len(nrow(dt))]
plot.dt <- dt[c(seq(1,3),seq(44,46),seq(75,77))]
plot.dt2 <- dt[c(1,2,44,45,75,77)]
plot.dt2 <- plot.dt2[, pathway:=firstup(pathway)]

ggplot(plot.dt2[gene.set %in% c("knee","hip")], aes(x=-log(`q-value`), y=pathway, fill=gene.set)) +
  geom_col() +
  facet_grid(rows = vars(gene.set),scales = "free_y",space = "free_y") +
  theme_bw() + 
  # scale_fill_manual(values=RColorBrewer::brewer.pal(8, "Paired")[c(3,2,4)]) + 
  scale_fill_manual(values = my_colors[c(2,3)]) +
  ylab("") +
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.position = "none")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60))

ggsave(paste0("/project_data/overlap_T2D_OA/paper/figures/pathway_analysis_pink.png"), width=9, height=3, dpi=320)

library(data.table)
library(openxlsx)
library(dplyr)

DEG.OA.all <- as.data.table(read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/overlap_T2D_OA/data/DEGs_OA.xlsx", startRow = 2))
DEG.OA <- DEG.OA.all[abs(logFC)>log(1.5, base=2) & adj.P.Val<=0.05]

DEG.T2D <- as.data.table(read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/overlap_T2D_OA/data/DEGs_T2D.xlsx", startRow = 3))

genes.dt <- fread("C:/Users/ana.arruda/Documents/OA_T2D/overlap_T2D_OA/list_genes.csv")
genes.dt <- unique(genes.dt, by="id")
hc.genes <- c("FTO", "GLIS3", "WSCD2", "IRX3", "TCF7L2", "APOE", "TMEM176A", "JAZF1")

OA <- DEG.OA[DEG.OA$GENEID %in% genes.dt$id] %>% .[, .(gene=GENENAME, ID=GENEID, direction=ifelse(direction=="up", "+", "-"))]
T2D <- DEG.T2D[DEG.T2D$Ensembl.gene %in% genes.dt$id] %>% .[, .(gene=Gene.name, ID=Ensembl.gene, direction=ifelse(logFC>0, "+", "-"))]
                                                                
# DEG.table <- cbind(data.table(gene=genes.dt$name, rbindlist(combined.result)) %>% .[, OA:=ifelse(rowSums(.SD)>0,1,0), .SDcols = 2:4] %>% .[, T2D:=ifelse(rowSums(.SD)>0,1,0), .SDcols = 5:7]
DEG.table <- data.table(gene=genes.dt$name, ID=genes.dt$id) %>% .[, `:=` (OA=ifelse(ID %in% OA$ID, 1, 0), T2D=ifelse(ID %in% T2D$ID, 1, 0))]
fwrite(DEG.table, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/DEG_table.csv")

library(dplyr)
library(httr)
library(data.table)
library(openxlsx)
library(ggplot2)

# --------------------- Get GO high-confidence genes ----------------------
GO.topgenes <- read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/Boer_et_al_Suppl-Table_10.xlsx") %>% as.data.table() %>% .[X18>=3 & X19!="Gene", X19] %>% unique()
length(GO.topgenes)

# --------------------- Get DIAMANTE high-confidence genes ----------------------
DIAMANTE.topgenes <- vector()
DIAMANTE.topgenes <- Filter(Negate(is.na), lapply(content(GET("https://kp4cd.org/egldata/dataset?dataset=egls&trait=t2d"))$data, 
                                                  function(i) {DIAMANTE.topgenes <- ifelse(i$Score>=3, c(DIAMANTE.topgenes, i$gene), DIAMANTE.topgenes); return(DIAMANTE.topgenes)})) %>% as.character()
length(DIAMANTE.topgenes)

# --------------------- MAIN ----------------------
genes.dt <- fread("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/list_genes.csv")
genes.dt <- unique(genes.dt, by="id")

genes.dt <- as.data.table(read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/overlap_T2D_OA/GWAS_comorbidity/list_genes_OA.T2D.xlsx"))
setnames(genes.dt, c("name", "id"))

OA <- GO.topgenes[GO.topgenes %in% genes.dt$name]
length(OA)
# OA <- genes.dt[name %in% DIAMANTE.topgenes] %>% .[, .(gene=name, ID=id)]
T2D <- DIAMANTE.topgenes[DIAMANTE.topgenes %in% genes.dt$name]
length(T2D)

HC.table <- data.table(gene=genes.dt$name, ID=genes.dt$id) %>% .[, `:=` (OA=ifelse(gene %in% OA, 1, 0), T2D=ifelse(gene %in% T2D, 1, 0))]
fwrite(HC.table, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/HC_table.csv")


  


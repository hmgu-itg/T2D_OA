library(data.table)
library(stringr)

read.data <- function(gene.set){
  dt <- fread(paste0("/project_data/overlap_T2D_OA/pathway_analysis/ConsensusPathDb/", gene.set, "_ORA_results.tab")) %>%
    .[, .(`p-value`,  `q-value`, pathway, source, external_id, members_input_overlap, size, effective_size)]
  dt.GO <- fread(paste0("/project_data/overlap_T2D_OA/pathway_analysis/ConsensusPathDb/", gene.set, "_GO_ORA_results.tab")) %>%
    .[, source:="GO"] %>%
    .[, .(`p-value`,  `q-value`, pathway=term_name, source, external_id=term_goid, members_input_overlap=members_input_overlap_geneids, size, effective_size)]
  dt <- rbind(dt, dt.GO) %>%
    .[, `:=`(GeneSet=gene.set, N_genes=str_count(members_input_overlap, ";")+1)] %>%
    dplyr::arrange(log(`q-value`))
  setcolorder(dt, c(colnames(dt)[9], colnames(dt)[1:8])) 
  return(dt)
}

dt <- data.table()
for (gene.set in c("hc", "hc_knee", "hc_hip")) {
  dt <- rbind(dt, read.data(gene.set))
  fwrite(dt, "/project_data/overlap_T2D_OA/pathway_analysis/ConsensusPathDb/high_conf_genes.csv")
}

dt <- data.table()
for (gene.set in c("likely", "likely_knee", "likely_hip")){
  dt <- rbind(dt, read.data(gene.set))
  fwrite(dt, "/project_data/overlap_T2D_OA/pathway_analysis/ConsensusPathDb/likely_genes.csv")
}





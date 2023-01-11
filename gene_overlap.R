library(data.table)
library(httr)
library(jsonlite)
library(xml2)
library(dplyr)

overlap_ensembl <- function(chr, start, end){
  server <- "https://grch37.rest.ensembl.org"
  ext <- paste0("/overlap/region/human/", chr, ":", start, "-", end, "?feature=gene")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  
  gene.dt <- as.data.table(fromJSON(toJSON(content(r)))) %>% .[, .(name=as.character(external_name), id=as.character(id))]

  return(gene.dt)
}

genes.dt <- data.table(name=character(), id=character(), OA.phen=character())

for (trait in c("AllOA", "KneeOA", "HipOA", "KneeHipOA", "TKR", "THR", "TJR")) {
  # credset.dt <- fread(paste0("/project_data/processed_data/GWAS", trait, "/credible_set2.csv"))
  credset.dt <- fread(paste0("/project_data/processed_data/GWAS", trait, "/credible_set.csv"))
  lead.dt <- credset.dt[credset.dt[, .I[which.max(SNP.PP4)], by=region]$V1] 
  lead.dt <- lead.dt[, `:=` (start=ifelse(POS-1000000<0, 0, POS-1000000), end=POS+1000000)]
  
  tmp <- lead.dt[, overlap_ensembl(CHR, start, end), by=seq_len(nrow(lead.dt))]
  genes.dt <- rbind(genes.dt, tmp, fill=TRUE)
  genes.dt <- unique(genes.dt, by=c('name', 'id'))
  genes.dt <- genes.dt[name %in% tmp$name, OA.phen:=ifelse(is.na(OA.phen), trait, paste(c(OA.phen, trait), collapse=",")), by=id]
}

genes.dt <- genes.dt[, joint:=ifelse(grepl("KneeOA|TKR", OA.phen), ifelse(grepl("\\bHipOA\\b|THR", OA.phen), "both", "knee"), ifelse(grepl("\\bHipOA\\b|THR", OA.phen), "hip", "both"))]

fwrite(genes.dt[, .(name,id, OA.phen, joint)], "/project_data/overlap_T2D_OA/list_genes_phen.csv")

egenes <- c("FAM150B", "JADE2", "JAZF1", "HIBADH", "TMEM176A", "FASTK", "MTMR9", "MSRA", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", "AF131215.8",
           "RP11-297N6.4", "SIGMAR1", "TCF7L2", "TSKU", "WSCD2", "UNG", "DIABLO", "ARL6IP4", "RP11-197N18.2", "TCTN2", "RILPL2", "CDK2AP1", "SERPINA1", "IRX3", "APOE", "ZNF226", "SNRPD2")

# dt <- fread("/project_data/overlap_T2D_OA/list_genes_phen.csv")
# dt <- dt[, joint:=ifelse(grepl("KneeOA|TKR", OA.phen), ifelse(grepl("\\bHipOA\\b|THR", OA.phen), "both", "knee"), ifelse(grepl("\\bHipOA\\b|THR", OA.phen), "hip", "both"))]


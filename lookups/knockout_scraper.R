library(data.table)
library(httr)
library(jsonlite)
library(xml2)
library(dplyr)
library(openxlsx)

########################### Functions ###############################
.simpleCap <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

search.IMPC <- function(gene){
  server <- "https://www.ebi.ac.uk"
  ext <- paste0("/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:", gene, "&rows=500")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  
  IMPC.dt <- as.data.table(fromJSON(toJSON(content(r)))$response$docs) 
  if(length(IMPC.dt)!=0){
    # IMPC.dt <- IMPC.dt[, .(phenotype=mp_term_name, effect=effect_size, pval=p_value, class=parameter_name, source="IMPC")]
    IMPC.dt <- IMPC.dt[, .(phenotype=as.character(mp_term_name), source="IMPC")]
    return(IMPC.dt)
  }
  else{
    return(IMPC.dt)
  }
}

mgi.ret <- function(gene, dt){
  if (is.na(dt[Input==gene][1]$Term)) {
    return(data.table())
  }
  else{
    return(dt[Input==gene, .(phenotype=Term, source="MGI")])
  }
}

search.RGD <- function(gene, specie=2, src="MP"){
  # species: 1=human, 2=mouse, 3=rat
  # source: mouse phenotype=MP, human pheotype=HP, diseases=DOID
  
  gene <- tolower(gene)
  
  server <- "https://rest.rgd.mcw.edu"
  r <- GET(paste(server, "/rgdws/genes/", gene, "/", specie, sep = ""), content_type("application/json"))
  stop_for_status(r)
  if (length(fromJSON(toJSON(content(r))))>1){
    rgdID <- fromJSON(toJSON(content(r)))$rgdId
    
    ext <- paste0("/rgdws/annotations/rgdId/", rgdID, "/", src)
    r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
    stop_for_status(r)
    
    if(length(fromJSON(toJSON(content(r))))!=0){
      phen <- fromJSON(toJSON(content(r)))$term
      RGD.dt <- data.table(phenotype=as.character(phen), source="RGD")
      return(RGD.dt)
    }
    else return(data.table())
  }
  else return(data.table())
}

knockout.lookup <- function(data, sources, terms){
  result <- data.table(t(rep(0, length(sources))))
  colnames(result) <- c(paste(sources, terms, sep="_"))
  if(nrow(data)!=0){
    pattern <- paste(get(terms), collapse="|")
    tmp <- data[, grepl(pattern, phenotype, ignore.case=TRUE), by=.(source, seq_len(nrow(data)))]
    for (src in sources) {
      if(sum(tmp[source==src]$V1)>0) result <- result[, (paste(src, terms, sep="_")):=1]
    }
  }
  return(result)
}

########################### Main ###############################
genes.dt <- fread("/project_data/processed_data/list_genes.csv")
genes.dt <- unique(genes.dt, by="id")
genes.dt <- genes.dt[name %in% c("PITX1", "TCF7L2")]

#------------------- add IMPC ----------------------#
impc.lst <- lapply(.simpleCap(genes.dt$name), search.IMPC)
names(impc.lst) <- genes.dt$name

#------------------- add MGI ----------------------#
# http://www.informatics.jax.org/batch/summary
# mgi <- as.data.table(read.xlsx("/project_data/processed_data/MGIBatchReport_20210907_110229.xlsx")) %>% .[, .(Input, Term)]
mgi <- fread("~/T2D_OA/MGIBatchReport_20220215_053516.csv") %>% .[, .(Input, Term)]
mgi.lst <- lapply(genes.dt$name, mgi.ret, dt=mgi)
# names(mgi.lst) <- unique(mgi$Input)[1:100]
names(mgi.lst) <- genes.dt$name

#------------------- add RGD ----------------------#
rgd.lst <- lapply(genes.dt$name, search.RGD)
names(rgd.lst) <- genes.dt$name

#------------------ Combine all ------------------ #
knockout.lst <- Map(rbind, Map(rbind, impc.lst, mgi.lst), rgd.lst)
names(knockout.lst) <- genes.dt$name

#------------------- check for OA and T2D terms ----------------------#
OA <- c("bone", "muscle", "skeleton", "osteo", "arthritis", "muscular", "joint", "body size", "growth", "skeletal", "stature", "height")
        # "scoliosis", "Joubert Syndrome", "Cutislaxa", "Boissel", "atrophy", "distal", "amyotrophic lateral sclerosis", "rhabdomyolysis", "polymyositis", "brachydactyly"
T2D <- c("insulin", "glycemia", "glucose", "diabetes", "pancreas", "beta-cell", "beta cell", "beta cells",
         "glucosuria", "body weight", "obesity", "BMI", "body mass", "body fat", "hyperglycemia", "pancreatic")

OA.result <- lapply(knockout.lst, knockout.lookup, sources=c("IMPC", "MGI", "RGD"), terms="OA")
T2D.result <- lapply(knockout.lst, knockout.lookup, sources=c("IMPC", "MGI", "RGD"), terms="T2D")
combined.result <- Map(cbind, OA.result, T2D.result)

#---------------- output result -------------------#
knockout.table <- cbind(data.table(gene=names(combined.result)), rbindlist(combined.result)) %>% .[, OA:=ifelse(rowSums(.SD)>0,1,0), .SDcols = 2:4] %>% .[, T2D:=ifelse(rowSums(.SD)>0,1,0), .SDcols = 5:7]
# fwrite(knockout.table, "/project_data/processed_data/knockout_table.csv")
fwrite(knockout.table, "~/T2D_OA/knockout_table.csv")

knockout.table.old <- fread("~/T2D_OA/knockout_table_old.csv")
setdiff(knockout.table.old, knockout.table)

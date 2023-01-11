library(data.table)
library(httr)
library(jsonlite)
library(xml2)
library(dplyr)
library(openxlsx)

########################### Functions ###############################
search.RGD.OMIM <- function(gene){
  gene <- tolower(gene)
  
  server <- "https://rest.rgd.mcw.edu"
  r <- GET(paste(server, "/rgdws/genes/", gene, "/1", sep = ""), content_type("application/json"))
  stop_for_status(r)
  if (length(fromJSON(toJSON(content(r))))>1){
    rgdID <- fromJSON(toJSON(content(r)))$rgdId
    
    ext <- paste0("/rgdws/annotations/rgdId/", rgdID, "/DOID")
    r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
    stop_for_status(r)
    
    
    if(length(fromJSON(toJSON(content(r))))!=0){
      OMIM.dt <- as.data.table(fromJSON(toJSON(content(r)))[fromJSON(toJSON(content(r)))$dataSrc=="OMIM", c("term", "qualifier")]) 
      phen <- OMIM.dt[as.character(qualifier)!="susceptibility", term]
      return(data.table(phenotype=phen, source="OMIM"))
    }
    else return(data.table())
  }
  else return(data.table())
}

search.OMIM <- function(gene){
  gene <- tolower(gene)
  
  server <- "https://api.omim.org"
  r <- GET(paste(server, "/api/geneMap/search?search=", gene, "&include=geneMap&apiKey=cBg3cWK0TxibU0wvtkwtMA&format=json&start=0&limit=100", sep = ""), content_type("application/json"))
  stop_for_status(r)
  if (length(fromJSON(toJSON(content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList)!=0){
    # OMIM.dt <- as.data.table(fromJSON(toJSON(content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList[[1]])
    # all.phen <- OMIM.dt$phenotypeMap.phenotype
    all.phen <- fromJSON(toJSON(content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList[[1]]$phenotypeMap$phenotype
    phen <- all.phen[sapply(all.phen, function(i){!grepl("\\{|\\[|\\?",i)})]
    if (length(phen)==0) return(data.table()) 
    else return(data.table(phenotype=phen, source="OMIM"))
  }
  else return(data.table())
}

omim.lookup <- function(data, source, terms){
  result <- data.table(t(rep(0, 1)))
  colnames(result) <- c(paste(source, terms, sep="_"))
  if(nrow(data)!=0){
    pattern <- paste(get(terms), collapse="|")
    tmp <- data[, grepl(pattern, phenotype, ignore.case=TRUE), by=seq_len(nrow(data))]
    if(tmp$V1) result <- result[, (paste(source, terms, sep="_")):=1]
  }
  return(result)
}

########################### Main ###############################
genes.dt <- fread("/project_data/processed_data/list_genes.csv")
genes.dt <- unique(genes.dt, by="id")
# genes.dt <- genes.dt[1:100]

# OA.genes.dt <- as.data.table(read.xlsx("list_genes_OA.T2D.xlsx"))
# setnames(OA.genes.dt, c("name", "id"))
# genes.dt <- OA.genes.dt
  
omim.lst <- sapply(genes.dt$name, search.OMIM)
names(omim.lst) <- genes.dt$name

#------------------- check for OA and T2D terms ----------------------#
OA <- c("bone", "muscle", "skeleton", "osteo", "arthritis", "muscular", "joint", "body size", "growth", "skeletal", "stature", "Hand-foot-uterus", "synostosis", 
        "MARTSOLF", "WARBURG", "LEUKODYSTROPHY", "height", "SQUALENE", "FINCA", "POLYDACTYLY")
# "scoliosis", "Joubert Syndrome", "Cutislaxa", "Boissel", "atrophy", "distal", "amyotrophic lateral sclerosis", "rhabdomyolysis", "polymyositis", "brachydactyly"
T2D <- c("insulin", "glycemia", "glucose", "diabetes", "pancreas", "beta-cell", "beta cell", "beta cells",
         "glucosuria", "body weight", "obesity", "BMI", "body mass", "body fat", "hyperglycemia", "pancreatic",
         "MARTSOLF", "ACIDURIA", "Aicardi-Goutières", "FINCA")

OA.result <- lapply(omim.lst, omim.lookup, source="OMIM", terms="OA")
T2D.result <- lapply(omim.lst, omim.lookup, source="OMIM", terms="T2D")
combined.result <- Map(cbind, OA.result, T2D.result)

#---------------- output result -------------------#
omim.table <- cbind(data.table(gene=names(combined.result)), rbindlist(combined.result))
# fwrite(omim.table, "/project_data/processed_data/omim_table.csv")
fwrite(omim.table, "~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/omim_table.csv")

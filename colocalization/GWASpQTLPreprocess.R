setwd("~/gwas_eqtl_colocalization")

library(data.table)
library(openxlsx)
library(stringr)
library(hash)

#########################################
#---------- Helper functions ------------
#########################################
sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
}) 

write_output <- function(var, dict, filename){
  # var can be one variable or a list
  fwrite(dict[[paste(var, collapse="_")]], filename) 
}

add_key_value <- function(hash, key, new.value){
  if(has.key(key, hash)){
    hash[key] <- append(hash[[key]], new.value)
  }
  else{
    .set(hash, key, new.value)
  }
}


###################################################
#---------- Problem specific functions ------------
###################################################
read_pQTL_pGenes <- function(tissue, source, dict_pQTL, dict_pGenes, delete.indels=FALSE){
  ##################################
  # Input can be a data table with:
  # col_1 == tissue
  # col_2 == source
  # eQTL = dictionary
  ##################################
  
  if (source == "FunGen") {
    if (!file.exists(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))) {
      print("FunGen preprocessed file does not exist")
      tmp <- fread(paste("/project_data/data_original/FunGen_pQTL_", tissue, ".FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz", 
                         sep=""), verbose=FALSE)
      
      dict_pQTL[[tissue]] <- tmp[, .(CHR=as.integer(sub(":.*", "", genotype_id)), 
                                     POS=as.integer(sub(".*:", "", genotype_id)), rsID=rep(".", nrow(tmp)),
                                     geneID=sub(".*_", "", phenotype_id), pval, beta=slope, se=slope_se, 
                                     MAF=as.numeric(str_match(INFO, "MAF=\\s*(.*?)\\s*;")[,2]),
                                     N=PERM_num_var, EA=ALT, NEA=REF)]
      dict_pGenes[[tissue]] <- tmp[pGene==1, .(CHR=as.integer(sub(":.*", "", genotype_id)), 
                                               POS=as.integer(sub(".*:", "", genotype_id)),
                                               geneID=sub(".*_", "", phenotype_id))]
      
      write_output(tissue, dict_pQTL, paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))
      write_output(tissue, dict_pGenes, paste0("/project_data/processed_data/pGenes_", tissue, ".csv"))
    }
    
    else{
      print("FunGen preprocessed file already exists")
      dict_pQTL[[tissue]] <- fread(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"), verbose=FALSE)
      dict_pGenes[[tissue]] <- fread(paste0("/project_data/processed_data/pGenes_", tissue, ".csv"), verbose=FALSE)
    }
  }

  else if (source == "InsPIRE"){
    if (!file.exists(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))) {
      print("InsPIRE preprocessed file does not exist")
      tmp <- fread(paste("/project_data/data_original/InsPIRE_", tissue, "_Gene_eQTLs_Nominal_Pvalues.txt.gz", 
                         sep=""), verbose=FALSE)
      tmp_pGenes <- fread(paste("/project_data/data_original/", tissue, "_independent_gene_eQTLs.txt", sep=""), verbose=FALSE)
      
      col_names <- colnames(tmp)[-1]
      tmp <- tmp[, !"FreqALT"]
      setnames(tmp, col_names)
      
      dict_pQTL[[tissue]] <- tmp[, .(CHR=SNPchr, POS=SNPposition, rsID=rsID, geneID=sub("\\..*", "",GeneID), 
                                     pval=Pvalue, beta=Slope, se=SE, MAF=ifelse(FreqREF<FreqALT, FreqREF, FreqALT), 
                                     N=NumSNPtest, EA=ALT, NEA=REF)]
      dict_pGenes[[tissue]] <- unique(tmp_pGenes, by="GeneID")[, .(CHR=ChrPheno, 
                                                                   POS=as.integer(StartPheno),
                                                                   geneID=sub("\\..*", "",GeneID))]
      # delete "MERGED"
      dict_pQTL[[tissue]] <- dict_pQTL[[tissue]][!rsID %like% "MERGED"] %>% .[!rsID %like% "rs", rsID:="."]
      
      if (delete.indels){
        dict_pQTL[[tissue]] <- dict_pQTL[[tissue]][!(EA %in% c("D", "I"))]
        dict_pQTL[[tissue]] <- dict_pQTL[[tissue]][nchar(EA)==1 & nchar(NEA)==1]
      }
      
      write_output(tissue, dict_pQTL, paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))
      write_output(tissue, dict_pGenes, paste0("/project_data/processed_data/pGenes_", tissue, ".csv"))
    }
    else {
      print("InsPIRE preprocessed file already exists")
      dict_pQTL[[tissue]] <- fread(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"), verbose=FALSE)
      dict_pGenes[[tissue]] <- fread(paste0("/project_data/processed_data/pGenes_", tissue, ".csv"), verbose=FALSE)
    }
  }
}

read_GWAS <- function(trait, source, dict, indep.signals=FALSE, delete.indels=FALSE){
  # read independent signals file
  if (indep.signals==TRUE && source == "GO"){
    dict[[trait]] <- as.data.table(
      read.xlsx("/project_data/data_original/GO_independent_signals.xlsx", startRow = 2))[PHENO == trait,
                                                                                          .(SNPID=SNV, CHR=as.integer(sub(":.*", "", `CHR:POS`)), 
                                                                                            POS = as.integer(sub(".*:", "", `CHR:POS`)))]
  }
  
  else if (indep.signals==TRUE && source == "DIAGRAM"){
    dict[[trait]] <- as.data.table(
      read.xlsx("/project_data/data_original/DIAGRAM_independent_signals.xlsx", 
                startRow = 2))[, .(SNPID=Index.variant, CHR=as.integer(Chromosome), 
                                   POS = as.integer(`Position.(Build.37.bp)`))]
  }
  
  # read GWAS summary statistics
  else if (indep.signals==FALSE){
    if (source == "GO") {
      if (!file.exists(paste0("/project_data/processed_data/GWAS_", trait, ".csv"))) {
        print("GO preprocessed file does not exist")
        dict[[trait]] <- fread(paste("/project_data/data_original/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz", sep=""), 
                               verbose=FALSE)[, .(CHR=CHR, POS=POS, rsID=SNP, pval=P, beta=BETA, se=SE, MAF=MAF, 
                                                  N=N, Ncases=NCASES, EA=toupper(EA), NEA=toupper(NEA))]
        # change rsID=<NA> to rsID=.
        dict[[trait]][is.na(rsID), rsID:="."]
        
        if (delete.indels){
          dict[[trait]] <- dict[[trait]][!(EA %in% c("D", "I"))]
        }
        write_output(trait, dict, paste0("/project_data/processed_data/GWAS_", trait, ".csv"))
      }
      
      else {
        print("GO preprocessed file already exists")
        dict[[trait]] <- fread(paste0("/project_data/processed_data/GWAS_", trait, ".csv"), verbose=FALSE)
        print(colnames(dict[[trait]]))
        # setnames(dict[[trait]], c("REF", "ALT"), c("EA", "NEA"))
        # write_output(trait, dict, paste0("/project_data/processed_data/GWAS_", trait, ".csv"))
      }
    }
    else if (source == "DIAGRAM"){
      if (!file.exists(paste0("/project_data/processed_data/GWAS_", trait, ".csv"))) {
        print("DIAGRAM preprocessed file does not exist")
        tmp <- fread(paste("/project_data/data_original/Mahajan.NatGenet2018b.", trait, ".European.txt.gz", sep=""), 
                     verbose=FALSE)
        dict[[trait]] <- tmp[, .(CHR=Chr, POS=Pos, rsID=rep(".", nrow(tmp)), pval=as.numeric(Pvalue), beta=as.numeric(Beta), 
                                 se=as.numeric(SE), MAF=fifelse(EAF < (1-EAF), EAF, (1-EAF)), N=Neff, 
                                 Ncases=as.integer(74124), EA=EA, NEA=NEA)]
        # if (delete.indels){
        #   dict[[trait]] <- dict[[trait]][!(EA %in% c("D", "I"))]
        #   dict[[trait]] <- dict[[trait]][nchar(EA)==1 & nchar(NEA)==1]
        # }
        write_output(trait, dict, paste0("/project_data/processed_data/GWAS_", trait, ".csv"))
      }
      
      else{
        print("DIAGRAM preprocessed file already exists")
        dict[[trait]] <- fread(paste0("/project_data/processed_data/GWAS_", trait, ".csv"), 
                               verbose=FALSE)[, `:=` (MAF=as.numeric(MAF), Ncases=as.integer(Ncases))]
        print(colnames(dict[[trait]]))
        # setnames(dict[[trait]], c("REF", "ALT"), c("EA", "NEA"))
        # write_output(trait, dict, paste0("/project_data/processed_data/GWAS_", trait, ".csv"))
      }
    }
  }
}

read_TSS <- function(tissue, source, dict, dict_pGenes){
  if (source=="FunGen"){
    dict[[tissue]] <- fread(paste("/project_data/data_original/TSS_FunGen.txt", sep=""), 
                            verbose=FALSE)[, .(CHR=`#Chr`, START=start, ID=sub(".*_", "", ID))]
  }
  else if (source=="InsPIRE"){
    dict[[tissue]] <- dict_pGenes[[tissue]][, .(CHR=CHR, START=POS, ID=geneID)]
  }
}

select_pGenes <- function(tissue, dict_signals, dict_pGenes, wd=100000) {
  pGenes <- data.table()
  for (i in 1:nrow(dict_signals)){
    pGenes <- rbind(pGenes, dict_pGenes[[tissue]][CHR %in% dict_signals$CHR[i] &
                                                  POS %between% .(dict_signals$POS[i]-wd, dict_signals$POS[i]+wd)])
  }
  dict_pGenes[[paste(tissue, "signals", sep="_")]] <- pGenes[, unique(geneID)]
}

select_pQTL <- function(tissue, dict_pGenes, dict_eQTL) {
  key <- paste(tissue, "signals", sep="_")
  dict_pQTL[[key]] <- rbind(dict_pQTL[[key]], dict_pQTL[[tissue]][geneID %in% dict_pGenes[[key]]])
}

select_pQTL2 <- function(tissue, dict_pQTL, credset.id) {
  dict_pQTL[[paste(tissue, "coloc", sep="_")]] <- dict_pQTL[[tissue]][ID %in% credset.id]
}

select_GWAS <- function(trait, dict_GWAS, TSS, wd=1000000) {
  GWAS <- data.table()
  for (i in 1:nrow(TSS)) {
    GWAS <- rbind(GWAS, dict_GWAS[[trait]][CHR %in% TSS$CHR[i] & POS %between% .(TSS$START[i]-wd, TSS$START[i]+wd)])
    GWAS <- unique(GWAS)
  }
  dict_GWAS[[paste(trait, "signals", sep="_")]] <- GWAS
}

select_GWAS2 <- function(trait, dict_GWAS, credset.id) {
  dict_GWAS[[paste(trait, "coloc", sep="_")]] <- dict_GWAS[[trait]][ID %in% credset.id]
}

add_ID <- function(dict, trait, listID){
  dict[[trait]] <- dict[[trait]][, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
  dict[[trait]][, add_key_value(listID, ID, trait), by=ID]
}

flip_alleles <- function(dict_GWAS, dict_pQTL, IDlist, logFile) {
  for (id in keys(IDlist)){
    if (length(IDlist[[id]])>1){
      ref <- IDlist[[id]][1]
      # dict_ref <- ifelse(has.key(ref, dict_pQTL), dict_pQTL[[ref]], dict_GWAS[[ref]])
      if(has.key(ref, dict_pQTL)) refEA <- dict_pQTL[[ref]][ID==id, EA][1]
      else refEA <- dict_GWAS[[ref]][ID==id, EA]
      lapply(IDlist[[id]][-1], function(i){
        # dict <- ifelse(i %in% TISSUE$tissue, dict_eQTL[[i]], dict_GWAS[[i]])
        if(has.key(i, dict_pQTL)) altEA <- dict_pQTL[[i]][ID==id, EA][1]
        else altEA <- dict_GWAS[[i]][ID==id, EA]
        if (altEA != refEA) {
          # print(paste0("FLIPPING beta of ID=", id, " and trait=", i, " with REF=", altEA, " using ref=", ref, " with REF=", refEA))
          cat(paste0("FLIPPING beta of ID=", id, " and trait=", i, " with REF=", altEA, " using ref=", ref, " with REF=", refEA), 
              file=logFile, append=TRUE, sep="\n")
          if(has.key(i, dict_pQTL)) dict_pQTL[[i]][ID==id, beta:=-beta]
          else dict_GWAS[[i]][ID==id, beta:=-beta]
        }
      })
    }
  }
}


#########################################
#--------------- Load data --------------
#########################################
# TISSUE <- data.table(tissue=c("HighGradeCartilage", "LowGradeCartilage",
#                               "Synovium", "Whole_Blood", "PancreaticIslets"),
#                      source=c(rep(c("FunGen"),3), "GTEx", "InsPIRE"))
TISSUE <- data.table(tissue=c("HighGradeCartilage", "LowGradeCartilage",
                              "PancreaticIslets"),
                     source=c(rep(c("FunGen"),2), "InsPIRE"))
TRAIT <- data.table(trait=c("T2D", "KneeOA"),
                    source=c("DIAGRAM", rep(c("GO"),1)))

#----------------------------------------
# Create dictionaries
#----------------------------------------
pQTL <- hash()
pGenes <- hash()
GWAS <- hash()
TSS <- hash()

#----------------------------------------
# Load eQTL and eGenes data
#----------------------------------------
print("Loading eQTL data")
TISSUE[, read_pQTL_pGenes(tissue=tissue, dict_pQTL=pQTL, dict_pGenes=pGenes, source=source, delete.indels=TRUE), by=tissue]

#----------------------------------------
# Load TSS for eGenes
#----------------------------------------
print("Loading TSS for eGenes data")
TISSUE[, read_TSS(tissue = tissue, source=source, dict = TSS, dict_eGenes = eGenes), by=tissue]

#----------------------------------------
# Load GWAS data
#----------------------------------------
print("Loading GWAS data")
TRAIT[, read_GWAS(trait=trait, dict=GWAS, source=source, indep.signals=FALSE), by = seq_len(nrow(TRAIT))]

#----------------------------------------
# Load colocalized GWAS regions
#----------------------------------------
print("Loading colocalized GWAS regions")
# regions <- fread("/project_data/processed_data/GWAS_colocalized_regions.csv")
colocalized_SNP <- fread("/project_data/processed_data/GWAS_colocalized_SNP.csv")

#----------------------------------------
# Get credible set
#----------------------------------------
credset.dt <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/credible_set.csv"))

################################################################################
#--------------- Select rows of GWASes in credible set --------------
################################################################################
TRAIT[, select_GWAS2(trait=trait)]


################################################################################
#--------------- Select eQTLs in colocalized GWAS regions --------------
################################################################################
#---------------------------------------------------------------------------
# Calculate overlap of eGenes within 100kb from GWAS colocalized SNPs
#---------------------------------------------------------------------------
print("Calculating overlap of eGenes colocalized GWAS regions")
TISSUE[, select_eGenes(tissue=tissue, signals=colocalized_SNP, dict_eGenes=eGenes), by=tissue]

#-----------------------------------------------------------
# Take all variants from the selected eGenes
#-----------------------------------------------------------
print("Selecting all variants from the eGenes in GWAS relevant regions")
TISSUE[, select_eQTL(tissue, dict_eGenes=eGenes, dict_eQTL=eQTL), by=tissue]


#############################################################
#--------------- Select GWAS signals --------------
#############################################################
#----------------------------------------
# Merge TSS from all eGenes
#----------------------------------------
print("Merging TSS from all relevant eGenes")
all_TSS <- TSS[[TISSUE$tissue[1]]][ID %in% eGenes[[paste(TISSUE$tissue[1], "signals", sep="_")]]]

for (tissue in TISSUE$tissue[-1]){
  all_TSS <- merge(all_TSS, TSS[[tissue]][ID %in% eGenes[[paste(tissue, "signals", sep="_")]]], all=T)
}

all_TSS <- unique(all_TSS, by="ID")

#---------------------------------------------------------------------------------------
# Extract GWAS variants within 1Mb from TSS of relevant eGenes from all analyzed tissues
#---------------------------------------------------------------------------------------
print("Extracting GWAS variants around TSS from relevant eGenes")
TRAIT[, select_GWAS(trait, dict_GWAS=GWAS, TSS=all_TSS), by=trait]


############################################
#--------------- Flip alleles --------------
############################################
#------------------------------------------------------------
# Add ID column and get list of IDs with corresponding traits
#------------------------------------------------------------
print("Adding ID......")

IDlist <- hash()

print("..... For GWAS datasets")
for (trait in TRAIT$trait){
  add_ID(GWAS, paste(trait, "signals", sep="_"), IDlist)
}

print("..... For eQTL datasets")
for (tissue in TISSUE$tissue){
  add_ID(eQTL, paste(tissue, "signals", sep="_"), IDlist)
}

#----------------------------------------
# Flip alleles per ID and log flipped
#----------------------------------------
print("Flipping alleles when needed")
logFile <- "/project_data/processed_data/flipped.txt"
flip_alleles(dict_GWAS=GWAS, dict_eQTL=eQTL, IDlist=IDlist, logFile=logFile)


########################################################################
#--------------- Output files for colocalization analysis --------------
########################################################################
print("Outputing GWAS files for colocalization")
for (trait in TRAIT$trait) {
  fwrite(GWAS[[paste(trait, "signals", sep="_")]], paste0("/project_data/processed_data/GWAS_", trait, "_coloc.csv"))
}
print("Outputing eQTL files for colocalization")
for (tissue in TISSUE$tissue) {
  fwrite(eQTL[[paste(tissue, "signals", sep="_")]], paste0("/project_data/processed_data/eQTL_", tissue, "_coloc.csv"))
}







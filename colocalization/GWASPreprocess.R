library(data.table)
library(openxlsx)
library(stringr)
library(hash)
library(ieugwasr)

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

annotate_rsids <- function(trait, trait_dt, rsid){
  # read data
  # rsid <- fread(pste0("/project_data/processed_data/Variants2rsID/", trait,".id.rsid.bed"), verbose=FALSE)
  # trait_dt <- fread(paste0("/project_data/processed_data/GWASKneeOA/GWAS_", trait,"_precoloc_regions"), verbose=FALSE)
  
  # rename columns
  colnames(rsid) <- c("CHR:POS", "CHR", "POS", "rsID", "A1", "A2", "ID")
  
  # remove indels
  new_rsid <- rsid[nchar(A1)==1]
  
  # separate A2 column when there is a komma
  new_rsid <- separate(new_rsid, "A2", c("A2a", "A2b"), sep=",")
  new_rsid <- melt(new_rsid, measure.vars=c("A2a","A2b"), variable.name="variable", value.name="A2")
  
  # add IDs to trait_rsid
  new_rsid[, ID:=paste(`CHR:POS`, sort_alleles(A1, A2), sep="_")]
  
  # merge rsid into trait data.table
  new_trait_dt <- merge(trait_dt[, .SD, .SDcols = !"rsID"], new_rsid[, .(ID, POS=as.integer(POS), CHR=as.integer(CHR), rsID)], by=c("ID", "CHR", "POS"), all.x=TRUE)

  # save anotated dataset
  fwrite(new_trait_dt[, `CHR:POS`:=NULL], paste0("/storage/hmgu/projects/OA_T2D/processed_data/GWASKneeOA/GWAS_", trait,"_rsid_precoloc_regions.csv"), verbose=FALSE)
}

get.indep.signals <- function(data, chr, pval.col, id.col, snp.col, pval=5e-08, r2=0.1, kb=1000){
  clumped.dt <- ld_clump(dplyr::tibble(rsid=data[,get(id.col)], pval=data[,get(pval.col)]),
                         plink_bin=genetics.binaRies::get_plink_binary(),
                         bfile=paste0("/project_data/StatIndependence/ukb_bgen2plink/chr", chr),
                         clump_p=pval,
                         clump_r2=r2,
                         clump_kb=kb)
  return(clumped.dt$rsid)
}

###################################################
#---------- Problem specific functions ------------
###################################################
read_GWAS <- function(trait, source, dict, indep.signals=FALSE, delete.indels=FALSE){
  # read independent signals file
  if (indep.signals==TRUE && source == "GO"){
    dict[[trait]] <- as.data.table(
      read.xlsx("/storage/hmgu/projects/OA_T2D/data_original/GO_independent_signals.xlsx", startRow = 2))[PHENO == trait,
                                                                                          .(SNPID=SNV, CHR=as.integer(sub(":.*", "", `CHR:POS`)), 
                                                                                            POS=as.integer(sub(".*:", "", `CHR:POS`)),
                                                                                            ID=paste(paste(as.integer(sub(":.*", "", `CHR:POS`)), as.integer(sub(".*:", "", `CHR:POS`)), sep=":"), sort_alleles(EA, NEA), sep="_"))]
  }
  
  else if (indep.signals==TRUE && source == "DIAGRAM"){
    dict[[trait]] <- as.data.table(read.xlsx("/storage/hmgu/projects/OA_T2D/data_original/DIAGRAM_independent_signals.xlsx", startRow = 2)) %>%
      .[Model=="BMI unadjusted"] %>%
      # .[`Primary/Secondary`=="Primary" & Model=="BMI unadjusted"] %>%
      .[, .(SNPID=Index.variant, CHR=as.integer(Chromosome), POS=as.integer(`Position.(Build.37.bp)`), ID=paste(paste(Chromosome, `Position.(Build.37.bp)`, sep=":"), sort_alleles(Risk.allele, Other.allele), sep="_"))]
  }
  
  # read GWAS summary statistics
  else if (indep.signals==FALSE){
    if (source == "GO") {
      if (!file.exists(paste0("/project_data/processed_data/GWAS_", trait, ".csv"))) {
        print("GO preprocessed file does not exist")
        dict[[trait]] <- fread(paste("/storage/hmgu/projects/OA_T2D/data_original/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz", sep=""), 
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
        tmp <- fread(paste0("/project_data/processed_data/GWAS_", trait, ".csv"), verbose=FALSE) %>% .[, ID:="."]
        dict[[trait]] <- tmp
        # setnames(dict[[trait]], c("REF", "ALT"), c("EA", "NEA"))
        # write_output(trait, dict, paste0("/project_data/processed_data/GWAS_", trait, ".csv"))
      }
    }
    else if (source == "DIAGRAM"){
      if (!file.exists(paste0("/project_data/processed_data/GWAS_", trait, ".csv"))) {
        print("DIAGRAM preprocessed file does not exist")
        tmp <- fread(paste("/storage/hmgu/projects/OA_T2D/data_original/Mahajan.NatGenet2018b.T2D.European.txt.gz", sep=""), 
                     verbose=FALSE)
        dict[[trait]] <- tmp[, .(CHR=Chr, POS=Pos, rsID=rep(".", nrow(tmp)), pval=as.numeric(Pvalue), beta=as.numeric(Beta), 
                                 se=as.numeric(SE), MAF=fifelse(EAF < (1-EAF), EAF, (1-EAF)), N=Neff, 
                                 Ncases=as.integer(74124), EA=EA, NEA=NEA, ID=rep(".", nrow(tmp)))]
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
        # setnames(dict[[trait]], c("REF", "ALT"), c("EA", "NEA"))
        # write_output(trait, dict, paste0("/project_data/processed_data/GWAS_", trait, ".csv"))
      }
    }
  }
}

select_GWAS <- function(trait, dict_GWAS, signals, wd=1000000) {
  GWAS <- data.table()
  for (i in 1:nrow(signals)) {
    GWAS <- rbind(GWAS, dict_GWAS[[trait]][CHR %in% signals$CHR[i] & POS %between% .(signals$POS[i]-wd, signals$POS[i]+wd)])
    GWAS <- unique(GWAS)
  }
  dict_GWAS[[paste(trait, "signals", sep="_")]] <- GWAS
}

add_ID <- function(dict, trait, listID){
  dict[[trait]] <- dict[[trait]][, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
  dict[[trait]][, add_key_value(listID, ID, trait), by=ID]
}

flip_alleles <- function(dict_GWAS, IDlist) {
  for (id in keys(IDlist)){
    if (length(IDlist[[id]])>1){
      ref <- IDlist[[id]][1]
      # dict_ref <- ifelse(has.key(ref, dict_eQTL), dict_eQTL[[ref]], dict_GWAS[[ref]])
      refEA <- dict_GWAS[[ref]][ID==id, EA]
      lapply(IDlist[[id]][-1], function(i){
        # dict <- ifelse(i %in% TISSUE$tissue, dict_eQTL[[i]], dict_GWAS[[i]])
        altEA <- dict_GWAS[[i]][ID==id, EA]
        if (altEA != refEA) {
          dict_GWAS[[i]][ID==id, beta:=-beta]
        }
      })
    }
  }
}

#########################################
#--------------- Load data --------------
#########################################
for (t in c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR")){
  TRAIT <- data.table(trait=c("T2D", t),
                      source=c("DIAGRAM", rep(c("GO"),1)))
  
  #----------------------------------------
  # Create dictionaries
  #----------------------------------------
  GWAS <- hash()
  GWAS_signals <- hash()
  
  #----------------------------------------
  # Load independent GWAS signals
  #----------------------------------------
  print("Loading independent GWAS signals")
  TRAIT[, read_GWAS(trait=trait, dict=GWAS_signals, source=source, indep.signals=TRUE), by=trait]
  
  #----------------------------------------
  # Load GWAS data
  #----------------------------------------
  print("Loading GWAS data")
  TRAIT[, read_GWAS(trait=trait, dict=GWAS, source=source, indep.signals=FALSE, delete.indels=TRUE), by=trait]
  
  ################################################################################
  #--------------- Define GWAS independent signals regions --------------
  ################################################################################
  #----------------------------------------
  # Merge GWAS signals from all traits
  #----------------------------------------
  print("Merging independent GWAS signals from all traits")
  indep_signals <- GWAS_signals[[TRAIT$trait[1]]]
  # for (i in 2:nrow(TRAIT)){
  for (trait in TRAIT$trait[-1]){
    indep_signals <- na.omit(unique(merge(indep_signals, GWAS_signals[[trait]], all=T)))
  }
  
  #--------------------------------------------
  # Define regions of overlapping  GWAS signals 
  #--------------------------------------------
  wd = 1e+6
  regions <- na.omit(indep_signals[, ':='(start=POS-wd, end=POS+wd)])
  GWAS_overlap_regions <- regions[, as.data.table(IRanges::reduce(IRanges::IRanges(start, end), min.gapwidth = 0L)), CHR]
  
  ################################################################################
  #--------------- Select GWAS that overlap GWAS regions --------------
  ################################################################################
  #-----------------------------------------------------------------
  # Extract GWAS variants within 1Mb from merged independent signals
  #-----------------------------------------------------------------
  print("Selecting GWAS signals around merged independent signals")
  TRAIT[, select_GWAS(trait, dict_GWAS=GWAS, signals=indep_signals), by=trait]
  
  ############################################
  #--------------- Flip alleles --------------
  ############################################
  #------------------------------------------------------------
  # Add ID column and get list of IDs with corresponding traits
  #------------------------------------------------------------
  print("Adding ID for GWAS datasets")
  
  IDlist <- hash()
  add_ID(GWAS, paste(TRAIT$trait[2], "signals", sep="_"), IDlist)
  # for (trait in TRAIT$trait){
  #   add_ID(GWAS, paste(trait, "signals", sep="_"), IDlist)
  # }
  
  #----------------------------------------
  # Flip alleles per ID and log flipped
  #----------------------------------------
  print("Flipping alleles when needed")
  flip_alleles(dict_GWAS=GWAS, IDlist=IDlist)
  
  ########################################################################
  #--------------- Output files for colocalization analysis --------------
  ########################################################################
  print("Outputing GWAS files for colocalization")
  for (trait in TRAIT$trait) {
    fwrite(GWAS[[paste(trait, "signals", sep="_")]], paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_", trait, "_precoloc_regions.csv"))
  }
  fwrite(GWAS_overlap_regions, paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_overlap_regions.csv"))
  fwrite(regions, paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_regions.csv"))
}
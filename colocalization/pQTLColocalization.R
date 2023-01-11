setwd("~/gwas_eqtl_colocalization/T2D_OA")

library(data.table)
library(stringr)
library(coloc)
library(moloc)
library(hyprcoloc)
library(hash)
library(dplyr)
library(tidyr)

sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
}) 

write_output <- function(var, dict, filename){
  # var can be one variable or a list
  fwrite(dict[[paste(var, collapse="_")]], filename) 
}

read_pQTL <- function(tissue, source, dict_pQTL, delete.indels=FALSE){
  if (source == "FunGen") {
    if (!file.exists(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))) {
    # if (!file.exists(paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL_", tissue, ".csv"))) {
      print("FunGen preprocessed file does not exist")
      # tmp <- fread(paste("/project_data/data_original/FunGen_pQTL_", tissue, ".FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz", sep=""), verbose=FALSE)
      tmp <- fread(paste("~/gwas_eqtl_colocalization/T2D_OA/data/FunGen_pQTL_", tissue, ".FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz", sep=""), verbose=FALSE)
      
      dict_pQTL[[tissue]] <- tmp[, .(CHR=as.integer(sub(":.*", "", genotype_id)), 
                                     POS=as.integer(sub(".*:", "", genotype_id)), rsID=rep(".", nrow(tmp)),
                                     geneID=sub(".*_", "", phenotype_id), pval, beta=slope, se=slope_se, 
                                     MAF=as.numeric(str_match(INFO, "MAF=\\s*(.*?)\\s*;")[,2]),
                                     N=PERM_num_var, EA=ALT, NEA=REF)]
      
      # write_output(tissue, dict_pQTL, paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))
      write_output(tissue, dict_pQTL, paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL_", tissue, ".csv"))
    }
    
    else{
      print("FunGen preprocessed file already exists")
      dict_pQTL[[tissue]] <- fread(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"), verbose=FALSE)
      # dict_pQTL[[tissue]] <- fread(paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL_", tissue, ".csv"), verbose=FALSE)
    }
  }
}


select_variants <- function(i, dict, credset, mQTL=FALSE) {
  credset.id <- credset[, SNP]
  if (length(credset.id)!=1) {
    dict[[paste0(i, "_coloc")]] <- dict[[i]][ID %in% credset.id]
  }
  else{
    if(mQTL==TRUE) {
      genes <- dict[[i]][ID==credset.id, geneID]
      # dict[[paste0(i, "_coloc")]] <- dict[[i]][geneID %in% genes & POS %between% c(credset[, POS]-1000000, credset[, POS]+1000000)]
      dict[[paste0(i, "_coloc")]] <- dict[[i]][geneID %in% genes & POS %between% c(credset[, POS]-1000000, credset[, POS]+1000000)]
    }
    else {
      dict[[paste0(i, "_coloc")]] <- dict[[i]][POS %between% c(credset[, POS]-1000000, credset[, POS]+1000000)]
    }
  }
}

#----------------------------------------
# Load eQTL data
#----------------------------------------
TISSUE <- data.table(tissue=c("HighGradeCartilage", "LowGradeCartilage"), # "PancreaticIslets"),
                     source=c(rep(c("FunGen"),2))) #, "InsPIRE"))
pQTL <- hash()
print("Loading eQTL data")
TISSUE[, read_pQTL(tissue=tissue, dict_pQTL=pQTL, source=source, delete.indels=TRUE), by=tissue]

###################################################
#--------------- Add ID and rsID for FunGen --------------
###################################################
# for (tissue in TISSUE$tissue){
#   # rsid <- fread(paste0("/project_data/processed_data/Variants2rsID/pQTL.", tissue,".simple.bed")) %>% .[nchar(V5)==1]
#   rsid <- fread(paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL.", tissue,".simple.bed")) %>% .[nchar(V5)==1]
#   rsid <- separate(rsid, "V6", c("V6.1", "V6.2", "V6.3"), sep =",")  %>% .[nchar(V6.1)==1]
#   rsid <- melt(rsid, id.vars=c("V1", "V2", "V3", "V4", "V5"),
#                    measure.vars=c("V6.1", "V6.2", "V6.3"), value.name="V6") %>% .[!is.na(V6)]
#   rsid <- rsid[, variable:=NULL]
#   rsid[, ID:=paste(paste(V2, V3, sep=":"), sort_alleles(V5, V6), sep="_")]
#   # fwrite(rsid, paste0("/project_data/processed_data/Variants2rsID/pQTL.", tissue,".rsid"))
#   fwrite(rsid, paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL.", tissue,".rsid"))
# 
#   # trait <- fread(paste0("/project_data/processed_data/pQTL_", tissue, ".csv"))
#   trait <- fread(paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL_", tissue, ".csv"))
#   trait[, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
#   tst <- merge(trait, rsid, by="ID", all.x=TRUE) %>% .[, rsID:=V4] %>% .[, .(CHR,POS,rsID,geneID,pval,beta,se,MAF,N,EA,NEA,ID)]
#   # fwrite(tst, paste0("/project_data/processed_data/pQTL_", tissue, "_rsid.csv"))
#   fwrite(tst, paste0("~/gwas_eqtl_colocalization/T2D_OA/data/processed_data/pQTL_", tissue, "_rsid.csv"))
# }

OA_trait <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR", "HandOA")
OA_case_con <- c(0.27, 0.187, 0.224, 0.115, 0.078, 0.125, 0.0778, 0.0739)

source("~/gwas_eqtl_colocalization/T2D_OA/scripts/PlotFunctions.R")
GRCh37_Genes <- read.delim(paste0("/project_data/data_original/UCSC_GRCh37_Genes_UniqueList.txt"), stringsAsFactors = FALSE, header = TRUE)
range <- 5e+05
PP=4
# ld.file <- fread(paste0("/project_data/processed_data/LDvariants/chr11_76506572.ld"))
ld.regions <- list("3965689"=4291928, "53501946"=53800954, "150521096"=150537635, "133414622"=133864599, "124468572"=123732769, # "124509177"=123450765,
                   "124509177"=123732769, "51180765"=50788778, "44938870"=45411941, "422144"=653575,"10808687"=9974824) # , "9974824"=10808687, "10808687"=9974824)

for (trait_idx in 1:7) {
  TRAIT <- data.table(trait=c("T2D", OA_trait[trait_idx]), 
                      source=c("DIAGRAM", "GO"), 
                      cas_con=c(0.089956, OA_case_con[trait_idx])) 
  
  ########################################################################
  #------------------------- Read analysis files -------------------------
  ########################################################################
  GWAS <- hash()
  
  for (trait in TRAIT$trait) {
    GWAS[[trait]] <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2],"/GWAS_", trait, "_precoloc_regions.csv"), verbose=FALSE)
  }
  setcolorder(GWAS[["T2D"]], colnames(GWAS[[TRAIT$trait[2]]]))
  
  GWAS_regions <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_regions.csv"))
  rowidx <- order(GWAS_regions[, "CHR"], GWAS_regions[, "POS"])
  GWAS_regions <- GWAS_regions[rowidx]
  
  # Read credible set and select regions with more than 1 SNP in credible set
  credset.dt <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/credible_set2.csv"))
  
  ###################################################
  #--------------- Add ID to GWASes --------------
  ###################################################
  # sort_alleles <- Vectorize(function(x,y) {
  #   paste(sort(c(x, y)), collapse = "_")
  # })
  # 
  # for (trait in TRAIT$trait[1]){
  #   GWAS[[trait]] <- GWAS[[trait]][, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
  #   fwrite(GWAS[[trait]], paste0("/project_data/processed_data/GWAS_", trait,"_id.csv"))
  # }

  ###########################################################################
  #----------- Run HyPrColoc analysis for each colocalized region -----------
  ###########################################################################
  res <- vector(mod ="list", length=length(unique(credset.dt$region)))
  names(res) <- unique(credset.dt$region)
  
  for (reg in unique(credset.dt$region)){
    # ----------- Select only variants included in the credible set -----------
    TRAIT[, select_variants(trait, GWAS, credset.dt[region==reg], mQTL=FALSE), by=trait]
    TISSUE[, select_variants(tissue, pQTL, credset.dt[region==reg], mQTL=TRUE), by=tissue]
    
    for (tissue in TISSUE$tissue){
      coloc.genes <- unique(pQTL[[paste0(tissue, "_coloc")]]$geneID)
      
      for (gene in unique(pQTL[[paste0(tissue, "_coloc")]]$geneID)){
        # ----------- Prepare beta and ses matrices -----------
        print("Preparing input matrices")
        ses <- GWAS[[paste0(TRAIT$trait[1], "_coloc")]][, .(SNPID=ID, rsID=rsID, se=se)]
        setnames(ses, "se", paste0("GWAS_", TRAIT$trait[1]))
        betas <- GWAS[[paste0(TRAIT$trait[1], "_coloc")]][, .(SNPID=ID, rsID=rsID, beta=beta)]
        setnames(betas, "beta", paste0("GWAS_", TRAIT$trait[1]))
        
        for (trait in TRAIT$trait[-1]) {
          ses <- merge(ses, GWAS[[paste0(trait, "_coloc")]][, .(SNPID=ID, se=se)], by="SNPID", allow.cartesian=T)
          setnames(ses, "se", paste0("GWAS_", trait))
          betas <- merge(betas, GWAS[[paste0(trait, "_coloc")]][, .(SNPID=ID, beta=beta)], by="SNPID", allow.cartesian=T)
          setnames(betas, "beta", paste0("GWAS_", trait))
        }
        
        ses <- merge(ses, pQTL[[paste0(tissue, "_coloc")]][geneID==gene, .(SNPID=ID, geneID=paste0("eQTL_", tissue, "_", geneID), se=se)], by="SNPID", allow.cartesian=T)
        ses <- dcast(ses, as.formula(paste(paste(colnames(ses)[c(-ncol(ses), -ncol(ses)+1)], collapse="+"), "~ geneID")), value.var = "se")
        betas <- merge(betas, pQTL[[paste0(tissue, "_coloc")]][geneID==gene, .(SNPID=ID, geneID=paste0("eQTL_", tissue, "_", geneID), beta=beta)], by="SNPID", allow.cartesian=T)
        betas <- dcast(betas, as.formula(paste(paste(colnames(betas)[c(-ncol(betas), -ncol(betas)+1)], collapse="+"), "~ geneID")), value.var = "beta")
        
        betas <- as.data.table(betas)
        ses <- as.data.table(ses)
        
        # ----------- Run HyPrColoc -----------
        print("Colocalization analysis using HyPrColoc")
        id <- betas$SNPID
        rsid <- betas$rsID
        betas_mat <- as.matrix(betas[, c('SNPID', 'rsID'):=NULL])
        rownames(betas_mat) <- id
        ses_mat <- as.matrix(ses[, c('SNPID', 'rsID'):=NULL])
        rownames(ses_mat) <- id
        binary.traits = c(rep(1,nrow(TRAIT)), rep(0,ncol(betas_mat)-nrow(TRAIT)))
        traits <- colnames(betas_mat)
        
        betas_mat <- na.omit(betas_mat)
        ses_mat <- na.omit(ses_mat)
        
        if (nrow(betas_mat)>1){
          res[[reg]] <- hyprcoloc(betas_mat, ses_mat, trait.names=traits, snp.id=id, bb.alg=FALSE,
                                  binary.outcomes=binary.traits, prior.1=1e-10, prior.2=0.7, snpscores=T)
          print(res[[reg]])
          if(res[[reg]]$results$posterior_prob>0.6){
            
            sel.snp <- res[[reg]]$results$candidate_snp[1]
            sel.chr <- as.integer(gsub(":.*", "", sel.snp))
            sel.pos <- as.integer(gsub("_.*", "", gsub(".*:", "", sel.snp)))
            sel.rsid <- GWAS[[TRAIT$trait[2]]][ID==sel.snp, rsID][1]
            
            ld.pos <- GWAS_regions[as.numeric(reg)]$POS
            ld.file <- fread(paste0("/project_data/processed_data/LDvariants/chr", sel.chr , "_", ifelse(ld.pos %in% names(ld.regions),
                                                                                                         as.numeric(unname(ld.regions[match(ld.pos, names(ld.regions))])),
                                                                                                         ld.pos), ".ld"))[SNP_A==sel.rsid | SNP_B==sel.rsid] %>% .[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
            
            # Get data
            OA.data.region <- GWAS[[TRAIT$trait[2]]][CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range) & rsID!="",
                                                     .(ID=ID, CHR=CHR, SNP=rsID, P=pval, BP=POS, logP=-log10(as.numeric(pval)))]
            OA.data.region <- OA.data.region[!duplicated(OA.data.region$SNP)]
            
            T2D.data.region <- GWAS[[TRAIT$trait[1]]][CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range),
                                                      .(ID=ID, CHR=CHR, SNP=rsID, P=pval, BP=POS, logP=-log10(as.numeric(pval)))] %>% .[logP==Inf, logP:=320]
            T2D.data.region <- unique(T2D.data.region, by="SNP")
            
            pQTL.data.region <- pQTL[[tissue]][CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range) & geneID==gene  & rsID!="",
                                               .(ID=ID, CHR=CHR, SNP=rsID, geneID, P=pval, BP=POS, logP=-log10(as.numeric(pval)))]
            
            # Get plot
            data.lst <- list(OA.data.region, T2D.data.region, pQTL.data.region)
            names(data.lst) <-  c(TRAIT$trait[2], "T2D", paste(tissue, gene, sep="_"))
            locus.zoom(data = data.lst,
                       # region = c(sel.chr, sel.pos-range/2, sel.pos+range/2),
                       offset_bp = range,                                            
                       genes.data = GRCh37_Genes,                                       
                       file.name = paste0(TRAIT$trait[2], "_pQTL_", reg, "_", tissue, "_", gene, ".png"),
                       snp=sel.rsid, ignore.lead=TRUE,
                       ld.file=ld.file,
                       pp=paste0("PP",PP),
                       pp.value=round(res[[reg]]$results$posterior_prob[1], digits=3),
                       nplots=TRUE)
          }
        }
      }
    }
  }
}
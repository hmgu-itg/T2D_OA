library(data.table)
library(stringr)
library(coloc)
library(hyprcoloc)
library(hash)
library(dplyr)
library(tidyr)

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

source("/project_data/overlap_T2D_OA/scripts/PlotFunctions.R")
GRCh37_Genes <- read.delim(paste0("/storage/hmgu/projects/OA_T2D/data_original/UCSC_GRCh37_Genes_UniqueList.txt"), stringsAsFactors = FALSE, header = TRUE)
range <- 5e+05
PP=4
# ld.regions <- list("3965689"=4291928, "53501946"=53800954, "150521096"=150537635, "133414622"=133864599, "124468572"=123732769, # "124509177"=123450765,
#                    "124509177"=123732769, "51180765"=50788778, "44938870"=45411941, "422144"=653575) # , "9974824"=10808687, "10808687"=9974824)

ld.regions <- list("219643649"=219584164, "219748818"=219584164, "150521096"=150537635, "219644224"=219584164,
                   "219741820"=219584164, "53501946"=53800954, "123450765"=123732769, "124509177"=123450765,
                   "124468572"=123450765, "9974824"=10808687, "114758349"=114699835, "53800200"=53800954, 
                   "3965689"=4291928, "44938870"=45411941, "133414622"=133864599) # 123732769

TISSUE <- data.table(tissue=c("HighGradeCartilage", "LowGradeCartilage", "Synovium", "PancreaticIslets"),
                     source=c(rep(c("FunGen"),1), "InsPIRE"))
eQTL <- hash()
for (tissue in TISSUE$tissue) {
  eQTL[[tissue]] <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/eQTL_", tissue, ".csv"), verbose=FALSE)
}

OA_trait <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR", "HandOA")
OA_case_con <- c(0.27, 0.187, 0.224, 0.115, 0.078, 0.125, 0.0778, 0.0739)

res.dt <- data.table()

for (trait_idx in c(6,7)) {  # seq(1,7)
  # trait_idx=4
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
  credset.dt <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/credible_set.csv"))

  ###########################################################################
  #----------- Run HyPrColoc analysis for each colocalized region -----------
  ###########################################################################
  # res <- vector(mod ="list", length=length(unique(credset.dt$region)))
  # names(res) <- unique(credset.dt$region)
  # 
  # sim.mat <- vector(mod ="list", length=length(unique(credset.dt$region)))
  # names(sim.mat) <- unique(credset.dt$region)
  
  for (reg in unique(credset.dt$region)){
    # ----------- Select only variants included in the credible set -----------
    TRAIT[, select_variants(trait, GWAS, credset.dt[region==reg], mQTL=FALSE), by=trait]
    TISSUE[, select_variants(tissue, eQTL, credset.dt[region==reg], mQTL=TRUE), by=tissue]

    for (tissue in TISSUE$tissue){
      # eQTL[[paste0(tissue, "_coloc")]] <- eQTL[[tissue]][ID %in% credset.dt[region==reg, SNP]]
      
      coloc.genes <- unique(eQTL[[paste0(tissue, "_coloc")]]$geneID)
      
      for (gene in unique(eQTL[[paste0(tissue, "_coloc")]]$geneID)){
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
        
        ses <- merge(ses, eQTL[[paste0(tissue, "_coloc")]][geneID==gene, .(SNPID=ID, geneID=paste0("eQTL_", tissue, "_", geneID), se=se)], by="SNPID", allow.cartesian=T)
        ses <- dcast(ses, as.formula(paste(paste(colnames(ses)[c(-ncol(ses), -ncol(ses)+1)], collapse="+"), "~ geneID")), value.var = "se")
        betas <- merge(betas, eQTL[[paste0(tissue, "_coloc")]][geneID==gene, .(SNPID=ID, geneID=paste0("eQTL_", tissue, "_", geneID), beta=beta)], by="SNPID", allow.cartesian=T)
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
          res <- hyprcoloc(betas_mat, ses_mat, trait.names=traits, snp.id=id, bb.alg=FALSE,
                           binary.outcomes=binary.traits, prior.1=1e-4, prior.c=0.7) #, snpscores=T)
          # print(res[[reg]])
          # res.sen = sensitivity.plot(betas_mat, ses_mat, trait.names = traits, snp.id=id, bb.alg=FALSE,
          #                            # reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9),
          #                            prior.c = c(0.1, 0.3, 0.5, 0.7, 0.9), equal.thresholds = FALSE, similarity.matrix = TRUE,
          #                            binary.outcomes=binary.traits, prior.1=1e-10)
          # sim.mat = res.sen[[2]]
          # sim.mat
          
          if(res$results$posterior_prob>0.7){
            res.dt <- rbind(res.dt, as.data.table(res$results) %>% .[, region:=reg])
            sel.snp <- res$results$candidate_snp[1]
            sel.chr <- as.integer(gsub(":.*", "", sel.snp))
            sel.pos <- as.integer(gsub("_.*", "", gsub(".*:", "", sel.snp)))
            sel.rsid <- GWAS[[TRAIT$trait[2]]][ID==sel.snp, rsID][1]
            
            ld.pos <- GWAS_regions[as.numeric(reg)]$POS
            ld.file <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/LDvariants/chr", sel.chr , "_", ifelse(ld.pos %in% names(ld.regions), as.numeric(unname(ld.regions[match(ld.pos, names(ld.regions))])), ld.pos), ".ld"))[SNP_A==sel.rsid | SNP_B==sel.rsid]
            ld.file <- ld.file[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
            
            # Get data
            OA.data.region <- GWAS[[TRAIT$trait[2]]][CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range) & rsID!="", .(ID=ID, CHR=CHR, SNP=rsID, P=pval, BP=POS, logP=-log10(as.numeric(pval)))]
            OA.data.region <- OA.data.region[!duplicated(OA.data.region$SNP)]
            
            T2D.data.region <- GWAS[[TRAIT$trait[1]]][CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range), .(ID=ID, CHR=CHR, SNP=rsID, P=pval, BP=POS, logP=-log10(as.numeric(pval)))] 
            T2D.data.region <- T2D.data.region[logP==Inf, logP:=320]            
            T2D.data.region <- unique(T2D.data.region, by="SNP")
            
            eQTL.data.region <- eQTL[[tissue]][CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range) & geneID==gene & rsID!="", .(ID=ID, CHR=CHR, SNP=rsID, geneID, P=pval, BP=POS, logP=-log10(as.numeric(pval)))]
            
            # Get plot
            data.lst <- list(OA.data.region, T2D.data.region, eQTL.data.region)
            names(data.lst) <-  c(TRAIT$trait[2], "T2D", paste(tissue, gene, sep="_"))
            locus.zoom(data = data.lst,
                       # region = c(sel.chr, sel.pos-range/2, sel.pos+range/2),
                       offset_bp = range,                                            
                       genes.data = GRCh37_Genes,                                       
                       file.name = paste0("/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/plots/eQTL/", TRAIT$trait[2], "_eQTL_", reg, "_", tissue, "_", gene, ".png"),
                       snp=sel.rsid, ignore.lead=TRUE,
                       ld.file=ld.file,
                       pp=paste0("PP",PP),
                       pp.value=round(res$results$posterior_prob[1], digits=3),
                       nplots=TRUE)
          }
        }
      }
    }
  }
}

fwrite(res.dt, "/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/gwas_eqtl_hyprcoloc.csv")

library(biomaRt)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

dt <- fread("/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/gwas_eqtl_hyprcoloc.csv")
dt[, `:=` (phen1=strsplit(unlist(strsplit(traits, split=","))[1], split="_")[[1]][2], 
           phen2=strsplit(unlist(strsplit(traits, split=","))[2], split="_")[[1]][2], 
           tissue=strsplit(unlist(strsplit(traits, split=","))[3], split="_")[[1]][2],
           gene=strsplit(unlist(strsplit(traits, split=","))[3], split="_")[[1]][3]), by=seq_len(nrow(dt))]
dt[, `:=`(traits=NULL, disease=ifelse(tissue=="PancreaticIslets", "T2D", "OA"))]
dt[, gene.name:=getBM(attributes='hgnc_symbol', filters='ensembl_gene_id', values=gene, mart=ensembl), by=seq_len(nrow(res.dt2))]

# TO DO: get rsid from chr:pos 
# https://support.bioconductor.org/p/133105/

fwrite(dt, "/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/gwas_eqtl_hyprcoloc.csv")
fwrite(dt[posterior_prob>=0.8], "/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/gwas_eqtl_hyprcoloc_08.csv")

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

####################################################
# --------------- Add ID to eQTL PI --------------
####################################################
# eQTL[["PancreaticIslets"]] <- eQTL[["PancreaticIslets"]][, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
# fwrite(eQTL[["PancreaticIslets"]], paste0("/project_data/processed_data/eQTL_PancreaticIslets_id.csv"))

###################################################
#--------------- Add rsID to FunGen --------------
###################################################
# for (tissue in TISSUE$tissue[c(2,3)]){
#   rsid <- fread(paste0("/project_data/processed_data/Variants2rsID/", tissue,".simple.bed")) %>% .[nchar(V5)==1]
#   rsid <- separate(rsid, "V6", c("V6.1", "V6.2", "V6.3"), sep =",")  %>% .[nchar(V6.1)==1]
#   rsid <- melt(rsid, id.vars=c("V1", "V2", "V3", "V4", "V5"),
#                    measure.vars=c("V6.1", "V6.2", "V6.3"), value.name="V6") %>% .[!is.na(V6)]
#   rsid <- rsid[, variable:=NULL]
#   rsid[, ID:=paste(paste(V2, V3, sep=":"), sort_alleles(V5, V6), sep="_")]
#   fwrite(rsid, paste0("/project_data/processed_data/Variants2rsID/", tissue,".rsid"))
#   
#   trait <- fread(paste0("/project_data/processed_data/eQTL_", tissue, ".csv"))
#   trait[, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
#   tst <- merge(trait, rsid, by="ID", all.x=TRUE) %>% .[, rsID:=V4] %>% .[, .(CHR,POS,rsID,geneID,pval,beta,se,MAF,N,EA,NEA,ID)]
#   fwrite(tst, paste0("/project_data/processed_data/eQTL_", tissue, "_rsid.csv"))
# }


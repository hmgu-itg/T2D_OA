# setwd("~/gwas_eqtl_colocalization")

library(data.table)
library(stringr)
library(coloc)
library(hash)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

OA_trait <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR") #, "HandOA")
OA_case_con <- c(0.27, 0.187, 0.224, 0.115, 0.078, 0.125, 0.0778) #, 0.0739)

for (trait_idx in seq(2,7)) {
  # trait_idx=5
  TRAIT <- data.table(trait=c("T2D", OA_trait[trait_idx]),
                      source=c("DIAGRAM", "GO"), 
                      cas_con=c(0.089956, OA_case_con[trait_idx])) 
  
  ########################################################################
  #------------------------- Read analysis files -------------------------
  ########################################################################
  GWAScoloc <- hash()
  
  for (trait in TRAIT$trait) {
    # if(!file.exists(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_", trait, "_precoloc_regions.csv"))){
    #   source("~/gwas_eqtl_colocalization/scripts/GWASPreprocess.R")
    # }
    GWAScoloc[[trait]] <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_", trait, "_precoloc_regions.csv"), verbose=FALSE)
    # GWAScoloc[[trait]] <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/GWAS", TRAIT$trait[2], "/GWAS_", trait, "_precoloc_regions.csv"), verbose=FALSE)
  }
  
  GWAS_regions <- fread(paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/GWAS_regions.csv"))
  # GWAS_regions <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/GWAS", TRAIT$trait[2], "/GWAS_regions.csv"))
  rowidx <- order(GWAS_regions[, "CHR"], GWAS_regions[, "POS"])
  GWAS_regions <- GWAS_regions[rowidx]
  
  data <- merge(GWAScoloc[[TRAIT$trait[1]]], GWAScoloc[[TRAIT$trait[2]]], by=c("ID", "CHR", "POS"), allow.cartesian = T)
  data[rsID.x=="", rsID.x:=rsID.y]
  data[rsID.y==".", rsID.y:=rsID.x]
  
  ########################################################################
  #----------------------- Run coloc.abf analysis ------------------------
  # https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html
  ########################################################################
  print("Colocalization analysis using coloc.abf")
  
  #------------------ coloc.abf locus-wise ----------------------
  res <- vector(mod ="list", length=nrow(GWAS_regions))
  names(res) <- seq(1:nrow(GWAS_regions))
  credset <- vector(mod ="list", length=nrow(GWAS_regions))
  names(credset) <- seq(1:nrow(GWAS_regions))
  
  for (i in 1:nrow(GWAS_regions)){
    data_tmp <- data[CHR==GWAS_regions[i]$CHR & POS %between%
                       c(GWAS_regions[i]$start, GWAS_regions[i]$end)]
    res[[i]] <- coloc.abf(dataset1=list(beta=data_tmp$beta.x, varbeta=(data_tmp$se.x)^2, type="cc", s=TRAIT$cas_con[1], N=nrow(GWAScoloc[[TRAIT$trait[1]]])),
                          dataset2=list(beta=data_tmp$beta.y, varbeta=(data_tmp$se.y)^2, type="cc", s=TRAIT$cas_con[2], N=nrow(GWAScoloc[[TRAIT$trait[2]]])))
  
    # res.abf[[i]] <- coloc.susie(dataset1=list(beta=data_tmp$beta.x, varbeta=(data_tmp$se.x)^2, type="cc", s=TRAIT$cas_con[1], N=nrow(GWAScoloc[[TRAIT$trait[1]]])),
    #                             dataset2=list(beta=data_tmp$beta.y, varbeta=(data_tmp$se.y)^2, type="cc", s=TRAIT$cas_con[2], N=nrow(GWAScoloc[[TRAIT$trait[2]]])))
    
    # Get 95% credible set for regions where PP4 > 0.8
    # https://chr1swallace.github.io/coloc/articles/a03_enumeration.html
    if (res[[i]]$summary[6]>0.8){
      o <- order(res[[i]]$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(res[[i]]$results$SNP.PP.H4[o])
      w <- which(cs > 0.95)[1]
      credset[[i]] <- data_tmp[as.numeric(gsub("(SNP.)","",res[[i]]$results[o,][1:w,]$snp)), ID]
    }
  }
  
  # output credible set as data.table for GWAS eQTL colocalization
  credset.dt <- data.table()
  length(Filter(Negate(is.null), credset))
  for(i in 1:length(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))])){
    reg <- names(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][i])
    credible.set <- Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][[i]]
    region.n <- length(credible.set)
    
    o <- order(res[[reg]]$results$SNP.PP.H4, decreasing=TRUE)
    
    credset.dt <- rbind(credset.dt, data.table(region=rep(reg, region.n), SNP=credible.set, 
                                               rsID=data[, rsID.x[match(credible.set, ID)]], 
                                               CHR=data[, CHR[match(credible.set, ID)]], 
                                               POS=data[, POS[match(credible.set, ID)]],
                                               PP4=rep(round(unname(res[[reg]]$summary[paste0("PP.H4.abf")]), digits=3), region.n), 
                                               SNP.PP4=res[[reg]]$results$SNP.PP.H4[o] %>% .[1:region.n]))
  }
  fwrite(credset.dt, paste0("/project_data/processed_data/GWAS", TRAIT$trait[2], "/credible_set.csv"))
  
  ###################################################################
  # --------------- Get LD between SNPs in each region --------------
  ###################################################################
  # for (i in 1:length(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))])){
    # region <- as.numeric(names(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][i]))
  #   region.data <- data[CHR==GWAS_regions[region]$CHR & POS %between% c(GWAS_regions[region]$POS - range, GWAS_regions[region]$POS + range) & rsID.x!=".", rsID.x]
  #   chr <- data[rsID.x %in% region.data, unique(CHR)]
  #   fwrite(as.list(region.data), paste0("/project_data/processed_data/LDvariants/", TRAIT$trait[2], "/region", region, ".chr", chr), sep="\n")
  # }
  
  ###############################################################
  # --------------- Plot regional association plot --------------
  ###############################################################
  source("/project_data/overlap_T2D_OA/scripts/PlotFunctions.R")
  GRCh37_Genes <- read.delim(paste0("/storage/hmgu/projects/OA_T2D/data_original/UCSC_GRCh37_Genes_UniqueList.txt"), stringsAsFactors = FALSE, header = TRUE)
  range <- 5e+05
  PP=4
  
  ld.regions <- list("219643649"=219584164, "219748818"=219584164, "150521096"=150537635, "219644224"=219584164,
                     "219741820"=219584164, "53501946"=53800954, "123450765"=123732769, "124509177"=123450765,
                     "124468572"=123450765, "9974824"=10808687, "114758349"=114699835, "53800200"=53800954, 
                     "3965689"=4291928, "44938870"=45411941, "133414622"=133864599) # 123732769
  
  for (i in 1:length(unique(Filter(Negate(is.null), credset)))){
    i=2
    region.nr <- names(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][i])
    sel.snp <- Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][[i]][1]
    sel.chr <- data[ID==sel.snp, CHR]
    sel.pos <- data[ID==sel.snp, POS]
    data.region <- data %>% filter(CHR==sel.chr, between(POS, sel.pos-range, sel.pos+range))
    
    if(data[ID==sel.snp, rsID.y]!=".") sel.rsid=data[ID==sel.snp, rsID.y]
    else if(data[ID==sel.snp, rsID.x]!=".") sel.rsid=data[ID==sel.snp, rsID.x]
    else print(paste0("ERROR: no rsID available for the lead SNP ", sel.snp))
    
    ld.pos <- GWAS_regions[as.numeric(region.nr)]$POS
    ld.file <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/LDvariants/chr", sel.chr , "_", 
                            ifelse(ld.pos %in% names(ld.regions), as.numeric(unname(ld.regions[match(ld.pos, names(ld.regions))])), ld.pos), ".ld"))[SNP_A==sel.rsid | SNP_B==sel.rsid] %>% .[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
    # ld.file <- fread("/storage/hmgu/projects/OA_T2D/processed_data/LDvariants/chr14_94844947.ld")[SNP_A==sel.rsid | SNP_B==sel.rsid] %>% .[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
    
    if(length(i[-1])!=0){
      credible.set <- as.vector(data[ID %in% Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][[i]][-1] & rsID.y!=".", rsID.y])
    }else{
      credible.set=NA
    }
    
    # Get data
    OA.data.region <- data[CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range) & rsID.x!=".",
                           .(CHR=CHR, SNP=rsID.x, P=pval.y, BP=POS, logP=-log10(as.numeric(pval.y)))]
    T2D.data.region <- data[CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range),
                            .(CHR=CHR, SNP=rsID.x, P=pval.x, BP=POS, logP=-log10(as.numeric(pval.x)))] %>% .[logP==Inf, logP:=320]
    T2D.data.region <- unique(T2D.data.region, by="SNP")
    
    # Regional association plot
    data.lst <- list(OA.data.region, T2D.data.region)
    names(data.lst) <-  c(TRAIT$trait[2], "T2D")
    locus.zoom(data = data.lst, 
               # region = c(sel.chr, sel.pos-range/2, sel.pos+range/2),
               offset_bp = range,                                            
               genes.data = GRCh37_Genes,                                       
               file.name = paste0("/project_data/ColocPlots/", TRAIT$trait[2], "_T2D_", region.nr, "_beta.png"),
               secondary.snp = credible.set,    
               snp=sel.rsid, ignore.lead=TRUE,
               ld.file=ld.file,
               pp=paste0("PP",PP),
               pp.value=round(unname(res[[region.nr]]$summary[paste0("PP.H",PP,".abf")]), digits=3),
               nplots=TRUE)
    
    # Plot PP4 of credible set
    if (res[[region.nr]]$summary[6]>0.8){
      o <- order(res[[region.nr]]$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(res[[region.nr]]$results$SNP.PP.H4[o])
      w <- which(cs > 0.95)[1]
      PP4 <- res[[region.nr]]$results$SNP.PP.H4[o] %>% .[1:length(credset[[region.nr]])]}
    
    plot.dt <- data.table(SNP=data[, rsID.x[match(credset[[region.nr]], ID)]],
                          POS=data[, POS[match(credset[[region.nr]], ID)]], 
                          PP4=PP4)
    
    locus.zoom(data = plot.dt[, .(CHR=sel.chr, SNP=SNP, PP4=PP4, BP=POS, logP=PP4)],   
               region = c(sel.chr, min(plot.dt$POS)-100, max(plot.dt$POS)+100),                                        
               genes.data = GRCh37_Genes,                                       
               file.name = paste0("/project_data/ColocPlots/", TRAIT$trait[2], "_PP4_", region.nr, "_beta.png"),         
               # secondary.label = TRUE,
               snp=sel.rsid, ignore.lead=TRUE,
               ld.file=ld.file,
               sig.type="PP4")
  }
}


##############################################################################
# -------------------- Compare credset from run1 and run2 -------------------
##############################################################################
for (trait in c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR")){
  trait <- "TKR"
  # credset1 <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/GWAS", trait,"/credible_set2.csv")) %>% .[PP4>0.8]
  # credset1 <- credset1[region!=397]
  credset2 <- fread(paste0("/project_data/processed_data/GWAS", trait, "/credible_set.csv"))
  all.equal(credset1[, .(SNP, CHR, POS, PP4, SNP.PP4)], credset2[, .(SNP, CHR, POS, PP4, SNP.PP4)])
  all.equal(credset1[, .(SNP, CHR, POS, SNP.PP4)], credset2[, .(SNP, CHR, POS, SNP.PP4)])
  all.equal(credset1[, .(SNP, CHR, POS)], credset2[, .(SNP, CHR, POS)])
}


##############################################################################
# --------------- Get info about lead variants from each region --------------
##############################################################################
for (i in names(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))])){
  print(paste0("REGION: ", i))
  credset[[i]]
  print(paste0("Length of credible set: ", length(credset[[i]])))
  print(paste0("rsID: ", data[ID==credset[[i]][2], rsID.x]))
  
  o <- order(res[[i]]$results$SNP.PP.H4,decreasing=TRUE)
  cs <- cumsum(res[[i]]$results$SNP.PP.H4[o])
  w <- which(cs > 0.95)[2]
  res[[i]]$results[o,][1:w,]$snp
  top.snp <- res[[i]]$results[o,][1:w,]$snp[2]
  print(paste0("Top SNP: ", top.snp))
  
  # data_tmp <- data[CHR==GWAS_regions[i]$CHR & POS %between%
  #                  c(GWAS_regions[i]$start, GWAS_regions[i]$end)]
  # data_tmp[as.numeric(gsub("(SNP.)","", top.snp))]
  print(paste0("SNP.PP.H4: ", as.data.table(res[[i]]$results) %>% .[snp==top.snp] %>% .[,SNP.PP.H4]))
  print("-----------------------------------")
}



##############################################################
# --------------------- Test coloc.susie ---------------------
# https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
##############################################################
source("~/gwas_eqtl_colocalization/scripts/susie.R")
credset <- credset.knee
res <- res.knee
data <- data.knee
GWAS_regions <- GWAS_regions.knee
i <- 105
data_tmp <- data[CHR==GWAS_regions[i]$CHR & POS %between%
                   c(GWAS_regions[i]$start, GWAS_regions[i]$end)]

# Get LD matrix
# fwrite(as.list(D1$snp), "/project_data/processed_data/LDvariants/D1.snps", sep="\n")
ld.test <- acast(ld.file, SNP_A~SNP_B, value.var="R2",fun.aggregate = list)
rownames(ld.test) <- ld.test$SNP_A
ld.test$SNP_A <- NULL
ld.mat <- as.matrix(ld.test)
for (i in 1:nrow(ld.mat)){
  for (j in 1:ncol(ld.mat)){
    if(identical(ld.mat[[i,j]], numeric(0))) {
      ld.mat[i,j]=ld.mat[j,i]}
    if(identical(ld.mat[[j,i]], numeric(0))){
      ld.mat[j,i]=ld.mat[j,i]}
  }
}

diag(ld.mat) <- 1

D1 <- data_tmp[, .(ID, CHR, position=POS, pvalues=pval.x, beta=beta.x, varbeta=(se.x)^2, MAF=MAF.x, N=N.x, Ncases=N.x, EA=EA.x, NEA=NEA.x, snp=rsID.x)]
D1 <- unique(D1, by="snp") %>% .[MAF==0, MAF:=0.000001]
D1 <- as.list(D1) 
D1$type="cc"
names(D1$beta) <- D1$snp
names(D1$varbeta) <- D1$snp
D1$LD <- ld.mat
check_dataset(D1)
S1 <- runsusie(D1)

D2 <- data_tmp[, .(ID, CHR, type="cc", position=POS, pvalues=pval.y, beta=beta.y, varbeta=(se.y)^2, MAF=MAF.y, N=N.y, Ncases=N.y, EA=EA.y, NEA=NEA.y, snp=rsID.y)]
D2 <- unique(D2, by="snp")
D2$type="cc"
names(D2$beta) <- D2$snp
names(D2$varbeta) <- D2$snp
D2$LD <- ld.mat
check_dataset(D2)
S2 <- runsusie(D2)



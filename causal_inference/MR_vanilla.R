library(TwoSampleMR)
library(data.table)
library(dplyr)
library(hash)
library(ieugwasr)
library(stringr)
library(ggplot2)
library(MendelianRandomization)

################################################
#------------------- Functions -----------------
################################################
sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
})

local.clump <- function(data, pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure), #id=data$id.exposure),
                         plink_bin=genetics.binaRies::get_plink_binary(),
                         bfile="/project_data/data_original/1kg_v3/EUR",
                         clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

get.filename <- function(trait, subset="all"){
  if(trait=="OA.T2D"){
    dt <- fread("/project_data/GWAS_OA_T2D/regenie/ukb_assoc/out_info03/OA.T2D.assoc.txt.gz") %>% 
      .[, .(rsID=ID,CHR=CHROM, POS=GENPOS, pval=PVAL, beta=BETA, se=SE, EAF=A1FREQ, N=31472, Ncases=6799, 
            EA=ALLELE1, NEA=ALLELE0, Ncontrols=24673)]
  }
  else if(trait=="T2D" & subset=="all"){
    dt <- fread("/home/nalu/T2D_OA/T2D_rsid.csv") %>%
      .[, .(rsID, CHR, POS, pval, beta, se, EAF, EA, NEA, Ncases=74124, N=898130, Ncontrols=824006)]
  }
  else if(trait=="T2D" & subset=="UKBB"){
    dt <- fread("/storage/hmgu/projects/OA_T2D/processed_data/T2D_UKBB_rsid.csv") %>% 
      .[, .(rsID, CHR, POS, pval, beta, se, EAF, EA, NEA, Ncases=19119, N=442817, Ncontrols=423698)]
  }
  else if(trait=="T2D" & subset=="noUKBB"){
    dt <- fread("/project_data/data/T2D.noUKBB.gz") %>% .[, .(rsID=SNP, CHR, POS, pval=P, beta=BETA, se=StdErr, EAF=Freq1, 
                                                              EA, NEA, Ncases=55005, N=455313, Ncontrols=400308)]
  }
  else if(trait %nin% c("T2D", "OA.T2D") & subset=="all") {
    dt <- fread(paste0("/storage/hmgu/projects/OA_T2D/data_original/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz")) %>%
      .[, .(CHR,POS,beta=BETA,rsID=SNP,se=SE,EAF,EA,NEA,pval=P,Ncases=NCASES,Ncontrols=NCONTROLS,N)] }
  else if(trait %nin% c("T2D", "OA.T2D") & subset=="UKBB"){
    dt <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/", trait, "_UKBB_rsid.csv"))
  }
  else if(trait %nin% c("T2D", "OA.T2D") & subset=="noUKBB") {
    dt <- fread(paste0("/project_data/data/GO_noUKBB/GO.FILTER.GW.", trait, ".noUKBB.13082019.txt.gz")) %>%
      .[, .(CHR,POS,beta=BETA,rsID=SNP,se=SE,EAF,EA,NEA,pval=P,Ncases=NCASES,Ncontrols=NCONTROLS,N)]
  }
  else print("ERROR: Subset not known")
  
  return(dt)
}

get.GO.filename <- function(trait, subset){
  if(subset=="all") {
    filename <- fread(paste0("/project_data/data_original/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz")) %>%
      .[, .(CHR,POS,beta=BETA,rsID=SNP,se=SE,EAF,EA,NEA,pval=P,Ncases=NCASES,Ncontrols=NCONTROLS,N)] }
  else if(subset=="UKBB"){
    filename <- fread(paste0("/project_data/processed_data/", phen, "_UKBB_rsid.csv"))
  }
  else if(subset=="noUKBB") {
    filename <- fread(paste0("/project_data/GO_noUKBB/GO.FILTER.GW.", trait, ".noUKBB.13082019.txt.gz")) %>%
      .[, .(CHR,POS,beta=BETA,rsID=SNP,se=SE,EAF,EA,NEA,pval=P,Ncases=NCASES,Ncontrols=NCONTROLS,N)]
  }
  else print("ERROR: Subset not known")

  return(filename)
}

get.data <- function(file, data.type){
  data <- format_data(type = data.type,
                      dat = file,
                      snp_col = "rsID",
                      beta_col = "beta",
                      se_col = "se",
                      effect_allele_col = "EA",
                      other_allele_col = "NEA",
                      eaf_col = "EAF",
                      pval_col = "pval",
                      chr_col="CHR",
                      pos_col="POS",
                      ncase_col="Ncases",
                      samplesize_col="N",
                      ncontrol_col="Ncontrols")
  return(data)
}

run.TwoSampleMR <- function(data, filename, folder="vanilla"){
  if(nrow(data)==0) print("No SNP in common between exposure and outcome :(")
  else if (nrow(data)==1) {
    print("Only one SNP in common between exposure and outcome, running Wald ratio method")
    res <- mr(data, method_list="mr_wald_ratio")
    res <- generate_odds_ratios(res)
    
    fstat <- (res$b)^2/(res$se)^2
    tmp.dt <- as.data.table(res) %>% .[, `:=`(Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MR")
    res <- mr(data, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    res <- generate_odds_ratios(res)
    
    #------------ Sensitivity analysis
    mr_scatter_plot(res, data)   #scatter plot
    ggsave(paste0("/home/nalu/T2D_OA/causal_inference/", folder, "/TwoSampleMR_", filename, ".jpg"))
    
    res$`Fstat_Inverse variance weighted` <- (res[res$method=="Inverse variance weighted",]$b)^2/(res[res$method=="Inverse variance weighted",]$se)^2
    res$`Fstat_MR Egger` <- (res[res$method=="MR Egger",]$b)^2/(res[res$method=="MR Egger",]$se)^2
    res$`Fstat_Weighted median` <- (res[res$method=="Weighted median",]$b)^2/(res[res$method=="Weighted median",]$se)^2
    het <- mr_heterogeneity(data)   #heterogeneity statistics
    plt <- mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy
    tmp.dt <- as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T)) %>% 
      .[, `:=`(plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))] %>%
      .[method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]
    
    mr_forest_plot(mr_leaveoneout(data))
    ggsave(paste0("/home/nalu/T2D_OA/causal_inference/", folder, "/LeaveOneOut_", filename, ".jpg"))
    mr_forest_plot(mr_singlesnp(data))
    ggsave(paste0("/home/nalu/T2D_OA/causal_inference/", folder, "/SingleSNP_", filename, ".jpg"))
  }
  return(tmp.dt)
}

output.wide.result <- function(dt){
  dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]
  dt <- dcast(dt, exposure + outcome + nsnp ~ method,
              value.var=c("b", "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95", "Q", "Q_df", "Q_pval", "plt.egger_intercept", "plt.pval", "plt.se", "p.adj.fdr"))
  # remove plt from IVW
  dt[ ,`:=`(`plt.egger_intercept_Inverse variance weighted`=NULL, `plt.pval_Inverse variance weighted`=NULL, `plt.se_Inverse variance weighted`=NULL)]
  # remove Q and plt from WM
  dt[ ,`:=`(`plt.egger_intercept_Weighted median`=NULL, `plt.pval_Weighted median`=NULL, `plt.se_Weighted median`=NULL,
            `Q_Weighted median`=NULL, `Q_df_Weighted median`=NULL, `Q_pval_Weighted median`=NULL)]
  # remove Q and plt from Wald ratio
  dt[ ,`:=`(`plt.egger_intercept_Wald ratio`=NULL, `plt.pval_Wald ratio`=NULL, `plt.se_Wald ratio`=NULL,
            `Q_Wald ratio`=NULL, `Q_df_Wald ratio`=NULL, `Q_pval_Wald ratio`=NULL)]
  
  setcolorder(dt, c(colnames(dt)[1:3], 
                    str_subset(colnames(dt), "Inverse variance weighted"),
                    str_subset(colnames(dt), "Wald ratio"),
                    str_subset(colnames(dt), "Weighted median"),
                    str_subset(colnames(dt), "MR Egger")))
  return(dt)
}

wrap.MR <- function(exposure.lst, outcome.lst, ivs=NULL, exp.phen, out.phen, folder, dt, dt.steiger=NULL, output.file){
  # exposure.lst: list of vectors (for each vector: trait & subset) e.g. list(c("KneeOA", "UKBB"), c("T2D", "noUKBB"))
  # outcome.lst: list of vectors (for each vector: trait & subset) e.g. list(c("KneeOA", "UKBB"), c("T2D", "noUKBB"))
  
  for (exp in exposure.lst){
    #----------- Read and clump exposure -----------
    exposure <- get.data(get.filename(trait=exp[1], subset=exp[2]), "exposure")
    exposure$exposure <- paste(exp, collapse=".")
    if(ivs) clump.exposure <- exposure[exposure$SNP %in% ivs,]
    else clump.exposure <- local.clump(data=exposure)
    
    for (out in outcome.lst) {
      #----------- Get outcome -----------
      outcome <- get.data(get.filename(trait=out[1], subset=out[2]), "outcome")
      outcome$outcome <- paste(out, collapse=".")
      
      #----------- Harmonize data -----------
      data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
      
      #----------- Run TwoSampleMR for unfiltered data -----------
      dt <- rbind(dt,run.TwoSampleMR(data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse=".")), folder=folder))
      
      if(!is.null(dt.steiger)) {
        #----------- Apply steiger filtering on data -----------
        filtered.data <- steiger_filtering(data)
        filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
        
        #----------- Run TwoSampleMR for Steiger-filtered data -----------
        dt.steiger <- rbind(dt.steiger,run.TwoSampleMR(filtered.data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse="."), "_Steiger"), folder=folder))
      }
    }
  }
  
  dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
  fwrite(output.wide.result(dt), paste0("/project_data/overlap_T2D_OA/causal_inference/", folder, "/", output.file, "_analysis_MR.csv"))
  
  if(!is.null(dt.steiger)){
    dt.steiger[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
    fwrite(output.wide.result(dt.steiger), paste0("/project_data/overlap_T2D_OA/causal_inference/", folder, "/Steiger_", output.file, "_analysis_MR.csv"))
  }   
}

#################################################################
#---------------------- OA phenotypes --------------------------#
#################################################################
OA.phen <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "THR", "TJR")

##############################################################################
#---------------------- Exposure: T2D; Outcome: OA --------------------------#
##############################################################################
dt <- data.table()
dt.steiger <- data.table()

wrap.MR(exposure.lst=list(c("T2D", "all")), 
        outcome.lst=mapply(append, as.list(c(OA.phen, "OA.T2D")), as.list(rep("all", length(OA.phen)+1)), SIMPLIFY=FALSE), 
        # ivs=T2D.ivs, 
        folder="vanilla",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="T2D_OA")

##############################################################################
#---------------------- Exposure: OA; Outcome: T2D --------------------------#
##############################################################################
dt <- data.table()
dt.steiger <- data.table()

wrap.MR(exposure.lst=mapply(append, as.list(OA.phen), as.list(rep("all", length(OA.phen))), SIMPLIFY=FALSE), 
        outcome.lst=list(c("T2D", "all"), c("OA.T2D", "UKBB")), 
        folder="vanilla",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="OA_T2D")


##############################################################################
#---------------------- Exposure: T2D; Outcome: OA --------------------------#
##############################################################################
# Read and clump T2D IVs (all independent SNPs of GWAS(T2D))
# T2D <- fread("/home/nalu/T2D_OA/T2D_rsid.csv")
# T2D$Ncases <- 74124
# T2D$N <- 898130
# T2D$Ncontrols <- 824006 
# 
# OA.T2D <- fread("/home/nalu/T2D_OA/GWAS_comorbidity/OA.T2D.assoc.txt.gz") %>% 
#   .[, .(rsID=ID,CHR=CHROM, POS=GENPOS, pval=PVAL, beta=BETA, se=SE, 
#         EAF=A1FREQ, N=31472, Ncases=6799, EA=ALLELE1, NEA=ALLELE0, 
#         Ncontrols=24673)]
# 
# exposure <- get.data(T2D, "exposure")
# exposure$exposure <- "T2D"
# clump.exposure <- local.clump(exposure)
# 
# for (phen in c("OA.T2D", OA.phen)) {
#   #----------- Read outcome summary stats -----------
#   if (phen == "OA.T2D") outcome <- get.data(OA.T2D, "outcome")
#   else outcome <- get.data(get.GO.filename(trait=phen, subset="all"), "outcome")
#   outcome$outcome <- phen
# 
#   #----------- Harmonize data -----------
#   data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
# 
#   #----------- Run TwoSampleMR for unfiltered data -----------
#   tmp.dt <- run.TwoSampleMR(data, filename=paste0("T2D_", phen), folder="vanilla")
#   dt <- rbind(dt,tmp.dt)
# 
#   #----------- Apply steiger filtering on data -----------
#   filtered.data <- steiger_filtering(data)
#   filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
# 
#   #----------- Run TwoSampleMR for Steiger-filtered data -----------
#   tmp.dt <- run.TwoSampleMR(filtered.data, filename=paste0("T2D_", phen, "_Steiger"), folder="vanilla")
#   dt.steiger <- rbind(dt.steiger,tmp.dt)
# }
# 
# ##############################################################################
# #---------------------- Exposure: OA; Outcome: T2D --------------------------#
# ##############################################################################
# dt <- data.table()
# dt.steiger <- data.table()
# 
# wrap.MR(exposure.lst=mapply(append, as.list(OA.phen), as.list(rep("all", length(OA.phen))), SIMPLIFY=FALSE), 
#         outcome.lst=list(c("T2D", "all"), c("OA.T2D", "UKBB")), 
#         folder="vanilla",
#         dt=dt,
#         dt.steiger=dt.steiger,
#         output.file="OA_T2D")
# 
# #----------- Read outcome summary stats -----------
# # outcome <- get.data(fread("/home/nalu/T2D_OA/T2D_rsid.csv"), "outcome")
# for (out.phen in c("T2D", "OA.T2D")){
#   outcome <- get.data(get(out.phen), "outcome")
#   outcome$outcome <- out.phen
#   
#   for (phen in OA.phen) {
#     #----------- Read and clump OA IVs (all independent SNPs of GWAS(OA)) -----------
#     exposure <- get.data(get.GO.filename(trait=phen, subset="all"), "exposure")
#     exposure$exposure <- phen
#     clump.exposure <- local.clump(exposure)
#   
#     #----------- Harmonize data -----------
#     data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
#   
#     #----------- Run TwoSampleMR -----------
#     tmp.dt <- run.TwoSampleMR(data, filename=paste0(phen, "_", out.phen), folder="vanilla")
#     dt <- rbind(dt,tmp.dt)
#   
#     #----------- Apply steiger filtering on data -----------
#     filtered.data <- steiger_filtering(data)
#     filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
#   
#     #----------- Run TwoSampleMR -----------
#     tmp.dt <- run.TwoSampleMR(filtered.data, filename=paste0(phen, "_", out.phen, "_Steiger"), folder="vanilla")
#     dt.steiger <- rbind(dt.steiger,tmp.dt)
#   }
# }
# 
# ##############################################################################
# #-------------------- Adjust p-values and write output-----------------------#
# ##############################################################################
# #--------------- Unfiltered data ----------------
# dt[exposure=="OA.T2D" | exposure=="T2D", p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt[outcome=="OA.T2D" | outcome=="T2D", p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt <- output.wide.result(dt)
# fwrite(dt, paste0("/home/nalu/T2D_OA/causal_inference/vanilla/Analysis_MR.csv"))
# 
# #--------------- Steiger-filtered data ----------------
# dt.steiger[exposure=="OA.T2D" | exposure=="T2D", p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt.steiger[outcome=="OA.T2D" | outcome=="T2D", p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt.steiger <- output.wide.result(dt.steiger)
# fwrite(dt.steiger, paste0("/home/nalu/T2D_OA/causal_inference/vanilla/Steiger_Analysis_MR.csv"))

# dt[exposure=="OA.T2D" | exposure=="T2D", `p.adj.fdr_Wald ratio`:=p.adjust(`pval_Wald ratio`, method = "fdr")]
# dt[outcome=="OA.T2D" | outcome=="T2D", `p.adj.fdr_Wald ratio`:=p.adjust(`pval_Wald ratio`, method = "fdr")]
# dt[exposure=="OA.T2D" | exposure=="T2D", `p.adj.fdr_MR Egger`:=p.adjust(`pval_MR Egger`, method = "fdr")]
# dt[outcome=="OA.T2D" | outcome=="T2D", `p.adj.fdr_MR Egger`:=p.adjust(`pval_MR Egger`, method = "fdr")]
# dt[exposure=="OA.T2D" | exposure=="T2D", `p.adj.fdr_Weighted median`:=p.adjust(`pval_Weighted median`, method = "fdr")]
# dt[outcome=="OA.T2D" | outcome=="T2D", `p.adj.fdr_Weighted median`:=p.adjust(`pval_Weighted median`, method = "fdr")]
# dt[exposure=="OA.T2D" | exposure=="T2D", `p.adj.fdr_Inverse variance weighted`:=p.adjust(`pval_Inverse variance weighted`, method = "fdr")]
# dt[outcome=="OA.T2D" | outcome=="T2D", `p.adj.fdr_Inverse variance weighted`:=p.adjust(`pval_Inverse variance weighted`, method = "fdr")]









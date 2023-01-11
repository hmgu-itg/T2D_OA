library(TwoSampleMR)
library(data.table)
library(dplyr)
library(hash)
library(ieugwasr)
library(stringr)
library(ggplot2)
library(MendelianRandomization)
library(Hmisc)

################################################
#------------------- Functions -----------------
################################################
sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
}) 

local.clump <- function(data, pval=5e-08){
  clumped.dt <- ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure),# id=data$ID),
                         plink_bin=genetics.binaRies::get_plink_binary(),
                         bfile="/reference_data/1kG/EUR_1kg_v3/EUR",
                         clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

get.GO.filename <- function(trait, subset){
  if(subset=="all") {
    # filename <- fread(paste0("/project_data/data_original/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz")) %>%
    filename <- fread(paste0("/storage/hmgu/projects/OA_T2D/data_original/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz")) %>%
      .[, .(CHR,POS,beta=BETA,rsID=SNP,se=SE,EAF,EA,NEA,pval=P,Ncases=NCASES,Ncontrols=NCONTROLS,N)] }
  else if(subset=="UKBB"){
    # filename <- fread(paste0("/project_data/processed_data/", phen, "_UKBB_rsid.csv"))
    filename <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/", trait, "_UKBB_rsid.csv"))
  }
  else if(subset=="noUKBB") {
    # filename <- fread(paste0("/project_data/GO_noUKBB/GO.FILTER.GW.", trait, ".noUKBB.13082019.txt.gz")) %>%
    filename <- fread(paste0("/storage/hmgu/projects/OA_T2D/GO_noUKBB/GO.FILTER.GW.", trait, ".noUKBB.13082019.txt.gz")) %>%
      .[, .(CHR,POS,beta=BETA,rsID=SNP,se=SE,EAF,EA,NEA,pval=P,Ncases=NCASES,Ncontrols=NCONTROLS,N)]
  }
  else print("ERROR: Subset not known")
  
  return(filename)
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

run.TwoSampleMR <- function(data, filename, folder){
  if(nrow(data)==0) {
    print("No SNP in common between exposure and outcome :(")
    tmp.dt <- data.table()
  }
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
    ggsave(paste0("/project_data/MendelianRandomization/", folder, "/TwoSampleMR_", filename, ".jpg"))
    
    # res$`Fstat_Inverse variance weighted` <- (res[res$method=="Inverse variance weighted",]$b)^2/(res[res$method=="Inverse variance weighted",]$se)^2
    # res$`Fstat_MR Egger` <- (res[res$method=="MR Egger",]$b)^2/(res[res$method=="MR Egger",]$se)^2
    # res$`Fstat_Weighted median` <- (res[res$method=="Weighted median",]$b)^2/(res[res$method=="Weighted median",]$se)^2
    het <- mr_heterogeneity(data)   #heterogeneity statistics
    plt <- mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy
    tmp.dt <- as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T)) %>% 
      .[, `:=`(plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))] %>%
      .[method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]

    mr_forest_plot(mr_leaveoneout(data))
    ggsave(paste0("/project_data/MendelianRandomization/", folder, "/LeaveOneOut_", filename, ".jpg"))
    mr_forest_plot(mr_singlesnp(data))
    ggsave(paste0("/project_data/MendelianRandomization/", folder, "/SingleSNP_", filename, ".jpg"))
  }
  return(tmp.dt)
}

wrap.MR <- function(exposure.lst, outcome.lst, ivs=NULL, exp.phen, out.phen, folder, dt, dt.steiger=NULL, output.file){
  # exposure.lst: list of vectors (for each vector: trait & subset) e.g. list(c("KneeOA", "UKBB"), c("T2D", "noUKBB"))
  # outcome.lst: list of vectors (for each vector: trait & subset) e.g. list(c("KneeOA", "UKBB"), c("T2D", "noUKBB"))
  
  for (exp in exposure.lst){
    #----------- Read and clump exposure -----------
    exposure <- get.data(get.filename(trait=exp[1], subset=exp[2]), "exposure")
    exposure$exposure <- paste(exp, collapse=".")
    if(!is.null(ivs)) clump.exposure <- exposure[exposure$SNP %in% ivs,]
    else clump.exposure <- local.clump(data=exposure)
    
    for (out in outcome.lst) {
      #----------- Get outcome -----------
      outcome <- get.data(get.filename(trait=out[1], subset=out[2]), "outcome")
      outcome$outcome <- paste(out, collapse=".")
      
      #----------- Harmonize data -----------
      data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
      
      #----------- Run TwoSampleMR for unfiltered data -----------
      dt <- rbind(dt,run.TwoSampleMR(data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse=".")), folder=folder), fill=TRUE)
      
      if(!is.null(dt.steiger)) {
        #----------- Apply steiger filtering on data -----------
        filtered.data <- steiger_filtering(data)
        filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
        
        #----------- Run TwoSampleMR for Steiger-filtered data -----------
        dt.steiger <- rbind(dt.steiger,run.TwoSampleMR(filtered.data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse="."), "_Steiger"), folder=folder), fill=TRUE)
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

#--------------- Read summary stats --------------
OA.phen <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "THR", "TJR")

# Read and clump T2D IVs
T2D <- fread("/project_data/data/T2D_rsid.csv")
T2D.ivs <- ld_clump(dplyr::tibble(rsid=T2D$rsID, pval=T2D$pval),
                    plink_bin=genetics.binaRies::get_plink_binary(),
                    bfile="/reference_data/1kG/EUR_1kg_v3/EUR", #"/project_data/data_original/1kg_v3/EUR",
                    clump_p=5e-8)$rsid
rm(T2D)

# # Read and clump T2D UKBB based on T2D IVs
# T2D.UKBB <- fread("/storage/hmgu/projects/OA_T2D/processed_data/T2D_UKBB_rsid.csv")
# T2D.UKBB$Ncases <- 19119
# T2D.UKBB$N <- 442817
# T2D.UKBB$Ncontrols <- 423698
# 
# # Read and clump T2D noUKBB based on T2D IVs
# T2D.noUKBB <- fread("/project_data/data/T2D.noUKBB.gz")
# T2D.noUKBB$Ncases <- 55005
# T2D.noUKBB$N <- 455313
# T2D.noUKBB$Ncontrols <- 400308
# 
# # Read OA.T2D
# OA.T2D <- fread("/project_data/GWAS_OA_T2D/regenie/ukb_assoc/out_info03/OA.T2D.assoc.txt.gz") %>% 
#   .[, .(rsID=ID,CHR=CHROM, POS=GENPOS, pval=PVAL, beta=BETA, se=SE, 
#         EAF=A1FREQ, N=31472, Ncases=6799, EA=ALLELE1, NEA=ALLELE0, 
#         Ncontrols=24673)]

##############################################################################
#----------------- Exposure: T2D.UKBB Outcome: GO.noUKBB --------------------#
##############################################################################
# dt <- data.table()
# dt.steiger <- data.table()
# 
# wrap.MR(exposure.lst=list(c("T2D", "UKBB")), 
#         outcome.lst=mapply(append, as.list(OA.phen), as.list(rep("noUKBB", length(OA.phen))), SIMPLIFY=FALSE),
#         ivs=T2D.ivs, 
#         folder="vanilla_extra",
#         dt=dt,
#         dt.steiger=dt.steiger,
#         output.file="T2D.UKBB_OA.noUKBB")

##############################################################################
#------------- Exposure: T2D.noUKBB Outcome: GO.UKBB, OA.T2D ----------------#
##############################################################################
# dt <- data.table()
# dt.steiger <- data.table()
# 
# wrap.MR(exposure.lst=list(c("T2D", "noUKBB")), 
#         outcome.lst=mapply(append, as.list(c(OA.phen, "OA.T2D")), as.list(rep("UKBB", length(OA.phen)+1)), SIMPLIFY=FALSE), 
#         ivs=T2D.ivs, 
#         folder="vanilla_extra",
#         dt=dt,
#         dt.steiger=dt.steiger,
#         output.file="T2D.noUKBB_OA.UKBB")

##############################################################################
#-------------- Exposure: GO.noUKBB Outcome: T2D.UKBB, OA.T2D ---------------#
##############################################################################
dt <- data.table()
dt.steiger <- data.table()

wrap.MR(exposure.lst=mapply(append, as.list(OA.phen), as.list(rep("noUKBB", length(OA.phen))), SIMPLIFY=FALSE),
        outcome.lst=list(c("T2D", "UKBB"), c("OA.T2D", "UKBB")),
        folder="vanilla_extra",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="OA.noUKBB_T2D.UKBB")

##############################################################################
#----------------- Exposure: GO.UKBB Outcome: T2D.noUKBB --------------------#
##############################################################################
dt <- data.table()
dt.steiger <- data.table()

# wrap.MR(exposure.lst=mapply(append, as.list(OA.phen), as.list(rep("UKBB", length(OA.phen))), SIMPLIFY=FALSE), 
wrap.MR(exposure.lst=mapply(append, as.list(OA.phen[4]), as.list(rep("UKBB", 1)), SIMPLIFY=FALSE), 
        outcome.lst=list(c("T2D", "noUKBB")), 
        folder="vanilla_extra",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="OA.UKBB_T2D.noUKBB")

# #----------- Get exposure and clumped exposure -----------
# exposure <- get.data(T2D.UKBB, "exposure")
# exposure$exposure <- "T2D.UKBB"
# clump.exposure <- exposure[exposure$SNP %in% T2D.ivs,]
# 
# for (phen in OA.phen) {
#   #----------- Read outcome summary stats -----------
#   outcome <- get.data(get.GO.filename(trait=phen, subset="noUKBB"), "outcome")
#   outcome$outcome <- phen
#   
#   #----------- Harmonize data -----------
#   data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
#   
#   #----------- Run TwoSampleMR for unfiltered data -----------
#   tmp.dt <- run.TwoSampleMR(data, filename=paste0(exp.phen, "_", phen, ".noUKBB"), folder="vanilla_extra")
#   dt <- rbind(dt,tmp.dt)
#   
#   #----------- Apply steiger filtering on data -----------
#   filtered.data <- steiger_filtering(data)
#   filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
#   
#   #----------- Run TwoSampleMR for Steiger-filtered data -----------
#   tmp.dt <- run.TwoSampleMR(filtered.data, filename=paste0(exp.phen, "_", phen, ".noUKBB_Steiger"), folder="vanilla_extra")
#   dt.steiger <- rbind(dt.steiger,tmp.dt)
# }
# 
# dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt <- output.wide.result(dt)
# fwrite(dt, paste0("/home/nalu/T2D_OA/causal_inference/vanilla_extra/T2D.UKBB_OA.noUKBB_analysis_MR.csv"))
# 
# ##############################################################################
# #-------------- Exposure: GO.noUKBB Outcome: T2D.UKBB, OA.T2D ---------------#
# ##############################################################################
# dt <- data.table()
# dt.steiger <- data.table()
# 
# wrap.MR(exposure.lst=mapply(append, as.list(OA.phen), as.list(rep("noUKBB", length(OA.phen))), SIMPLIFY=FALSE), 
#         outcome.lst=list(c("T2D", "UKBB"), c("OA.T2D", "UKBB")), 
#         folder="vanilla_extra",
#         dt=dt,
#         dt.steiger=dt.steiger,
#         output.file="OA.noUKBB_T2D.UKBB")
# 
# for (phen in OA.phen) {
#   #----------- Read exposure summary stats -----------
#   exposure <- get.data(get.GO.filename(trait=phen, subset="noUKBB"), "exposure")
#   exposure$exposure <- phen
#   
#   #----------- Clump exposure summary stats -----------
#   clump.exposure <- local.clump(data=exposure)
#   
#   for (out.phen in c("OA.T2D", "T2D.UKBB")) {
#       #----------- Read outcome summary stats -----------
#       if (out.phen == "OA.T2D") outcome <- get.data(OA.T2D, "outcome")
#       else outcome <- get.data(T2D.UKBB, "outcome")
#       outcome$outcome <- out.phen
#       
#       #----------- Harmonize data -----------
#       data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
#       
#       #----------- Run TwoSampleMR for unfiltered data -----------
#       tmp.dt <- run.TwoSampleMR(data, filename=paste0(phen, ".noUKBB_", out.phen), folder="vanilla_extra")
#       dt <- rbind(dt,tmp.dt)
#       
#       #----------- Apply steiger filtering on data -----------
#       filtered.data <- steiger_filtering(data)
#       filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
#       
#       #----------- Run TwoSampleMR for Steiger-filtered data -----------
#       tmp.dt <- run.TwoSampleMR(filtered.data, filename=paste0(phen, ".noUKBB_", out.phen, "_Steiger"), folder="vanilla_extra")
#       dt.steiger <- rbind(dt.steiger,tmp.dt)
#   }
# }
# 
# dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# fwrite(dt, paste0("/home/nalu/T2D_OA/causal_inference/vanilla_extra/long_OA.noUKBB_T2D.UKBB_analysis_MR.csv"))
# dt2 <- output.wide.result(unique(dt))
# fwrite(dt2, paste0("/home/nalu/T2D_OA/causal_inference/vanilla_extra/OA.noUKBB_T2D.UKBB_analysis_MR.csv"))
# 
# ##############################################################################
# #------------- Exposure: T2D.noUKBB Outcome: GO.UKBB, OA.T2D ----------------#
# ##############################################################################
# dt <- data.table()
# dt.steiger <- data.table()
# 
# wrap.MR(exposure.lst=list(c("T2D", "noUKBB")), 
#         # outcome.lst=mapply(append, as.list(c(OA.phen, "OA.T2D")), as.list(rep("UKBB", length(OA.phen)+1)), SIMPLIFY=FALSE), 
#         outcome.lst=mapply(append, as.list(OA.phen[c(1,2)]), as.list(rep("UKBB", 2)), SIMPLIFY=FALSE),
#         ivs=T2D.ivs, 
#         folder="vanilla_extra",
#         dt=dt,
#         dt.steiger=dt.steiger,
#         output.file="T2D.noUKBB_OA.UKBB")
# 
# #----------- Get exposure and clumped expoure -----------
# exposure <- format_data(type = "exposure",
#                         dat = T2D.noUKBB,
#                         snp_col = "SNP",
#                         beta_col = "BETA",
#                         se_col = "StdErr",
#                         effect_allele_col = "EA",
#                         other_allele_col = "NEA",
#                         eaf_col = "Freq1",
#                         pval_col = "P",
#                         chr_col="CHR",
#                         pos_col="POS",
#                         ncase_col="Ncases",
#                         samplesize_col="N",
#                         ncontrol_col="Ncontrols")
# 
# 
# exposure$exposure <- "T2D.noUKBB"
# clump.exposure <- exposure[exposure$SNP %in% T2D.ivs,]
# 
# for (phen in c(OA.phen, "OA.T2D")) {
#   #----------- Read outcome summary stats -----------
#   if (phen == "OA.T2D") outcome <- get.data(OA.T2D, "outcome")
#   else outcome <- get.data(get.GO.filename(trait=phen, subset="UKBB"), "outcome")
#   outcome$outcome <- phen
#   
#   #----------- Harmonize data -----------
#   data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
#   
#   #----------- Run TwoSampleMR for unfiltered data ----------- <- <- <- <- 
#   tmp.dt <- run.TwoSampleMR(data, filename=paste0("T2D.noUKBB_", phen, ".UKBB"), folder="vanilla_extra")
#   dt <- rbind(dt,tmp.dt)
#   
#   #----------- Apply steiger filtering on data -----------
#   filtered.data <- steiger_filtering(data)
#   filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
# 
#   #----------- Run TwoSampleMR for Steiger-filtered data -----------
#   tmp.dt <- run.TwoSampleMR(filtered.data, filename=paste0("T2D.noUKBB_", phen, ".UKBB_Steiger"), folder="vanilla_extra")
#   dt.steiger <- rbind(dt.steiger,tmp.dt)
# }
# 
# #--------------- Unfiltered data ----------------
# dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt <- output.wide.result(dt)
# fwrite(dt, paste0("/project_data/overlap_T2D_OA/causal_inference/vanilla_extra/T2D.noUKBB_OA.UKBB_analysis_MR.csv"))
# 
# #--------------- Steiger-filtered data ----------------
# dt.steiger[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt.steiger <- output.wide.result(dt.steiger)
# fwrite(dt.steiger, paste0("/project_data/overlap_T2D_OA/causal_inference/vanilla_extra/Steiger_T2D.noUKBB_OA.UKBB_analysis_MR.csv"))
# 
# ##############################################################################
# #----------------- Exposure: GO.UKBB Outcome: T2D.noUKBB --------------------#
# ##############################################################################
# dt <- data.table()
# dt.steiger <- data.table()
# 
# wrap.MR(exposure.lst=mapply(append, as.list(OA.phen), as.list(rep("UKBB", length(OA.phen))), SIMPLIFY=FALSE), 
#         outcome.lst=list(c("T2D", "noUKBB")), 
#         folder="vanilla_extra",
#         dt=dt,
#         dt.steiger=dt.steiger,
#         output.file="OA.UKBB_T2D.noUKBB")
# 
# for (phen in c(OA.phen)) {
#   #----------- Read exposure summary stats -----------
#   exposure <- get.data(get.GO.filename(trait=phen, subset="UKBB"), "exposure")
#   exposure$exposure <- phen
#   
#   #----------- Clump exposure summary stats -----------
#   clump.exposure <- local.clump(data=exposure)
#   
#   #----------- Read outcome summary stats -----------
#   outcome <- format_data(type = "outcome",
#                          dat = T2D.noUKBB,
#                          snp_col = "SNP",
#                          beta_col = "BETA",
#                          se_col = "StdErr",
#                          effect_allele_col = "EA",
#                          other_allele_col = "NEA",
#                          eaf_col = "Freq1",
#                          pval_col = "P",
#                          chr_col="CHR",
#                          pos_col="POS",
#                          ncase_col="Ncases",
#                          samplesize_col="N",
#                          ncontrol_col="Ncontrols")
#   outcome$outcome <- "T2D.noUKBB"
#   
#   #----------- Harmonize data -----------
#   data <- subset(harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
#   
#   #----------- Run TwoSampleMR for unfiltered data ----------- <- <- <- <- 
#   tmp.dt <- run.TwoSampleMR(data, filename=paste0("T2D.noUKBB_", phen, ".UKBB"), folder="vanilla_extra")
#   dt <- rbind(dt,tmp.dt)
#   
#   #----------- Apply steiger filtering on data -----------
#   filtered.data <- steiger_filtering(data)
#   filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
#   
#   #----------- Run TwoSampleMR for Steiger-filtered data -----------
#   tmp.dt <- run.TwoSampleMR(filtered.data, filename=paste0("T2D.noUKBB_", phen, ".UKBB_Steiger"), folder="vanilla_extra")
#   dt.steiger <- rbind(dt.steiger,tmp.dt)
# }
# 
# #--------------- Unfiltered data ----------------
# dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt <- output.wide.result(dt)
# fwrite(dt, paste0("/project_data/overlap_T2D_OA/causal_inference/vanilla_extra/OA.UKBB_T2D.noUKBB_analysis_MR.csv"))
# 
# #--------------- Steiger-filtered data ----------------
# dt.steiger[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
# dt.steiger <- output.wide.result(dt.steiger)
# fwrite(dt.steiger, paste0("/project_data/overlap_T2D_OA/causal_inference/vanilla_extra/Steiger_OA.UKBB_T2D.noUKBB_analysis_MR.csv"))




library(TwoSampleMR)
library(data.table)
library(dplyr)
library(hash)
library(biomaRt)
library(ieugwasr)
library(stringr)
library(ggplot2)
library(readxl)

################################################
#------------------- Functions -----------------
################################################
sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
}) 

local.clump <- function(data, pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure, id=data$id.exposure),
                         plink_bin=genetics.binaRies::get_plink_binary(),
                         bfile="/project_data/data_original/1kg_v3/EUR",
                         clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

get.outcome <- function(dt, ivs=NULL){
  outcome <- format_data(dt,
                         snps=ivs,
                         type="outcome",
                         snp_col = "rsID",
                         beta_col = "beta",
                         se_col = "se",
                         effect_allele_col = "EA",
                         other_allele_col = "NEA",
                         eaf_col = "EAF",
                         pval_col = "pval",
                         chr_col="CHR",
                         pos_col="POS",
                         gene_col="geneID")
  return(outcome)
}

output.wide.result <- function(dt){
  # dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]
  dt <- dcast(dt, exposure + outcome + nsnp ~ method,
              value.var=c("b", "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95", "Q", "Q_df", "Q_pval", "plt.egger_intercept", "plt.pval", "plt.se", "p.adj.fdr", "Fstat"))
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

read.fat <- function(trait){
  if(trait=="BMI") {
    dt <- fread("/project_data/data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>% 
      .[, .(CHR, POS, rsID=sub(":.*", "", SNP), pval=P, beta=BETA, se=SE, MAF=as.double(NA), N=as.integer(N), Ncases=as.integer(NA), EA=Tested_Allele, NEA=Other_Allele, ID=as.character(NA))]
    
  }
  else if(trait=="WHR"){
    dt <- fread("/project_data/data/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>% 
      .[, .(CHR, POS, rsID=sub(":.*", "", SNP), pval=P, beta=BETA, se=SE, MAF=as.double(NA), N=as.integer(N), Ncases=as.integer(NA), EA=Tested_Allele, NEA=Other_Allele, ID=as.character(NA))]
    
  }
  else if(trait=="body_mass"){
    dt <- fread("/project_data/data/whole_body_fat_mass_neale_rsid.txt.gz") %>% 
      .[, .(CHR, POS, rsID=SNP, pval, beta, se, MAF=as.double(NA), N=as.integer(n_complete_samples), Ncases=as.integer(NA), EA=minor_allele, NEA=major_allele, ID=as.character(NA))] %>%
      .[rsID!=""]
    
  }
  else if(trait=="fat_per"){
    dt <- fread("/project_data/data/body_fat_percentage_neale_rsid.txt.gz") %>% 
      .[, .(CHR, POS, rsID=SNP, pval, beta, se, MAF=as.double(NA), N=as.integer(n_complete_samples), Ncases=as.integer(NA), EA=minor_allele, NEA=major_allele, ID=as.character(NA))] %>%
      .[rsID!=""]
    
  }
  return(dt)
}

######################################################
#--------------- Add indep. eQTL column --------------
######################################################
# add column to eQTL datasets: indep.eQTL (=T/F)
# my.func <- function(snp){
#   clean.dt <- dt[eQTL==1]
#   clumped.dt <- ld_clump(dplyr::tibble(rsid=clean.dt$rsID, clean.dt$pval, clean.dt$ID),
#                          plink_bin=genetics.binaRies::get_plink_binary(),
#                          bfile="/reference_data/1kG/EUR_1kg_v3/EUR")
#   if(rsID %in% clumped.dt$rsid) return(TRUE)
#   else return(FALSE)
# }
# #dt[, indep.eQTL:=my.func(snp=rsID, dt=dt), by="geneID"]
# 
# for (t in tissues){
#   t <- tissues[1]
#   dt <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/eQTL_", t, "_complete.csv"))
#   dt[, indep.eQTL:=FALSE]
# 
#   # for each gene: check if there are eQTLs, if yes clump them
#   for (g in hc.genes$gene.id){
#     g = "ENSG00000002933"
#     clean.dt <- dt[geneID==g & eQTL==1]
#     if(nrow(clean.dt)>0){
#       clumped.dt <- ld_clump(dplyr::tibble(rsid=clean.dt$rsID, clean.dt$pval, clean.dt$ID),
#                              plink_bin=genetics.binaRies::get_plink_binary(),
#                              bfile="/reference_data/1kG/EUR_1kg_v3/EUR",
#                              clump_p=5e-8)
#       dt[rsID %in% clumped.dt$rsid & geneID==g, indep.eQTL:=TRUE]
#     }
#   }
#   # dt[rsID %in% clumped.dt$rsid, indep.eQTL:=TRUE]
#   #dt[eQTL==1, indep.eQTL:=local.clump(gene=geneID, snp=rsID), by=seq_len(nrow(dt))]
#   
#   dt[eQTL==1, .N]
#   dt[indep.eQTL==TRUE, .N]
#   fwrite(dt, paste0("/project_data/data/eQTL_", t, ".csv"))
# }
# 
# t <- "PancreaticIslets"
# dt <- fread(paste0("/storage/hmgu/projects/OA_T2D/processed_data/eQTL_", t, "_EAF.csv"))
# pancreatic.islets.ivs <- fread("/storage/hmgu/projects/OA_T2D/data_original/PancreaticIslets_independent_gene_eQTLs.txt") %>% .[, .(gene.name=GeneName, eQTL=SNPid, gene.id=sub("\\..*", "", GeneID))]
# 
# dt[, indep.eQTL:=FALSE]
# for (g in hc.genes$gene.id){
#    dt[rsID %in% pancreatic.islets.ivs[gene.id==g, eQTL] & geneID==g, indep.eQTL:=TRUE]
# }
# fwrite(dt, paste0("/project_data/data/eQTL_", t, ".csv"))

# dt[, indep.eQTL:=ifelse(rsID %in% pancreatic.islets.ivs[gene.id==geneID, eQTL], TRUE, FALSE), by=seq_len(nrow(dt))]
# dt[, indep.eQTL2:=ifelse(rsID %in% pancreatic.islets.ivs$eQTL, TRUE, FALSE), by="geneID"]
# dt[indep.eQTL==TRUE,.N]
# dt[indep.eQTL2==TRUE,.N]
# dt[indep.eQTL3==TRUE,.N]
# dt[eQTL==1, .N]

###################################################
#--------------- Read summary stats --------------
###################################################
hc.genes <- as.data.table(read_excel("/project_data/overlap_T2D_OA/genes_table.xlsx", 5, skip=1)[c(1,2)]) %>% setnames(c("gene.symbol", "gene.id"))
tissues <- c("Synovium", "HighGradeCartilage", "LowGradeCartilage", "PancreaticIslets")

##################################################################################################
#---- Step 1: exposure=BMI, outcome=gene expression, IVs=BMI ivs excluding independent eQTLs ----
##################################################################################################
#----------- Read BMI summary stats -----------
# bmi <- extract_instruments("ieu-b-40")
# bmi.ivs <- bmi$SNP
# bmi$exposure <- "bmi"
# 
# whr <- extract_instruments("ieu-a-72")
# whr.ivs <- whr$SNP
# whr$exposure <- "whr"

whole_body_fat_mass <- extract_instruments("ukb-a-265")
whole_body_fat_mass.ivs <- whole_body_fat_mass$SNP
whole_body_fat_mass$exposure <- "whole_body_fat_mass"

body_fat_per <- extract_instruments("ukb-a-264")
body_fat_per.ivs <- body_fat_per$SNP
body_fat_per$exposure <- "body_fat_per"

get.exposure <- function(dt, ivs=NULL){
  exposure <- format_data(dt,
                          snps=ivs,
                          type="exposure",
                          snp_col = "rsID",
                          beta_col = "beta.combined",
                          se_col = "se.combined",
                          effect_allele_col = "A1.combined",
                          other_allele_col = "A2.combined",
                          eaf_col = "frqA1.combined",
                          pval_col = "pval.combined",
                          chr_col="Chr.ref.males",
                          pos_col="Pos.ref.males")
  return(exposure)
}

# bmi <- get.exposure(fread("/project_data/data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>%
#                       .[, rsID:=sub(":.*", "", SNP)])
# bmi.ivs <- fread("/project_data/data/bmi.giant-ukbb.meta.1.merged.indexSnps.combined.parsed.txt") %>% 
#   .[,  rsID:=sub(":.*", "", SNP)] %>% .[, rsID]
bmi <- get.exposure(fread("/project_data/data/bmi.giant-ukbb.meta.1.merged.indexSnps.combined.parsed.txt") %>%
                      .[, rsID:=sub(":.*", "", SNP)])
bmi.ivs <- bmi$SNP
bmi$exposure <- "bmi"

whr <- get.exposure(fread("/project_data/data/whr.giant-ukbb.meta.1.merged.indexSnps.combined.parsed.txt") %>%
                      .[, rsID:=sub(":.*", "", SNP)])
whr.ivs <- whr$SNP
whr$exposure <- "whr"

# cur.tissue=tissues[1]
# exp=bmi
# cur.gene=hc.genes$gene.id[1]

dt <- data.table()
for (exp in c("bmi", "whr", "whole_body_fat_mass", "body_fat_per")){  # bmi
  exp.dt <- data.table()
  exp=get(exp)
  for (cur.tissue in tissues){
    tissues.dt <- fread(paste0("/project_data/data/eQTL_", cur.tissue, ".csv")) %>% .[rsID %in% exp$SNP & indep.eQTL==FALSE, ]
    for (cur.gene in hc.genes$gene.id){
      cur.gene.symbol <- hc.genes[gene.id==cur.gene, gene.symbol]
      tissue.dt <- tissues.dt[geneID==cur.gene]
      # data <- harmonise_data(exposure_dat=exp, outcome_dat=outcome)
      
      #----------- Read outcome summary stats -----------
      if (nrow(tissue.dt[geneID==cur.gene])<1) print(paste0("No ", unique(exp$exposure), " risk variant for gene ", cur.gene.symbol, " in ", cur.tissue))
      else{
        outcome <- get.outcome(tissue.dt[geneID==cur.gene])
        outcome$outcome <- paste(cur.tissue, cur.gene.symbol, sep="_")
        
        #----------- Harmonize data -----------
        if(nrow(exp[exp$SNP %in% outcome$SNP,])==0) print("No valid SNP in common between exposure and outcome, not able to run MR :(")
        else {
          data <- harmonise_data(exposure_dat=exp, outcome_dat=outcome)
        
          if(nrow(data)==0) print("No valid SNP in common between exposure and outcome, not able to run MR :(")
          else {
            data <- subset(data, mr_keep==TRUE)
            if(nrow(data)==0) print("No valid SNP in common between exposure and outcome, not able to run MR :(")
            else if(nrow(data)==1) {
              #----------- Run TwoSampleMR wald ratio method -----------
              print("Only one SNP left: running wald ratio method")
              
              snp <- data$SNP
              res <- mr(data, method_list="mr_wald_ratio")
              res <- generate_odds_ratios(res)
              
              fstat <- (res$b)^2/(res$se)^2
              tmp.dt <- as.data.table(res) %>% .[, `:=`(gene=cur.gene.symbol,tissue=cur.tissue,SNP=snp,Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
              exp.dt <- rbind(exp.dt,tmp.dt)
            }
            else{
              #----------- Run TwoSampleMR -----------
              snps <- data$SNP
              res <- mr(data, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
              res <- generate_odds_ratios(res)
              
              #------------ Sensitivity analysis
              fstat <- mean((res$b)^2/(res$se)^2)
              het <- mr_heterogeneity(data)   #heterogeneity statistics
              plt <- mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy 
              tmp.dt <- as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T)) %>% 
                .[, `:=`(gene=cur.gene.symbol,tissue=cur.tissue,SNP=paste(snps, collapse=", "),Fstat=fstat,plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))] %>%
                .[exposure==paste(cur.tissue, cur.gene.symbol, sep="_") & method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]
              exp.dt <- rbind(exp.dt,tmp.dt)
              
              #------------ Sensitivity plots
              #mr_scatter_plot(res, data)   #scatter plot
              #ggsave(paste0("/project_data/overlap_T2D_OA/causal_inference/geneexp/TwoSampleMR_", trait, "_", cur.tissue, "_", cur.gene, ".jpg"))
            }
          }
        }
      }
    }
  }
  
  dt <- rbind(dt, exp.dt)
  # dt$exposure=unique(exp$exposure)
  
  # ############## Calculate adjusted p-values
  exp.dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
  exp.dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]

  # fwrite(dt, "/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step1_", unique(exp$exposure), "_geneexp.csv")
  fwrite(exp.dt[pval<0.05], paste0("/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step1_", unique(exp$exposure), "_geneexp_significant.csv"))
  fwrite(output.wide.result(exp.dt), paste0("/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step1_", unique(exp$exposure), "_geneexp_wide.csv"))
}

############## Calculate adjusted p-values
dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]

# fwrite(dt, "/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step1_", unique(exp$exposure), "_geneexp.csv")
fwrite(dt[pval<0.05], paste0("/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step1_adiposity_geneexp_significant_19092022.csv"))
fwrite(output.wide.result(dt), "/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step1_adiposity_geneexp_wide_19092022.csv")

##################################################################################################
#------------ Step 2: exposure=gene expression, outcome=OA/T2D, IVS=independent eQTLs -----------
##################################################################################################
get.exposure <- function(dt, ivs=NULL){
  exposure <- format_data(dt,
                          snps=ivs,
                          type="exposure",
                          snp_col = "rsID",
                          beta_col = "beta",
                          se_col = "se",
                          effect_allele_col = "EA",
                          other_allele_col = "NEA",
                          eaf_col = "EAF",
                          pval_col = "pval",
                          chr_col="CHR",
                          pos_col="POS",
                          gene_col="geneID")
  return(exposure)
}

get.outcome <- function(trait, ivs=NULL){
  if(trait=="T2D"){
    outcome <- read_outcome_data(filename = "/project_data/data/T2D_rsid.csv",
                                 snps = ivs,
                                 sep = ",",
                                 snp_col = "rsID",
                                 beta_col = "beta",
                                 se_col = "se",
                                 effect_allele_col = "EA",
                                 other_allele_col = "NEA",
                                 eaf_col = "EAF",
                                 pval_col = "pval")
                                 #chr_col = "CHR",
                                 # pos_col = "POS")
  }
  else{
    outcome <- read_outcome_data(filename = paste0("/project_data/data/GO.FILTER.GW.", trait,".FULL.09052019.txt.gz"),
                                 snps = ivs,
                                 sep = "\t",
                                 snp_col = "SNP",
                                 #chr_col = "CHR",
                                 #pos_col = "POS",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "EA",
                                 other_allele_col = "NEA",
                                 eaf_col = "EAF",
                                 pval_col = "P",
                                 ncase_col = "NCASES",
                                 ncontrol_col = "NCONTROLS",
                                 samplesize_col = "N")
  }
  return(outcome)
}

tissues <- c("Synovium", "HighGradeCartilage", "LowGradeCartilage", "PancreaticIslets")
tissues.lst <- lapply(tissues[1:4], function(i) {fread(paste0("/project_data/data/eQTL_", i, ".csv")) %>% .[indep.eQTL==TRUE]})
names(tissues.lst) <- tissues

dt <- data.table()
for (trait in c("T2D", "AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "THR", "TJR")) {
  #----------- Read outcome summary stats -----------
  outcome <- get.outcome(trait)
  outcome$outcome <- trait 
  
  for (cur.tissue in names(tissues.lst)){
    for (cur.gene in hc.genes$gene.id){
      cur.gene.symbol <- hc.genes[gene.id==cur.gene, gene.symbol]
      tissue.dt <- tissues.lst[[cur.tissue]][geneID==cur.gene]

      #----------- Read exposure summary stats -----------
      if (nrow(tissue.dt)<1) print(paste0("No eQTL for gene ", cur.gene.symbol, " in ", cur.tissue))
      else{
        exposure <- get.exposure(tissue.dt)
        exposure$exposure <- paste(cur.tissue, cur.gene.symbol, sep="_")
        
        #------------ Harmonize data -----------                                                                                                                                                                                                        ---------- Harmonize data -----------
        data <- harmonise_data(exposure_dat=exposure, outcome_dat=outcome)
        
        if(nrow(data)==0) print("No valid SNP in common between exposure and outcome, not able to run MR :(")
        else {
          data <- subset(data, mr_keep==TRUE)
          if(nrow(data)==0) print("No valid SNP in common between exposure and outcome, not able to run MR :(")
          else if(nrow(data)==1) {
            #----------- Run TwoSampleMR wald ratio method -----------
            print("Only one SNP left: running wald ratio method")
            filtered.data <- steiger_filtering(data, get_r_from_lor=TRUE)
            filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
            if(is.na(filtered.data$SNP)) print("Steiger filtering data is NA")
            if(nrow(filtered.data)==0) print("Steiger filtering data has no rows")
            
            snp <- data$SNP
            res <- mr(data, method_list="mr_wald_ratio")
            res <- generate_odds_ratios(res)
            
            fstat <- (res$b)^2/(res$se)^2
            tmp.dt <- as.data.table(res) %>% .[, `:=`(gene=cur.gene.symbol,tissue=cur.tissue,SNP=snp,Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
            dt <- rbind(dt,tmp.dt)
          }
          else{
            #----------- Run TwoSampleMR -----------
            snps <- data$SNP
            res <- mr(data, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
            res <- generate_odds_ratios(res)
            
            #------------ Sensitivity analysis
            fstat <- mean((res$b)^2/(res$se)^2)
            het <- mr_heterogeneity(data)   #heterogeneity statistics
            plt <- mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy 
            tmp.dt <- as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T)) %>% 
              .[, `:=`(gene=cur.gene.symbol,tissue=cur.tissue,SNP=paste(snps, collapse=", "),Fstat=fstat,plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))] %>%
              .[exposure==paste(cur.tissue, cur.gene.symbol, sep="_") & method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]
            dt <- rbind(dt,tmp.dt)
            
            #------------ Sensitivity plots
            # mr_scatter_plot(res, data)   #scatter plot
            # ggsave(paste0("/home/nalu/T2D_OA/causal_inference/geneexp/TwoSampleMR_", cur.tissue, "_", cur.gene, "_", trait, ".jpg"))
            # mr_forest_plot(mr_leaveoneout(data))  #leave-one-out analysis
            # ggsave(paste0("/home/nalu/T2D_OA/causal_inference/geneexp/LeaveOneOut_", cur.tissue, "_", cur.gene, "_", trait, "_T2D.jpg"))
            # mr_forest_plot(mr_singlesnp(data))
            # ggsave(paste0("/home/nalu/T2D_OA/causal_inference/geneexp/SingleSNP_", cur.tissue, "_", cur.gene, "_", trait, "_T2D.jpg"))
          }
        }
      }
    }
  }
}

############## Calculate adjusted p-values
dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]

fwrite(dt, paste0("/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step2_geneexp_trait.csv"))
fwrite(dt[pval<0.05], paste0("/project_data/overlap_T2D_OA/causal_inference/2step_analysis/step2_geneexp_trait_significant.csv"))

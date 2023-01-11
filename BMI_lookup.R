# https://academic.oup.com/hmg/article/28/1/166/5098227?login=false#165124482
# https://zenodo.org/record/1251813

library(data.table)
library(ieugwasr)
library(R.utils)
library(readr)
library(stringr)

get.neale.gwas <- function(filename, gwas.name){
  dt <- gunzip(paste0("/project_data/data/", filename, ".tsv.bgz"), paste0("/project_data/data/", filename, ".tsv"))
  dt <- as.data.table(read_tsv(paste0("/project_data/data/", filename, ".tsv")))
  setDT(dt)[, c("CHR", "POS", "major_allele") := tstrsplit(variant, ":")[c(1,2,3)]]
  dt <- dt[nchar(major_allele)<2 & CHR!="X"]  # remove indels and chromosome X
  dt <- dt[, `:=`(CHR=as.integer(CHR), POS=as.integer(POS))]
  fwrite(dt, paste0("/project_data/data/", gwas.name, "_neale.txt.gz"))
}

#########################################################
#--------------------- Read sumstats --------------------
#########################################################
t2d.dt <- fread("/project_data/data/T2D_rsid.csv") %>% 
  .[, .(CHR, POS, rsID, P=pval, BETA=beta, SE=se, EA, NEA)]

#------------- BMI -------------
# BMI <- fread("/project_data/data/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
BMI <- fread("/project_data/data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>% .[, SNP:=sub(":.*", "", SNP)]
# BMI.neale <- fread("/project_data/data/BMI_neale.txt.gz")

#------------- WHR ---------------
WHR <- fread("/project_data/data/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>% .[, SNP:=sub(":.*", "", SNP)]

#------------- Whole body fat mass ------------------
# get.neale.gwas("23100_irnt.gwas.imputed_v3.both_sexes", "whole_body_fat_mass")
fat_mass <- fread("/project_data/data/whole_body_fat_mass_neale.txt.gz")

#------------- Body fat percentage ------------------
#get.neale.gwas("23099_irnt.gwas.imputed_v3.both_sexes", "body_fat_per")
body_fat <- fread("/project_data/data/body_fat_per_neale.txt.gz")

#########################################################
#--------------------- Run lookup --------------------
#########################################################
OA.phen <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR")

for (phen in OA.phen) {
  phen.dt <- fread(paste0("/project_data/data/GO.FILTER.GW.", phen, ".FULL.09052019.txt.gz")) %>%
    .[,.(CHR, POS, rsID=SNP, EA, NEA, BETA, SE, P)]
  
  credset <- fread(paste0("/project_data/processed_data/GWAS", phen, "/credible_set.csv"))
  dt <- merge(credset, BMI, all.x=TRUE, by.x=c("CHR", "POS", "rsID"), by.y=c("CHR", "POS", "SNP")) %>%
    .[, .(CHR, POS, rsID, region, SNP, PP4, SNP.PP4, BMI.EA=Tested_Allele, BMI.NEA=Other_Allele, BMI.beta=BETA, BMI.se=SE, BMI.p=P)]
  
  dt <- merge(dt, WHR, all.x=TRUE, by.x=c("CHR", "POS", "rsID"), by.y=c("CHR", "POS", "SNP")) %>%
    .[, .(CHR, POS, rsID, region, SNP, PP4, SNP.PP4, BMI.EA, BMI.NEA, BMI.beta, BMI.se, BMI.p, WHR.EA=Tested_Allele, WHR.NEA=Other_Allele, WHR.beta=BETA, WHR.se=SE, WHR.p=P)]
  
  dt <- merge(dt, fat_mass[, .(CHR, POS, fat_mass.EA=minor_allele, fat_mass.NEA=major_allele, fat_mass.beta=beta, fat_mass.se=se, fat_mass.p=pval)], all.x=TRUE, by=c("CHR", "POS"))
  dt <- merge(dt, body_fat[, .(CHR, POS, body_fat.EA=minor_allele, body_fat.NEA=major_allele, body_fat.beta=beta, body_fat.se=se, body_fat.p=pval)], all.x=TRUE, by=c("CHR", "POS"))
  
  res.dt <- merge(dt, t2d.dt, by=c("CHR", "POS"), all.x=TRUE)
  res.dt <- merge(res.dt, phen.dt, by=c("CHR", "POS"), all.x=TRUE, suffixes=c(".t2d", ".oa"))
  
  # Flip beta to have concordant EA
  res.dt <- res.dt[, BETA.oa:=ifelse(EA.oa!=EA.t2d, -BETA.oa, BETA.oa)]
  res.dt <- res.dt[, BMI.beta:=ifelse(BMI.EA!=EA.t2d, -BMI.beta, BMI.beta)]
  res.dt <- res.dt[, WHR.beta:=ifelse(WHR.EA!=EA.t2d, -WHR.beta, WHR.beta)]
  res.dt <- res.dt[, fat_mass.beta:=ifelse(fat_mass.EA!=EA.t2d, -fat_mass.beta, fat_mass.beta)]
  res.dt <- res.dt[, body_fat.beta:=ifelse(body_fat.EA!=EA.t2d, -body_fat.beta, body_fat.beta)]
  
  clean.dt <- res.dt[, `:=` (rsID.x=NULL, rsID.y=NULL, BMI.NEA=NULL, BMI.EA=NULL, WHR.EA=NULL, fat_mass.NEA=NULL, fat_mass.EA=NULL, WHR.NEA=NULL, body_fat.NEA=NULL, body_fat.EA=NULL, SNP=NULL)]
  fwrite(res.dt, paste0("/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/credset_bmi_", phen, ".csv"))
  fwrite(clean.dt, paste0("/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/credset_bmi_", phen, "_clean.csv"))
}


#---------------- Table per region ----------------
dt <- fread(paste0("/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/credset_bmi_AllOA_clean.csv")) %>% 
  .[, `:=` (BMI.direction=ifelse(BMI.beta>0, "+", "-"), WHR.direction=ifelse(WHR.beta>0, "+", "-"), fat_mass.direction=ifelse(fat_mass.beta>0, "+", "-"), body_fat.direction=ifelse(body_fat.beta>0, "+", "-"), T2D.direction=ifelse(BETA.t2d>0, "+", "-"), OA.direction=ifelse(BETA.oa>0, "+", "-"))] %>%
  .[, .(CHR, POS, rsID, BMI.p, BMI.direction, WHR.p, WHR.direction, fat_mass.p, fat_mass.direction, body_fat.p, body_fat.direction, T2D.direction, OA.direction)]
for (phen in OA.phen){
  tmp.dt <- fread(paste0("/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/credset_bmi_", phen, "_clean.csv")) %>% 
    .[, `:=` (BMI.direction=ifelse(BMI.beta>0, "+", "-"), WHR.direction=ifelse(WHR.beta>0, "+", "-"), fat_mass.direction=ifelse(fat_mass.beta>0, "+", "-"), body_fat.direction=ifelse(body_fat.beta>0, "+", "-"), T2D.direction=ifelse(BETA.t2d>0, "+", "-"), OA.direction=ifelse(BETA.oa>0, "+", "-"))] %>%
    .[, .(CHR, POS, rsID, BMI.p, BMI.direction, WHR.p, WHR.direction, fat_mass.p, fat_mass.direction, body_fat.p, body_fat.direction, T2D.direction, OA.direction)]
  # check if direction of effect is consistent across OA phenotypes
  dt <- merge(tmp.dt, dt, by=colnames(dt)[1:12], all=TRUE, suffixes=c("", paste0(".", phen)))
  dt <- dt[, OA.direction:=ifelse(!is.na(OA.direction), OA.direction, get(paste0("OA.direction.", phen))), by=seq_len(nrow(dt))]
  dt[, (paste0("OA.direction.", phen)):=NULL]
  # dt <- unique(rbind(dt, tmp.dt), by="rsID")
}

direction.check <- function(dt, colnames){
  tmp.dt <- dt[!is.na(get(colnames[1])) & !is.na(get(colnames[2]))]
  return(identical(tmp.dt[, get(colnames[1])], tmp.dt[, get(colnames[2])]))
}

for (phen1 in c("", ".KneeOA", ".KneeHipOA", ".HipOA", ".TKR", ".TJR", ".THR")){
  for (phen2 in c("", ".KneeOA", ".KneeHipOA", ".HipOA", ".TKR", ".TJR", ".THR")){
   if (direction.check(dt, c(paste0("OA.direction", phen1), paste0("OA.direction", phen2)))==FALSE){
     return(print(paste0(phen1, " & ", phen2)))
   }
  }
}


dt <- dt[order(CHR,POS)]
fwrite(dt, "/project_data/overlap_T2D_OA/gwas_eqtl_colocalization/regionwise_credset_bmi.csv")

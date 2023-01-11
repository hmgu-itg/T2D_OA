library(data.table)
library(GenomicSEM)
library(parallel)

# https://rdrr.io/github/MichelNivard/GenomicSEM/man/ldsc.html 
# https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS

permute <- function(dt.path, out.path){
  dt <- fread(dt.path)
  idx <- sample(1:nrow(dt[A1!=""]), replace=FALSE)
  # idx <- sample(1:nrow(dt), replace=FALSE)
  dt[A1!="", Z:=Z[idx]]
  # dt[, Z:=Z[idx]]
  fwrite(dt, out.path, sep="\t")
  # return(dt)
}

perm_corr <- function(t, s, n){
  # sample snpid column
  permute(paste0("/project_data/LDscore/", t, ".ldsc.sumstats.gz"), paste0("/project_data/LDscore/permuted.", t, ".", n, ".ldsc.sumstats.gz"))
  
  permute(paste0("/project_data/LDscore/T2D.ldsc.sumstats.gz"), paste0("/project_data/LDscore/permuted.T2D.", n, ".ldsc.sumstats.gz"))
  
  # permute("/project_data/LDscore/AllOA.ldsc.sumstats.gz","/project_data/LDscore/permuted.AllOA.1.ldsc.sumstats.gz")
  # permute("/project_data/LDscore/T2D.ldsc.sumstats.gz","/project_data/LDscore/permuted.T2D.1.ldsc.sumstats.gz")
  
  # run ldsc
  rg <- NA
  tryCatch({
    # ldsc <- GenomicSEM::ldsc(traits=c(paste0("/project_data/LDscore/permuted.T2D.", n, ".ldsc.sumstats.gz"), paste0("/project_data/LDscore/permuted.", t, ".", n, ".ldsc.sumstats.gz")), 
    ldsc <- GenomicSEM::ldsc(traits=c(paste0("/project_data/LDscore/T2D.ldsc.sumstats.gz"), paste0("/project_data/LDscore/permuted.", t, ".", n, ".ldsc.sumstats.gz")), 
                             sample.prev=c(0.089956, s),
                             population.prev=c(0.1, 0.0873),
                             ld="/project_data/LDscore/eur_w_ld_chr/",
                             wld="/project_data/LDscore/eur_w_ld_chr/",
                             trait.names=c("T2D", t),
                             ldsc.log="/project_data/LDscore/test_R",
                             stand=TRUE)
    rg <- as.numeric(ldsc$S_Stand[1,2])
  }, error=function(e){print("Error in ldsc")})
  
  # remove permutations and log files
  file.remove(paste0("/project_data/LDscore/permuted.", t, ".", n, ".ldsc.sumstats.gz"))
  file.remove(paste0("/project_data/LDscore/permuted.T2D.", n, ".ldsc.sumstats.gz"))
  
  # # remove log file
  file.remove("/project_data/LDscore/test_R_ldsc.log")
  
  return(rg)
}

corr <- function(t, s){
  ldsc <- GenomicSEM::ldsc(traits=c(paste0("/project_data/LDscore/T2D.ldsc.sumstats.gz"), paste0("/project_data/LDscore/", t, ".ldsc.sumstats.gz")), 
                           sample.prev=c(0.089956, s),
                           population.prev=c(0.1, 0.0873),
                           ld="/project_data/LDscore/eur_w_ld_chr/",
                           wld="/project_data/LDscore/eur_w_ld_chr/",
                           trait.names=c("T2D", "AllOA"),
                           ldsc.log="/project_data/LDscore/test_R",
                           stand=TRUE)
  # remove log file
  file.remove("/project_data/LDscore/test_R_ldsc.log")
  
  rg <- ldsc$S_Stand[1,2]
  return(rg)
}

traits.dt <- data.table(trait=c("T2D","AllOA", "KneeOA", "KneeHipOA", "HipOA", "TKR", "TJR", "THR"),
                        sample.prev=c(0.089956,0.27, 0.187, 0.224, 0.115, 0.078, 0.125, 0.078))

for (row in 2:nrow(traits.dt)){ 
  OA.trait <- traits.dt[row, trait]
  OA.sample.prev <- traits.dt[row, sample.prev]
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(33)
  out <- mclapply(1:20000, function(i){perm_corr(t=OA.trait, s=OA.sample.prev, n=i)}, mc.cores=30, mc.preschedule=TRUE)

  result.dt <- data.table(perm=numeric(), file1=character(), file2=character(), rg=numeric(), se=numeric())
  for (i in 1:length(out)){
    if (is.na(out[[i]])) result.dt <- rbind(result.dt, list(i, paste(t, i, sep="."), paste("T2D", i, sep="."), NA, NA))
    else result.dt <- rbind(result.dt, list(i, paste(t, i, sep="."), paste("T2D", i, sep="."), as.numeric(out[[i]]), NA))
  }
  result.dt
  rm(out)
  
  fwrite(result.dt, paste0("/project_data/LDscore/", OA.trait, ".corr.perm.R.csv"))
}


n1=seq(1,80)
n2=seq(1,250)
n=(n1-1)*250+n2


out <- mclapply(1:300, function(i){corr(t=OA.trait, s=OA.sample.prev)}, mc.cores=30, mc.preschedule=TRUE)
for (i in 1:length(out)){
  if (class(out[[i]])=="try-error") out[[i]] <- c(i, paste(t, i, sep="."), paste("T2D", i, sep="."), NA, NA)
  else out[[i]] <- c(i, paste(t, i, sep="."), paste("T2D", i, sep="."), as.numeric(out[[i]]), NA)
}
test <- as.data.table(do.call(rbind, out))
colnames(test) <- c("perm", "file1", "file2", "rg", "se")
test



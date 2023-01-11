# https://dgarcia-eu.github.io/SocialDataScience/5_SocialNetworkPhenomena/056_PermutationTests/PermutationTests

library(data.table)
library(ggplot2)

result <- data.table(trait=character(), n_tests=numeric(), empirical_p=numeric())

for (trait in c("AllOA", "KneeOA", "HipOA", "KneeHipOA", "TJR", "TKR", "THR")) { 
  obs.dt <- fread(paste0("/project_data/LDscore/out/", trait, ".correlation.out"))
  
  stats.dt <- fread(paste0("/project_data/LDscore/out/", trait, "2.correlation.perm.out"), fill=TRUE)
  stats.dt[!is.na(rg), .N]
  
  stats.dt <- rbind(stats.dt, fread(paste0("/project_data/LDscore/out/", trait, ".correlation.perm.out"), fill=TRUE))
  stats.dt[!is.na(rg), .N]

  stats.dt <- rbind(stats.dt, fread(paste0("/project_data/LDscore/out/", trait, "3.correlation.perm.out"), fill=TRUE))
  stats.dt[!is.na(rg), .N]
  
  # Calculate empirical p-value
  p <- (sum(stats.dt[!is.na(rg), rg]>=obs.dt$rg)+1)/stats.dt[!is.na(rg),.N]
  result <- rbind(result, as.list(c(trait, stats.dt[!is.na(rg), .N], p)))
  
  # Plot histogram + true value
  ggplot(stats.dt[!is.na(rg) & rg<=1 & rg>=-1], aes(rg)) +
    geom_histogram() +
    geom_vline(aes(xintercept=obs.dt$rg), color="red", size=0.5) +
    # ggtitle(paste0("Permutation test for genetic correlation between ", trait, " and T2D")) +
    labs(title=trait,
         x="genetic correlation",
         y="") +
    annotate("text", x = -0.25, y = 1000, label=paste0("P-value = ", round(p,4)), size=5) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(paste0("/project_data/LDscore/", trait, ".T2D.perm.correlation.png"), width=6, height=6, dpi=320)
}





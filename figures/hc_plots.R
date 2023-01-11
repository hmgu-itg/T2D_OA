library(dplyr)
library(httr)
library(data.table)
library(openxlsx)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(tidyr)
library(webshot)
library(networkD3)
library(eulerr)

########################### Functions ###########################
get.phen <- function(gene.name, gene.OA.phen){
  if (grepl("KneeOA|TKR", gene.OA.phen)){
    if (grepl("\\bHipOA|THR", gene.OA.phen)) gene.phen <- "both"
    else gene.phen <- "knee"
  }
  else if (grepl("\\bHipOA|THR", gene.OA.phen)) gene.phen <- "hip"
  else gene.phen <- "both"
  return(gene.phen)
}

my.venn <- function(data.lst, categories, filename){
  # https://venn.bio-spring.info/
  if (length(data.lst)==2) my.color=c("green", "magenta")
  if (length(data.lst)==3) my.color=c("red", "green", "blue")
  if (length(data.lst)==4) my.color=c("red", "green", "blue", "yellow")
  
  venn.diagram(
    x = data.lst,
    category.names = categories,
    col = "transparent",
    fill = my.color,
    alpha = 0.30,
    print.mode=c("raw","percent"),
    filename = paste0("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/", filename, ".png"),
    imagetype="png",
    # cex=2,
    # cat.cex = 2,
    # cat.pos = c(-30, 30),
    # cat.dist = c(0.035, 0.035),
    cat.default.pos = "outer",
    output=TRUE)
}

my.sankey <- function(data.dt, filename){
  # input: genes.dt[, .(source=..., target=..., group=...)]
  
  links <- data.dt %>% .[, `:=`(value=rep(1,nrow(data.dt)), group=as.factor(group), linkgroup=as.factor(linkgroup))]
  nodes <- rbind(unique(links[, .(names=source, group=as.factor(linkgroup))]), unique(links[, .(names=target, group)])) %>% as.data.frame()
  links <- as.data.frame(links)
  
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  my.color <- 'd3.scaleOrdinal() .range(["mediumseagreen", "lightskyblue", "yellowgreen"])'
  p <- sankeyNetwork(Links=links, Nodes=nodes, Source="IDsource", Target="IDtarget", Value="value", NodeID="names", 
                     fontSize=28, colourScale=my.color, NodeGroup="group", LinkGroup="group")
  saveNetwork(p, paste0("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/", filename, ".html"))
  webshot(paste0("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/", filename, ".html"),
          paste0("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/", filename, ".png"), vwidth = 1000, vheight = 900)
}

# --------------- Get high scoring genes ---------------------
potential.genes <- as.data.table(read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/genes_table.xlsx", sheet="potential_eg", startRow=2))
colnames(potential.genes) <- c("Gene", "Ensembl.ID", "Missense", "eGene.OA", "eGene.T2D",
                               "pGene.OA", "pGene.T2D", "DEG.OA", "DEG.T2D", "mouse.OA", 
                               "mouse.T2D", "OMIM.OA", "OMIM.T2D", "score", "HC.OA", 
                               "HC.T2D", "HC", "score2")
likely.genes <- potential.genes[score>=1]
hc.genes <- potential.genes[score>=3 | HC>0]

fwrite(likely.genes, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/likely_effector_genes.csv")
fwrite(hc.genes, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/high_confidence_genes.csv")

# --------------- Plot scoring distribution ---------------------
ggplot(data=potential.genes, aes(x=as.factor(score), fill=score)) +
  geom_bar() + 
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.6) +
  labs(x="Score", y="Number of genes", title="Potential effector genes' scoring")
ggsave("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/scoring_dist_potential.png")

ggplot(likely.genes, aes(x=as.factor(score), fill=score)) +
  geom_bar() + 
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.6) +
  labs(x="Score", y="Number of genes", title="Likely effector genes' scoring (score>=1)")
ggsave("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/scoring_dist_likely.png")

ggplot(hc.genes, aes(x=as.factor(score), fill=score)) +
  geom_bar() + 
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.6) +
  labs(x="Score", y="Number of genes", title="High confidence effector genes' scoring (score>=3)")
ggsave("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/scoring_dist_hc.png")

# --------------- Plot novel genes (not hc for T2D nor OA) ---------------------
ggplot(data=potential.genes, aes(fill=as.factor(HC), x=score)) + 
  geom_bar(position="stack") + 
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.6) +
  labs(x="Score", y="Number of genes", title="Potential effector genes' scoring")
ggsave("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/HC_scoring_dist_potential.png")

ggplot(likely.genes, aes(fill=as.factor(HC), x=score)) + 
  geom_bar(position="stack") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5) +
  labs(x="Score", y="Number of genes", title="Likely effector genes' scoring (score>=1)")
ggsave("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/HC_scoring_dist_likely.png")

ggplot(hc.genes, aes(fill=as.factor(HC), x=score)) + 
  geom_bar(position="stack") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.6) +
  labs(x="Score", y="Number of genes", title="High confidece effector genes' scoring (score>=3 | HC>0)")
ggsave("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/HC_scoring_dist_hc.png")

########################### Per phenotype ###########################
genes.phen.dt <- fread("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/list_genes_phen.csv")
likely.genes.phen.dt <- genes.phen.dt[name %in% likely.genes$Gene]
hc.genes.phen.dt <- genes.phen.dt[name %in% hc.genes$Gene]

likely.genes.phen.dt <- likely.genes.phen.dt[, phen:=get.phen(name, OA.phen), by=seq_len(nrow(likely.genes.phen.dt))]
hc.genes.phen.dt <- hc.genes.phen.dt[, phen:=get.phen(name, OA.phen), by=seq_len(nrow(hc.genes.phen.dt))]

#-------------------- Output list of genes -----------------
write(potential.genes$Gene, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/potential.txt")

write(likely.genes.phen.dt$name, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/likely.txt")
write(likely.genes.phen.dt[phen=="knee" | phen=="both", name], "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/likely_knee.txt")
write(likely.genes.phen.dt[phen=="hip" | phen=="both", name], "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/likely_hip.txt")

write(hc.genes.phen.dt$name, "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/hc.txt")
write(hc.genes.phen.dt[phen=="knee" | phen=="both", name], "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/hc_knee.txt")
write(hc.genes.phen.dt[phen=="hip" | phen=="both", name], "C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/pathway_analysis/hc_hip.txt")

#---------------- Venn plots for regions -----------------
knee <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
hip <- c(9,10,11,12,13,14,15,16,17,18)
my.venn(data.lst=list(knee, hip),
        categories=c("Knee OA", "Hip OA"),
        filename="region_kneeXhip")

# EULERR
# https://cran.r-project.org/web/packages/eulerr/vignettes/introduction.html
plot(euler(list(knee=knee, hip=hip)), quantities=TRUE)

#---------------- Venn plots for hc genes -----------------
my.venn(data.lst=list(hc.genes.phen.dt[phen=="knee" | phen=="both"]$name, hc.genes.phen.dt[phen=="hip" | phen=="both"]$name),
        categories=c("Knee OA", "Hip OA"),
        filename="hc_kneeXhip")
# EULERR
plot(euler(list(knee=hc.genes.phen.dt[phen=="knee" | phen=="both"]$name, hip=hc.genes.phen.dt[phen=="hip" | phen=="both"]$name)), quantities=TRUE)

#---------------- Venn plots for likely genes -----------------
my.venn(data.lst=list(likely.genes.phen.dt[phen=="knee" | phen=="both"]$id,likely.genes.phen.dt[phen=="hip" | phen=="both"]$id),
        categories=c("Knee OA", "Hip OA"),
        filename="likely_kneeXhip")
# EULERR
plot(euler(list(knee=knee, hip=hip)), quantities=TRUE)

#---------------- Sankey plots for hc genes -----------------
my.sankey(data.dt=as.data.table(unnest(hc.genes.phen.dt[, .(source=name, target=phen, linkgroup=phen)], cols=target)) 
          %>% .[, `:=` (target=ifelse(target=="both", "knee&hip", target))] 
          %>% .[,group:=target],
          filename="hc_sankey_OA")

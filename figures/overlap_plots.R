library(data.table)
library(openxlsx)
library(VennDiagram)
library(webshot)
library(networkD3)
library(dplyr)
library(tidyr)
library(eulerr)

setwd("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization")

####################### FUNCTIONS ########################
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
    filename = paste0("~/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/", filename, ".png"),
    imagetype="png",
    output=TRUE)
}

my.sankey <- function(data.dt, filename){
  # input: genes.dt[, .(source=..., target=..., group=...)]
  
  links <- data.dt %>% .[, `:=`(value=rep(1,nrow(data.dt)), group=as.factor(group), linkgroup=as.factor(linkgroup))]
  nodes <- rbind(unique(links[, .(names=source, group=as.factor(linkgroup))]), unique(links[, .(names=target, group)])) %>% as.data.frame()
  links <- as.data.frame(links)
  
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # if(length(unique(links$target))==4)  my.color <- 'd3.scaleOrdinal() .range(["paleturquoise", "lightsalmon", "lightgray", "yellowgreen"])'
  # if(length(unique(links$target))==3)  my.color <- 'd3.scaleOrdinal() .range(["paleturquoise", "lightsalmon", "yellowgreen"])' #"palegreen"
  if(length(unique(nodes$group))==4)  my.color <- 'd3.scaleOrdinal() .range(["lightsalmon", "mediumseagreen", "lightgray", "lightskyblue"])'
  if(length(unique(nodes$group))==3)  my.color <- 'd3.scaleOrdinal() .range(["mediumseagreen", "lightskyblue", "yellowgreen"])'
  
  p <- sankeyNetwork(Links=links, Nodes=nodes, Source="IDsource", Target="IDtarget", Value="value", NodeID="names", 
                     fontSize=28, colourScale=my.color, NodeGroup="group", LinkGroup="group")
  saveNetwork(p, paste0("comparison_plots/", filename, ".html"))
  webshot(paste0("comparison_plots/", filename, ".html"),
          paste0("comparison_plots/", filename, ".png"), vwidth = 1000, vheight = 900)
}

########################### DATA #############################
genes <- c("FAM150B", "JADE2", "JAZF1", "HIBADH", "TMEM176A", "FASTK", "MTMR9", "MSRA", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", "AF131215.8",
           "RP11-297N6.4", "SIGMAR1", "TCF7L2", "TSKU", "WSCD2", "UNG", "DIABLO", "ARL6IP4", "RP11-197N18.2", "TCTN2", "RILPL2", "CDK2AP1", "SERPINA1", "IRX3", "APOE", "ZNF226", "SNRPD2")
hc.genes <- c("JADE2", "JAZF1", "HIBADH", "TCF7L2", "WSCD2", "DIABLO", "IRX3", "APOE")
hc.genes2 <- c("JADE2", "TCF7L2", "WSCD2", "DIABLO", "IRX3", "APOE")

knee <- c("FAM150B", "JADE2", "TMEM176A", "FASTK", "TCF7L2", "WSCD2", "UNG", "DIABLO", "ARL6IP4", "RP11-197N18.2", "TCTN2", "RILPL2", "CDK2AP1", "SIGMAR1", "IRX3")
hip <- c("JAZF1", "HIBADH", "IRX3")
kneehip <- c("JAZF1", "HIBADH", "TMEM176A", "MTMR9", "MSRA", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", "AF131215.8",
              "RP11-297N6.4", "SIGMAR1", "TSKU", "WSCD2", "UNG", "DIABLO", "ARL6IP4", "RP11-197N18.2", "TCTN2", "RILPL2", "CDK2AP1", 
              "SERPINA1", "IRX3", "APOE", "ZNF226", "SNRPD2")
both <- c("MTMR9", "MSRA", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", "AF131215.8",
          "RP11-297N6.4", "TSKU", "SERPINA1", "IRX3", "APOE", "ZNF226", "SNRPD2")
# knee <- c(knee, both)
# hip <- c(hip, both)

HG <- c("HIBADH", "MTMR9", "TCF7L2", "WSCD2", "UNG", "DIABLO")
LG <- c("RP11-568K15.1","MTMR9", "SIGMAR1", "TCF7L2", "ZNF226")
Syn <- c("JAZF1","TMEM176A", "FASTK", "MTMR9", "MSRA", "TSKU", "SNRPD2")
PI <- c("FAM150B", "JADE2", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", 
        "AF131215.8", "RP11-297N6.4", "TCF7L2", "ARL6IP4", "RP11-197N18.2", "TCTN2", "RILPL2", "CDK2AP1", "IRX3")
missense <- c("APOE", "SERPINA1", "WSCD2")
cart <- c(HG, LG)

OAphen <- vector(mode="list", length=length(genes))
tissue <- vector(mode="list", length=length(genes))
disease <- vector(mode="character", length=length(genes))

for (i in 1:length(genes)){
  if (genes[i] %in% knee) OAphen[[i]] <- c(OAphen[[i]], "Knee")
  if (genes[i] %in% hip) OAphen[[i]] <- c(OAphen[[i]], "Hip")
  if (genes[i] %in% kneehip)OAphen[[i]] <- c(OAphen[[i]], "Both")
  
  if (genes[i] %in% cart) tissue[[i]] <- c(tissue[[i]], "Cartilage")
  if (genes[i] %in% PI) tissue[[i]] <- c(tissue[[i]], "Pancreatic islets")
  if (genes[i] %in% Syn) tissue[[i]] <- c(tissue[[i]], "Synovium")
  if (genes[i] %in% missense) tissue[[i]] <- c(tissue[[i]], "Missense")
  
  disease[i] <- "Missense"
  if (genes[i] %in% cart | genes[i] %in% Syn) disease[i] <- "OA"
  if(genes[i] %in% PI) {
    if (disease[i]=="Missense") disease[i] <- "T2D"
    else disease[i] <- "T2D&OA"
  }
}

genes.dt <- data.table(genes=genes, OAphen=OAphen, tissue=tissue, disease=disease) %>%
    .[, phen:=ifelse(grepl("Knee", OAphen) & !grepl("Hip", OAphen), "Knee", ifelse(grepl("Hip", OAphen) & !grepl("Knee", OAphen), "Hip", "Knee&Hip"))]
hc.genes.dt <- genes.dt[genes %in% hc.genes]
# genes.dt <- fread("~/OA_T2D/GWAS_eQTL_colocalization/table_genes.csv")[, tissue:=(list(tissue))]

########################### VENN PLOTS #############################
# --------------------- Plot per phenotype overlap ----------------------
my.venn(data.lst=list(genes.dt[grepl("Knee", phen)]$genes,genes.dt[grepl("Hip", phen)]$genes),
        categories=c("Knee OA", "Hip OA"),
        filename="kneeXhip")

# --------------------- Plot per tissue overlap ----------------------
my.venn(data.lst=list(genes.dt[grepl("Pancreatic islets", tissue)]$genes, 
                      genes.dt[grepl("Cartilage", tissue)]$genes, 
                      genes.dt[grepl("Synovium", tissue)]$genes),
        categories = c("Pancreatic islets", "Cartilage", "Synovium"),
        filename="overlap_tissue")

# --------------------- Plot combination of phenotype and tissue overlap ----------------------
my.venn(data.lst=list(genes.dt[grepl("Knee", OAphen) & grepl("Pancreatic islets", tissue)]$genes, 
                       genes.dt[grepl("Knee", OAphen) & grepl("Cartilage", tissue)]$genes, 
                       genes.dt[grepl("Knee", OAphen) & grepl("Synovium", tissue)]$genes),
         categories=c("Pancreatic islets", "Cartilage", "Synovium"),
         filenam="overlap_knee")


my.venn(data.lst=list(genes.dt[grepl("Hip", OAphen) & grepl("Pancreatic islets", tissue)]$genes, 
                       genes.dt[grepl("Hip", OAphen) & grepl("Cartilage", tissue)]$genes, 
                       genes.dt[grepl("Hip", OAphen) & grepl("Synovium", tissue)]$genes),
         categories=c("Pancreatic islets", "Cartilage", "Synovium"),
         filenam="overlap_hip")

my.venn(data.lst=list(genes.dt[grepl("Both", OAphen) & grepl("Pancreatic islets", tissue)]$genes, 
                       genes.dt[grepl("Both", OAphen) & grepl("Cartilage", tissue)]$genes, 
                       genes.dt[grepl("Both", OAphen) & grepl("Synovium", tissue)]$genes),
         categories=c("Pancreatic islets", "Cartilage", "Synovium"),
         filenam="overlap_kneehip")

# --------------------- Plot for high confidence genes ----------------------
my.venn(data.lst=list(hc.genes.dt[grepl("Knee", phen)]$genes, hc.genes.dt[grepl("Hip", phen)]$genes),
        categories=c("Knee OA", "Hip OA"),
        filename="hc_kneeXhip")

both <- c("reg1", "reg2", "reg9", "reg10","reg13", "reg16", "reg18", "reg19")
regions_knee <-c(both, "reg4", "reg5", "reg6", "reg8", "reg11", "reg12", "reg14", "reg15") 
regions_hip <- c(both, "reg3", "reg7")
my.venn(data.lst=list(regions_knee, regions_hip),
        categories=c("Knee OA", "Hip OA"),
        filename="region_kneeXhip")

# -------------------- EULERR
plot(euler(list(knee=hc.genes.dt[grepl("Knee", phen)]$genes, hip=hc.genes.dt[grepl("Hip", phen)]$genes)), quantities=TRUE)

########################### SANKEY PLOTS #############################
my.sankey(data.dt=as.data.table(unnest(hc.genes.dt[, .(source=genes, target=tissue, linkgroup=disease)], cols=target))
          %>% .[, group:= ifelse(target=="Pancreatic islets", "T2D", ifelse(target %in% c("Cartilage", "Synovium"), "OA", "Missense"))],
          filename="hc2_sankey_tissue")

my.sankey(data.dt=as.data.table(unnest(hc.genes.dt[, .(source=genes, target=OAphen, linkgroup=phen)], cols=target)) 
          %>% .[, `:=` (target=ifelse(target=="Both", "Knee&Hip", target))] 
          %>% .[,group:=target],
          filename="hc2_sankey_OA")



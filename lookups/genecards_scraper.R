library(rvest)
library(lubridate)
library(dplyr)
library(httr)
library(decapitated)
library(xml2)
library(jsonlite)
library(data.table)

################################################################################
# --------------------------------- FUNCTIONS ----------------------------------
################################################################################
groupPhenotype <- function(my.gene){
  gene.dt <- getEnsemblePhenotypes(my.gene)
  gene.dt$group <- lapply(gene.dt$description, getGroup)
  # gene.dt[, Ngroup:=.N, keyby=group]
  # gene.dt[, Ngroup := sapply(group, FUN = function(x) .N)]
  # gene.dt$Ngroup <- sapply(gene.dt$group, table)
  # gene.dt$Ngroup <- gene.dt$group
  phen.count <- table(unlist(gene.dt$group))
  print(phen.count)
  names(which.max(phen.count))
  return(gene.dt)
}

getEnsemblePhenotypes <- function(my.gene){
 # http://grch37.rest.ensembl.org/documentation/info/phenotype_gene 
  server <- "http://grch37.rest.ensembl.org"
  ext <- paste0("/phenotype/gene/homo_sapiens/", my.gene,"?include_associated=1")
  
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  
  stop_for_status(r)

  gene.dt <- data.table()
  gene.dt <- as.data.table(fromJSON(toJSON(content(r)))) %>% # .[!duplicated(description)]
    .[description!="ClinVar: phenotype not specified" & description!="Annotated by HGMD", .(description, Variation, source)] %>%
    .[, `:=`(description=as.character(description), Variation=as.character(Variation), source=as.character(source), group=rep("", nrow(gene.dt)))] 
  gene.dt[, Nvar:=.N, keyby = description]
  gene.dt <- gene.dt[!duplicated(description)]
  return(gene.dt)
}

getGroup <- function(phen){
  group=NULL
  if (grepl("body mass|bmi", phen, ignore.case=TRUE)){group=c(group, "BMI")}
  if (grepl("cancer|carcinoma|tumor|tumour|adenocarcinoma|lymphoma|sarcoma|melanoma|adenoma|myeloma|leukemia|medulloblastoma|astrocytoma|neuroblastoma|meningioma|glioma|nodule|neoplasm|MYELODYSPLASTIC", phen, ignore.case=TRUE)){group=c(group, "Cancer & tumor")}
  if (grepl("obesity|adiposity|adipose|fat|obese|lean mass", phen, ignore.case=TRUE)){group=c(group, "Adiposity")}
  if (grepl("weight|height|growth|size|circumference|waist|hip|anthropometric", phen, ignore.case=TRUE)){group=c(group, "Anthropometric")}
  if (grepl("bipolar|autism|neuroticism|extraversion|adventurousness|risk|mental health|neurociticism|depressive|bulimia", phen, ignore.case=TRUE)){group=c(group, "Psychological traits")}
  if (grepl("dementia|alzheimer|multiple sclerosis|schizophrenia|brain|subcortical|carpal tunnel|cerebrospinal|cognitive decline|TOURETTE|Parkinson's|antisaccade task|Amyotrophic lateral sclerosis", phen, ignore.case=TRUE)){group=c(group, "Neuronal traits")}
  if (grepl("sexual|sex|menarche", phen, ignore.case=TRUE)){group=c(group, "Sexual behaviour/trait")}
  if (grepl("asthma|lung|vital capacity|expiratory|FEV1|respiratory|pulmonary|adrenergics|inhalants", phen, ignore.case=TRUE)){group=c(group, "Respiratoty diseases & traits")}
  if (grepl("blood|hemoglobin|platelet|hematocrit|red cell|reticulocyte|reticulocytes|mean corpuscular volume|white cells|factor vii", phen, ignore.case=TRUE)){group=c(group, "Blood traits")}
  if (grepl("eosinophil|neutrophil|basophils|monocytes|white cells|immune|reactive protein|Lymphocyte|allergy|allergic|inflammatory|crohn's disease|inflammation|sKawasaki disease|HIV-1|Systemic sclerosis|hepatitis", phen, ignore.case=TRUE)){group=c(group, "Immune system & inflammation")}
  if (grepl("cholesterol|triglyceride|apolipoprotein|lipoproteins|HYPERLIPOPROTEINEMIA", phen, ignore.case=TRUE)){group=c(group, "Cholesterol traits")}  
  if (grepl("macular|ocular|retinal|retinitis|blind|visual|myopia|refractive|STARGARDT|electroretinogram|glaucoma|corneal|aniridia|Astigmatism|PI|peripheral iridotomy", phen, ignore.case=TRUE)){group=c(group, "Eye traits")}  
  if (grepl("diet|meat|vegetable|juice|fruit|beverage|fish|pork|dietary|lamb", phen, ignore.case=TRUE)){group=c(group, "Diet-related traits")}  
  if (grepl("smoke|smoking|alcohol|nicotine|coffee|coffeine", phen, ignore.case=TRUE)){group=c(group, "Drug-related traits")}  
  if (grepl("bone|osteo|arthritis|scoliosis|rhabdomyolysis|joint|polymyositis|brachydactyly", phen, ignore.case=TRUE)){group=c(group, "Musculoskeletal traits")}  
  if (grepl("lupus", phen, ignore.case=TRUE)){group=c(group, "Autoimmune diseases")}  
  if (grepl("cardiovascular|hypertension|heart|electrocardiogram|echocardiography|electrocardiography|pulse|carotid|stroke|ATHEROSCLEROSIS|Thrombosis|thromboembolism|calcium channel blockers|infarction|artery", phen, ignore.case=TRUE)){group=c(group, "Cardiovascular diseases")}  
  if (grepl("leptin|metabolic|Glycogen storage disease type X", phen, ignore.case=TRUE)){group=c(group, "Metabolism")}  
  if (grepl("Heschl's gyrus|ototoxicity|hearing|deafness", phen, ignore.case=TRUE)){group=c(group, "Auditaury system")}  
  if (grepl("lifespan|birth|longevity", phen, ignore.case=TRUE)){group=c(group, "Aging-related traits")}  
  if (grepl("dentures", phen, ignore.case=TRUE)){group=c(group, "Dental-related traits")} 
  if (grepl("kidney|urinary|glomerular filtration rate|Creatinine|diuretics|Dialysis", phen, ignore.case=TRUE)){group=c(group, "Renal diseases")}  
  if (grepl("hair|eyebrow|chin", phen, ignore.case=TRUE)){group=c(group, "Morphological traits")}  
  if (grepl("sleep|INSOMNIA|resting|chronotype", phen, ignore.case=TRUE)){group=c(group, "Sleep-related traits")}  
  if (grepl("thyroid|polycystic ovary syndrome", phen, ignore.case=TRUE)){group=c(group, "Hormonal diseases & traits")}  
  if (grepl("uric|urate|gout", phen, ignore.case=TRUE)){group=c(group, "?")}  
  if (grepl("glucose|magnesium|metabolite|lipids|iron|Cell Adhesion Molecules|potassium", phen, ignore.case=TRUE)){group=c(group, "?")}  
  if (grepl("beta-cell|diabetes|insulin|glucosuria", phen, ignore.case=TRUE)){group=c(group, "Diabetes")}   
  if (grepl("income|intellectual|acne|TCTN2-Related Disorders|exercise test|educational|intelligence|physical activity", phen, ignore.case=TRUE)){group=c(group, "Other")}  
  if (grepl("Joubert Syndrome|Meckel|Cutis laxa|Mucopolysaccharidosis-Plus Syndrome|fraxe", phen, ignore.case=TRUE)){group=c(group, "Rare genetic disorders")}  
  if (is.null(group)){group="None"}
  return(group)
}

getFunction <- function(my.gene){
  url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", my.gene,"&keywords=", my.gene)
  url2 <- paste0("https://www.ebi.ac.uk/gwas/genes/", my.gene)
  url3 <- paste0("http://grch37.ensembl.org/Homo_sapiens/Gene/Phenotype?db=core;g=", my.ensg,";r=", my.pos)
  
  x <- chrome_read_html(url) # chrome_read_html(url, render=FALSE)
  y <- chrome_read_html(url2)
  z <- chrome_read_html(url3)
  
  x2 <- html_node(x, '#function-phenotypes-from-gwas .ng-scope') %>% 
        html_table()
  
  y2 <- html_node(y, '#efotrait_panel .panel-body') %>% 
        html_table()
  
  z2 <- html_node(z, '#GenePhenotypeVariation') %>% 
        html_table()
}


################################################################################
# ----------------------------------- MAIN -------------------------------------
################################################################################
my.gene <- "FTO"
# work_dir <- "~/decapitated"

genes.lst <- c("RASAL2", "RP11-568K15.1", "TMEM18", "FAM150B", "JADE2", "SARB1", "JAZF1", "HIBADH", "ABP1", "NUB1", "TMEM176A", "FASTK", "MTMR9", "MSRA",
               "SOX7", "NEIL2", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", "AF131215.8",
               "RP11-297N6.4", "GLIS3", "RCL1", "TCF7L2", "RP11-57H14.2", "TSKU", "WSCD2", "DIABLO", "VPS33A", "RSRC2", "ARL6IP4", "RP11-197N18.2", "TCTN2", 
               "MPHOSPH9", "RILPL2", "CDK2AP1", "SERPINA1", "FTO", "IRX3", "APOE")
library(rvest)
library(lubridate)
library(dplyr)
library(httr)
library(decapitated)
library(xml2)
library(jsonlite)
library(data.table)
library(openxlsx)
library(ggplot2)
library(ggvenn)
library(ggVennDiagram)
library(VennDiagram)

################################################################################
# --------------------------------- FUNCTIONS ----------------------------------
################################################################################
getGroupedPhenotype <- function(my.gene){
  gene.dt <- getEnsemblePhenotypes(my.gene)
  
  if (length(gene.dt)!=0) {
    gene.dt$gene <- my.gene
    gene.dt$group <- lapply(gene.dt$description, getGroup)
    gene.dt$refGroup <- lapply(gene.dt$group, refineGroup)
    # phen.count <- table(unlist(gene.dt$group))
    phen.count <- table(unlist(gene.dt$refGroup))
    
    gene.summary <- list(gene=my.gene, topGroup=names(which.max(phen.count)), 
                         groupCount=max(phen.count), 
                         #score=max(phen.count)*gene.dt[group %like% names(which.max(phen.count)), sum(Nvar)])
                         score=max(phen.count)*gene.dt[refGroup %like% names(which.max(phen.count)), sum(Nvar)])
      
    return(list(gene.dt,  gene.summary))
  } else{
    print(paste0("No phenotype for gene ", my.gene, " available :("))
  }
}

getEnsemblePhenotypes <- function(my.gene){
 # http://grch37.rest.ensembl.org/documentation/info/phenotype_gene 
  server <- "http://grch37.rest.ensembl.org"
  ext <- paste0("/phenotype/gene/homo_sapiens/", my.gene,"?include_associated=1;non_specified=0;content-type=application/json")
  
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  
  stop_for_status(r)

  gene.dt <- data.table()
  gene.dt <- as.data.table(fromJSON(toJSON(content(r)))) # .[!duplicated(description)]
  gene.dt
  
  if (length(gene.dt)!=0) {
    gene.dt <- gene.dt[, .(description=as.character(description), Variation=as.character(Variation), 
                           source=as.character(source), group=rep("", nrow(gene.dt)))] %>%
      .[, Nvar:=.N, keyby = description] %>% .[!duplicated(description)]
  }
  
  return(gene.dt)
}

getGroup <- function(phen){
  group=NULL
  if (grepl("body mass|bmi", phen, ignore.case=TRUE)){group=c(group, "BMI")}
  if (grepl("cancer|carcinoma|tumor|tumour|adenocarcinoma|lymphoma|sarcoma|melanoma|adenoma|myeloma|leukemia|medulloblastoma|astrocytoma|neuropathy|neuroblastoma|meningioma|glioma|nodule|neoplasm|MYELODYSPLASTIC", phen, ignore.case=TRUE)){group=c(group, "Cancer & tumor")}
  if (grepl("obesity|adiposity|adipose|fat|obese|lean mass", phen, ignore.case=TRUE)){group=c(group, "Adiposity")}
  if (grepl("aneurysm|Moyamoya disease", phen, ignore.case=TRUE)){group=c(group, "Cerebrovascular disorders")}
  if (grepl("weight|height|growth|size|circumference|waist|hip|anthropometric", phen, ignore.case=TRUE)){group=c(group, "Anthropometric")}
  if (grepl("bipolar|autism|neuroticism|extraversion|adventurousness|risk|mental health|neurociticism|depressive|bulimia|feeling|mood|worry", phen, ignore.case=TRUE)){group=c(group, "Psychological traits")}
  if (grepl("dementia|alzheimer|multiple sclerosis|schizophrenia|brain|subcortical|carpal tunnel|cerebrospinal|cognitive|TOURETTE|Parkinson's|antisaccade task|Amyotrophic lateral sclerosis|oxidative phosphorylation deficiency|amyloid beta|Lewy body|Neuritic|memory|Neurofibrillary tangles|fraxe|Joubert Syndrome|Meckel|Cutis laxa|amyotrophic lateral sclerosis|Microcephaly", phen, ignore.case=TRUE)){group=c(group, "Neuronal traits")}
  if (grepl("sexual|sex|menarche", phen, ignore.case=TRUE)){group=c(group, "Sexual behaviour/trait")}
  if (grepl("asthma|lung|vital capacity|expiratory|FEV1|respiratory|pulmonary|adrenergics|inhalants|Joubert Syndrome|Mucopolysaccharidosis-Plus Syndrome|ALPHA-1-ANTITRYPSIN DEFICIENCY|Cystic Fibrosis", phen, ignore.case=TRUE)){group=c(group, "Respiratoty diseases & traits")}
  if (grepl("blood|hemoglobin|platelet|hematocrit|red cell|reticulocyte|reticulocytes|mean corpuscular volume|white cells|factor vii|plasma|Alanine transaminase|serum|Lipoprotein", phen, ignore.case=TRUE)){group=c(group, "Blood traits")}
  if (grepl("eosinophil|neutrophil|basophils|monocytes|white cells|immune|reactive protein|Lymphocyte|allergy|allergic|inflammatory|crohn's disease|inflammation|sKawasaki disease|HIV-1|Systemic sclerosis|hepatitis|Cerebral amyloid deposition|angioedema|immunodeficiency|Kawasaki disease|Hyper-IgM syndrome type 5", phen, ignore.case=TRUE)){group=c(group, "Immune system & inflammation")}
  if (grepl("cholesterol|triglyceride|apolipoprotein|lipoproteins|HYPERLIPOPROTEINEMIA|HMG CoA reductase inhibitors|Dysbetalipoproteinemia|atorvastatin", phen, ignore.case=TRUE)){group=c(group, "Cholesterol traits")}  
  if (grepl("macular|ocular|retinal|retinitis|blind|visual|myopia|refractive|STARGARDT|electroretinogram|glaucoma|corneal|aniridia|Astigmatism|PI|peripheral iridotomy", phen, ignore.case=TRUE)){group=c(group, "Eye traits")}  
  if (grepl("diet|meat|vegetable|juice|fruit|beverage|fish|pork|dietary|lamb|beef|sugar|gut", phen, ignore.case=TRUE)){group=c(group, "Diet-related traits")}  
  if (grepl("smoke|smoking|alcohol|nicotine|coffee|coffeine", phen, ignore.case=TRUE)){group=c(group, "Drug-related traits")}  
  if (grepl("bone|osteo|arthritis|scoliosis|rhabdomyolysis|joint|polymyositis|brachydactyly|strength|muscular|distal|atrophy|amyotrophic lateral sclerosis|Joubert Syndrome|Cutis laxa|Boissel|Mucopolysaccharidosis-Plus Syndrome|gout|uric|urate", phen, ignore.case=TRUE)){group=c(group, "Musculoskeletal traits")}  
  if (grepl("lupus", phen, ignore.case=TRUE)){group=c(group, "Autoimmune diseases")}  
  if (grepl("cardiovascular|hypertension|heart|electrocardiogram|echocardiography|electrocardiography|pulse|carotid|stroke|beta blocking|ATHEROSCLEROSIS|Thrombosis|thromboembolism|Dysbetalipoproteinemia|coronary|calcium channel blockers|warfarin|antithrombotic|infarction|artery", phen, ignore.case=TRUE)){group=c(group, "Cardiovascular diseases")}  
  if (grepl("Heschl's gyrus|ototoxicity|hearing|deafness", phen, ignore.case=TRUE)){group=c(group, "Auditaury system")}  
  if (grepl("lifespan|birth|longevity|Mortality", phen, ignore.case=TRUE)){group=c(group, "Aging-related traits")}  
  if (grepl("dentures", phen, ignore.case=TRUE)){group=c(group, "Dental-related traits")} 
  if (grepl("kidney|urinary|glomerular filtration rate|Creatinine|Urolithiasis|diuretics|Dialysis|Meckel|Nephronophthisis", phen, ignore.case=TRUE)){group=c(group, "Renal diseases")}  
  if (grepl("hair|eyebrow|chin|morphology|skin|Microcephaly", phen, ignore.case=TRUE)){group=c(group, "Morphological traits")}  
  if (grepl("sleep|INSOMNIA|resting|chronotype|morning", phen, ignore.case=TRUE)){group=c(group, "Sleep-related traits")}  
  if (grepl("thyroid|polycystic ovary syndrome", phen, ignore.case=TRUE)){group=c(group, "Hormonal diseases & traits")}  
  if (grepl("CHYLOMICRON RETENTION DISEASE|Boissel|Mucopolysaccharidosis-Plus Syndrome|ALPHA-1-ANTITRYPSIN DEFICIENCY|Cystic Fibrosis|Sea-blue|Lipoprotein glomerulopathy|Lipoproteins|apolipoprotein|Metabolite|glucose|lipids|adiponectin|Glycogen storage disease type X|leptin|metabolic|Glycogen storage disease type X", phen, ignore.case=TRUE)){group=c(group, "Metabolic traits")}  
  if (grepl("glucose|magnesium|metabolite|lipids|sugars|iron|Cell Adhesion Molecules|potassium|Clinical laboratory measurements|Glycemic traits pleiotropy|sulfonamide|urea derivatives response", phen, ignore.case=TRUE)){group=c(group, "?")}  
  if (grepl("beta-cell|diabetes|insulin|glucosuria", phen, ignore.case=TRUE)){group=c(group, "Diabetes")}   
  if (grepl("income|intellectual|acne|TCTN2|TCF7L2|APOE|exercise test|educational|intelligence|exercises|Sasang|physical activity|chromosomal aberration frequency|Sucrose liking|quantitative traits|hemorrhage|magnesium|itch intensity", phen, ignore.case=TRUE)){group=c(group, "Other")}  
  if (grepl("Joubert Syndrome|Meckel|Cutis laxa|CHYLOMICRON RETENTION DISEASE|Mucopolysaccharidosis-Plus Syndrome|Boissel|Nephronophthisis|ALPHA-1-ANTITRYPSIN DEFICIENCY|genetic diseases|genetic disease|Cystic Fibrosis|Sea-blue|Lipoprotein glomerulopathy|Kawasaki disease|Hyper-IgM syndrome type 5", phen, ignore.case=TRUE)){group=c(group, "(Rare) genetic disorders")}  
  if (is.null(group)){group="None"}
  return(group)
}

refineGroup <- function(group){
  refGroup=NULL
  if(group %in% c("BMI", "Adiposity", "Anthropometric", "Cholesterol traits", "Cardiovascular diseases", "Diabetes", "Metabolic")){refGroup=c(refGroup, "Metabolic traits")}
  if(group %in% c("Neuronal traits", "Psychological traits")){refGroup=c(refGroup, "Neuro-psychiatric traits")}
  if(group %in% c("Cancer & tumor")){refGroup=c(refGroup, "Cancer & tumor")}
  if(group %in% c("Musculoskeletal traits", "Dental-related traits")){refGroup=c(refGroup, "Musculoskeletal traits")}
  if(group %in% c("Immune system & inflammation", "Blood traits","Autoimmune diseases")){refGroup=c(refGroup, "Immune and hematopoietic system & inflammation")}
  if(group %in% c("Eye traits", "Auditaury system", "Diet-related traits", "Sleep-related traits", "Other", "Drug-related traits", "Morphological traits", "(Rare) genetic disorders")){refGroup=c(refGroup, "Others")}
  if(group %in% c("Hormonal diseases & traits", "Renal diseases", "None", "Aging-related traits", "Respiratoty diseases & traits", "Cerebrovascular disorders")){refGroup=c(refGroup, "?")}
  return(refGroup)
}

################################################################################
# MAIN 
################################################################################
genes.lst <- c("RASAL2", "RP11-568K15.1", "FAM150B", "JADE2", "JAZF1", "HIBADH", "TMEM176A", "FASTK", "MTMR9", "MSRA",
               "SOX7", "NEIL2", "FAM167A", "AF131216.5", "RP11-981G7.2", "RP11-981G7.6", "RP1L1", "AF131215.9", "AF131215.2", "AF131215.8",
               "RP11-297N6.4", "GLIS3", "SIGMAR1", "TCF7L2", "TSKU", "WSCD2", "UNG", "DIABLO", "ARL6IP4", "RP11-197N18.2", "TCTN2",
               "RILPL2", "CDK2AP1", "SERPINA1", "FTO", "IRX3", "APOE", "ZNF226", "SNRPD2")
length(genes.lst)

hc.genes <- c("JADE2", "JAZF1", "HIBADH", "TCF7L2", "WSCD2", "DIABLO", "IRX3", "APOE") # "FTO"
genes.lst <- hc.genes

genes.dt <- data.table()
genes.summary  <- data.table()
for (g in genes.lst){
  if(length(getGroupedPhenotype(g))==2){
    tmp <- getGroupedPhenotype(g)
    genes.summary <- rbind(genes.summary, tmp[[2]])
    genes.dt <- rbind(genes.dt, tmp[[1]])
  }
}

genes.summary <- genes.summary[order(-score)]
print(genes.summary[, .(gene, topGroup)])

# --------------------- Plot grouping stats ----------------------
ggplot(as.data.table(table(genes.summary$topGroup)), aes(x=reorder(V1, -N), y=N, fill=N))+
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Likely effector genes count of top phenotype", 
       x="Phenotype group",
       y="Gene count")
ggsave("~/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/hc_genes_phenotypes.png")

ggplot(as.data.table(table(unlist(genes.dt$group))), aes(x=reorder(V1, -N), y=N, fill=N))+
# ggplot(as.data.table(table(unlist(genes.dt$refGroup))), aes(x=reorder(V1, -N), y=N, fill=N))+
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Likely effector genes count of phenotype grouping", 
       x="Phenotype group",
       y="Gene count")
ggsave("~/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/hc_genes_top_phenotypes_raw.png")

################################################################################
# Gene overlap 
################################################################################

# --------------------- Get GO genes ----------------------
GO.genes <- read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/Boer_et_al_Suppl-Table_10.xlsx") %>% as.data.table() %>% .[X18>=1 & X19!="Gene", X19] %>% unique()
length(GO.genes)
in.GO <- intersect(genes.lst, GO.genes)
length(in.GO)

# --------------------- Get GO high-confidence genes ----------------------
GO.topgenes <- read.xlsx("C:/Users/ana.arruda/Documents/OA_T2D/GWAS_eQTL_colocalization/Boer_et_al_Suppl-Table_10.xlsx") %>% as.data.table() %>% .[X18>=3 & X19!="Gene", X19] %>% unique()
length(GO.topgenes)
in.GO.top <- intersect(hc.genes, GO.topgenes)
length(in.GO.top)

# --------------------- Get DIAMANTE genes ----------------------
# https://t2d.hugeamp.org/method.html?trait=t2d&dataset=egls
DIAMANTE.genes <- vector()
DIAMANTE.genes <- Filter(Negate(is.na), lapply(content(GET("https://kp4cd.org/egldata/dataset?dataset=egls&trait=t2d"))$data, 
                         function(i) {DIAMANTE.genes <- ifelse(i$Score>=1, c(DIAMANTE.genes, i$gene), DIAMANTE.genes); return(DIAMANTE.genes)})) %>% as.character()

length(DIAMANTE.genes)
in.DIAMANTE <- intersect(genes.lst, DIAMANTE.genes)
length(in.DIAMANTE)

# --------------------- Get DIAMANTE high-confidence genes ----------------------
DIAMANTE.topgenes <- vector()
DIAMANTE.topgenes <- Filter(Negate(is.na), lapply(content(GET("https://kp4cd.org/egldata/dataset?dataset=egls&trait=t2d"))$data, 
                                               function(i) {DIAMANTE.topgenes <- ifelse(i$Score>=3, c(DIAMANTE.topgenes, i$gene), DIAMANTE.topgenes); return(DIAMANTE.topgenes)})) %>% as.character()
length(DIAMANTE.topgenes)
in.DIAMANTE.top <- intersect(hc.genes, DIAMANTE.topgenes)
length(in.DIAMANTE.top)

# --------------------- Plot overlap ----------------------
venn.diagram(
  x = list(DIAMANTE.topgenes, GO.topgenes, hc.genes),
  category.names = c("T2D", "OA", "T2D+OA"),
  col = "transparent",
  fill = c("green", "red", "blue"),
  alpha = 0.30,
  print.mode=c("raw","percent"),
  filename = "~/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/hc_overlap_T2D_OA.png",
  imagetype="png",
  output=TRUE)

venn.diagram(
  x = list(DIAMANTE.genes, GO.genes, genes.lst),
  category.names = c("T2D", "OA", "T2D+OA"),
  col = "transparent",
  fill = c("green", "red", "blue"),
  alpha = 0.30,
  print.mode=c("raw","percent"),
  filename = "~/OA_T2D/GWAS_eQTL_colocalization/comparison_plots/overlap_T2D_OA.png",
  imagetype="png",
  output=TRUE)

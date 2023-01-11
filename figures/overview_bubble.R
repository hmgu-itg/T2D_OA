library(data.table)
library(ggplot2)
library(RColorBrewer)

regions.dt <- data.table(region=c(rep(1,5), rep(2,4), 3, 4, rep(5,2), 6, rep(7,4), rep(8,4), 9, rep(10,7), rep(11,4), 12, 13, rep(14,4), rep(15,2), 16, rep(17,7), 18),
                         chr=c(rep(1,5), rep(2,4), 4, 4, rep(5,2), 6, rep(7,4), rep(7,4), 8, rep(9,7), rep(9,4), 10, 11, rep(12,4), rep(12,2), 14, rep(16,7), 19),
                         pos=c(rep(219741820,5), rep(c(422144, 417167),2), 145621328, 45186139, rep(133864599,2), 50788778, 28192280, 28200097, 28192280, 28192280, rep(150537635,4),
                               10772644, rep(4291928,7), rep(34074476,4), 114758349, 76505202, rep(108629780,4), 123450765, 123681222, 94838142, 53816275, rep(53800200,3),
                               53809247, 53816275, 53828066, 45411941),
                         OA.phen=c("All", "KneeHip", "THR", "TKR", "TJR", "All", "TKR", "THR", "TJR", "THR", "Knee", "Knee", "KneeHip", "Knee", "KneeHip", "Hip", "THR", "TJR", 
                                   "Knee", "KneeHip", "TKR", "TJR", "KneeHip", "All", "Knee", "KneeHip", "Hip", "THR", "TKR", "TJR", "All", "KneeHip", "TKR", "TJR", "Knee", "KneeHip", 
                                   "All", "Knee", "KneeHip", "TJR", "Knee", "TKR", "TJR", "All", "Knee", "Hip", "KneeHip", "THR", "TKR", "TJR", "All"),
                         PP4=c(0.988,0.857,0.986,0.87,0.958,0.965,0.808,0.85,0.869,0.99,0.81,0.936,0.89,0.836,0.979,0.976,0.965,0.948,0.993,.918,0.95,0.932,0.839,0.989,0.998,1,0.903,0.968,
                               0.933,0.999,0.958,0.946,0.978,0.814,0.936,0.883,0.983,0.999,0.999,0.904,0.822,0.831,0.992,0.935,0.935,0.936,0.932,0.925,0.922,0.934,0.997),
                         credset=c(7,11,6,8,7,28,127,94,126,2,5,4,4,20,5,5,5,5,2,2,7,9,14,1,1,1,1,1,1,1,1,1,1,1,3,23,2,2,2,2,130,110,2,59,37,44,60,62,50,54,1))
# regions.dt[, OA.phen.group:=ifelse(OA.phen %in% c("All", "KneeHip", "TJR"), "both", ifelse(OA.phen %in% c("Knee", "TKR"), "knee", "hip"))]
regions.dt[, OA.phen.group:=ifelse(OA.phen %in% c("Knee", "TKR"), ifelse(OA.phen %in% c("Hip", "THR"), "both", "knee"), ifelse(OA.phen %in% c("Hip", "THR"), "hip", "both"))]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

ggplot(regions.dt, aes(credset, PP4, fill=OA.phen.group, size=credset)) +
  geom_point(color="black", pch=21) +
  # geom_point(color="black", pch=21, alpha=0.8) +
  labs(y="PP4", 
       x="Number of variants in 95% credible set", 
       # title="Overview of colocalized regions between T2D and OA", 
       fill="OA phenotype",
       size="Number of variants \n in 95% credible set") +
  theme_bw() +
  # scale_fill_manual(values=cbPalette) +
  scale_fill_brewer(palette="Paired") +
  scale_size_continuous(range = c(2, 10)) +
  guides(fill=guide_legend(override.aes = list(size=8)))
ggsave("/project_data/regions_figure.png", width=7.3, height=4.4, dpi=320)

my_colors <- RColorBrewer::brewer.pal(9, "RdPu")[c(3,5,8)]

ggplot(regions.dt, aes(credset, PP4, fill=OA.phen.group, size=credset)) +
  geom_point(color="black", pch=21) +
  labs(y="PP4", 
       x="Number of variants in 95% credible set", 
       fill="OA phenotype",
       size="Number of variants \n in 95% credible set") +
  theme_bw() +
  # scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values = my_colors) +
  scale_size_continuous(range = c(2, 10), guide="none") +
  guides(fill=guide_legend(override.aes = list(size=8))) +
  theme(legend.position = c(0.88, 0.69)) +
  annotate("text", label="Size = number of variants in 95% credible set", x=88, y=1)
# ggsave("/project_data/regions3_figure.png", width=7.7, height=4.7, dpi=320)
ggsave("/project_data/overlap_T2D_OA/paper/figures/coloc_regions_bubbles_pink.png", width=5.405, height=3.68, dpi=320, units="in")

##############################################################################
#---------------------- Distance to lead variant in bp ----------------------#
##############################################################################
# Look at distance in base pairs between the gene closest to the lead variant and the HC gene in each region
library(biomaRt)
library(IRanges)
library(GenomicRanges)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version="GRCh37")

hc.genes <- as.data.table(getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position"),
                                filters="external_gene_name",
                                values=as.list(fread("/project_data/overlap_T2D_OA/hc_genes.txt", header=FALSE))$V1,mart=ensembl)) %>% .[nchar(chromosome_name)<3]
hc.genes <- hc.genes[order(external_gene_name)]
hc.genes$region <- c(18,15,11,1,17,10,7,17,5,9,13,18,8,18,8,12,14,8,14)

merged.dt <- merge(regions.dt, hc.genes, by="region", all.x=TRUE, allow.cartesian=TRUE)
setnames(merged.dt, c(colnames(merged.dt)[1:7], "hc.id", "hc.name", "hc.chr", "hc.start", "hc.end"))
merged.dt[, hc.chr:=NULL]
merged.dt[, dist.lead:=ifelse(pos>hc.end, pos-hc.end, ifelse(pos<hc.start, hc.start-pos, 0))]

nearest.gene <- function(chr, pos){
  df.start=ACME::findClosestGene(paste0("chr",chr), pos, genome="hg19", position="txStart")
  df.end=ACME::findClosestGene(paste0("chr",chr), pos, genome="hg19", position="txEnd")
  if(abs(df.start$Distance) < abs(df.end$Distance)) df <- df.start
  else df <- df.end
  return(unique(df$geneName))
}

merged.dt[, nearest.name:=nearest.gene(chr, pos), by=seq_len(nrow(merged.dt))]
merged.dt[nearest.name=="LYPLAL1-AS1", nearest.name:="LYPLAL1"]

nearest.genes.dt <- as.data.table(getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position"),
                                  filters="external_gene_name",
                                  values=unique(merged.dt$nearest.name),mart=ensembl)) %>% .[nchar(chromosome_name)<3]
setnames(nearest.genes.dt, c("nearest.id", "nearest.name", "nearest.chr", "nearest.start", "nearest.end"))
nearest.genes.dt[, nearest.chr:=NULL]

merged.dt <- merge(merged.dt, nearest.genes.dt, by="nearest.name", all.x=TRUE, allow.cartesian=TRUE)
merged.dt[, dist.nearest:=ifelse(pos>hc.end, pos-hc.end, ifelse(pos<hc.start, hc.start-pos, 0))]
merged.dt[dist.nearest!=0, dist.nearest:=ifelse(hc.end<nearest.start, nearest.start-hc.end, hc.start-nearest.end)]
merged.dt[, hc.nearest:=ifelse(dist.nearest==0, "yes", "no")]
# setorder(merged.dt, c(colnames(merged.dt)[2:12], colnames(merged.dt)[1], colnames(merged.dt)[13:17]))

ggplot(merged.dt[!is.na(external_gene_name)], aes(dist.hc, PP4, color=hc.nearest)) +
  geom_point() +
  labs(x="Distance in bp to hc gene",
       color="Is the hc gene the \n nearest gene?") +
  theme_bw() +
  scale_color_brewer(palette="Set1")

ggplot(merged.dt[!is.na(hc.name)], aes(x=dist.lead, y=hc.name, color=hc.nearest)) +
  geom_point() +
  theme_bw() + 
  # scale_color_brewer(palette="Paired") +
  scale_color_manual(values=RColorBrewer::brewer.pal(8, "Paired")[c(4,2)]) +
  ylab("") +
  xlab("Distance in bp to lead variant") +
  labs(color="Is the HC gene the nearest gene?") +
  theme(legend.position="top")
ggsave("/project_data/overlap_T2D_OA/paper/figures/nearest_gene_lead_variant.png", width=6, height=4, dpi=320)

# Distance to nearest gene in bp
ggplot(merged.dt[!is.na(hc.name)], aes(x=dist.nearest, y=hc.name, color=hc.nearest)) +
  geom_point() +
  theme_bw() + 
  # scale_color_brewer(palette="Paired") +
  scale_color_manual(values=RColorBrewer::brewer.pal(8, "Paired")[c(4,2)]) +
  ylab("") +
  xlab("Distance in bp to nearest gene") +
  labs(color="Is the HC gene the nearest gene?") +
  theme(legend.position="top")
ggsave("/project_data/overlap_T2D_OA/paper/figures/nearest_gene.png", width=6, height=4, dpi=320)




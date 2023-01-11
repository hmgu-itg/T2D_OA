library("msigdbr")
library("dplyr")
library("clusterProfiler")
library("wesanderson")
library("stringr")
library("ggplot2")
library("cowplot")
library(data.table)
library(org.Hs.eg.db)

#--------------------- Load databases of interest
#KEGG
m_t2g_C2_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()

#REACTOME
m_t2g_C2_REACTOME <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()

# WIKIPATHWAYS
m_t2g_C2_WIKIPATHWAYS <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()

#GO (Here I am interested in the subcategories Molecular function and Biological Process, the third category which is very general is Cellular Component CC)
m_t2g_C5_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()
m_t2g_C5_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()
m_t2g_C5_CC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()
m_t2g_C5_ALL <- rbind(m_t2g_C5_BP,m_t2g_C5_MF,m_t2g_C5_CC)

#HPO
m_t2g_C5_HPO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO") %>%  dplyr::select(gs_name,gene_symbol) %>% as.data.frame()

#--------------------- Function to run ORA
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

run_ORA <- function(gene_list, universe, gene_set, min_size, max_size, pval, GOonly=FALSE, GOlevel=NA) {
  
  #gene_list: character vector of gene symbols of interest
  #universe: character vector of all genes - eg. all genes of an analysis, all homo sapiens genes or all genes expressed in a specific tissue
  #gene_set: a data.frame (generated above using the msigdbr function having as columns 
  #gs_name: temrs in a gene set 
  #gene_symbol: gene symbol contained in the specific term
  #min_size: minimal size of genes annotated for testing in a term
  #max_size: maximum size of genes annotated for testing
  #pval: pvalue cutoff   
  
  #Returns:
  #The output is a data.frame with columns:
  #ID: TermID in MgSigDB
  #Description: same as ID
  #GeneRatio:k/n
  #k = size of the overlap of the vector of gene names you input with the specific geneset (eg A KEGG pathway term)
  #n = size of the overlap of the vector of gene names you input with all the members of the collection of genesets (eg the KEGG collection)
  #BgRatio:M/N  
  #M = size of the geneset (eg size of a KEGG pathway)
  #N = size of all of the unique genes in the collection of genesets (example the KEGG collection)
  #pvalue: pvalue: Raw p-value 
  #p.adjust: FWER-adjusted p-value
  #qvalue: FDR q-value for each gene set 
  #geneID: gene symbols from your list of genes of interest found within a term (eg. a specific KEGG pathway)
  #Count: number of genes from your list included in the term (eg. a specific KEGG pathway)
  
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    
    #Exclude duplicate gene names
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  
  if(!GOonly) {
    #Perform ORA using the enricher function from clusterProfiler
    #enrichr() performs a hypergeometric test comparing the set of "significant" genes against the "universe" (or background) genes.
    OraRes <- enricher(gene = gene_list,
                       universe = universe,
                       minGSSize = min_size,
                       maxGSSize = max_size,
                       pvalueCutoff = pval,
                       pAdjustMethod = "fdr",
                       TERM2GENE =  gene_set)
    output = OraRes %>%
      as.data.frame() %>%
      dplyr::filter(p.adjust < !!pval) %>%
      dplyr::arrange(log(p.adjust))
  }
  else{
    OraRes <- enrichGO(gene = gene_list,
                       keyType = "SYMBOL",
                       ont = "BP",
                       universe =  universe,
                       OrgDb = 'org.Hs.eg.db',
                       minGSSize = min_size,
                       maxGSSize = max_size,
                       pvalueCutoff = pval,
                       pAdjustMethod = "fdr")
    
    # output <- OraRes %>%
    output <- gofilter(OraRes, level=GOlevel) %>%
      # output <- dropGO(OraRes, level=5) %>%
      as.data.frame() %>%
      dplyr::filter(p.adjust < !!pval) %>%
      dplyr::arrange(log(p.adjust))
  }
  
  return(output)
}

plot.ORA <- function(ora_all_datasets, my.set, max_n=10){
  # plot_terms_ora <- c(ora_all_datasets[[1]][1:max_n,]$ID, ora_all_datasets[[2]][1:max_n,]$ID, ora_all_datasets[[3]][1:max_n,]$ID)
  plot_terms_ora <- c(ora_all_datasets[[1]][1:max_n,]$ID)
  
  options(repr.plot.width = 10, repr.plot.height = 9)
  
  #Here I plot the most significant terms top 10 of categories (GO, REACTOME, KEGG)
  ggplot(data = subset(ora_all_datasets_df, ID %in% plot_terms_ora),
         aes(x = -log(p.adjust),
             y = reorder(term_to_plot_ora,-log(p.adjust)),
             fill = Dataset))+geom_col()+ 
    facet_grid(rows = vars(Dataset),scales = "free_y",space = "free_y") +
    theme_bw()+ 
    scale_fill_manual(values=wes_palette(n=3, name="Moonrise2")) + 
    ylab("") +theme(axis.text  = element_text(colour = "black", size = 15),
                    axis.title = element_text(colour = "black", size = 15),
                    legend.position = "none")+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
  
  ggsave(paste0("/project_data/overlap_T2D_OA/pathway_analysis/", names(my.set), ".png"))
}

my.ORA <- function(min_, max_, lst, my.database, GOonly=FALSE, GOlevel=NA) {
  pathway.df <- data.table()
  
  for (my.set in names(lst)){
    ora_all_datasets <- lapply(my.database, function(i){run_ORA(gene_list =  lst[[my.set]],                                                              
                                                                universe = , #all.genes$all, 
                                                                gene_set = i, 
                                                                min_size = min_, 
                                                                max_size = max_, 
                                                                pval = 0.05,
                                                                GOonly=GOonly,
                                                                GOlevel=GOlevel)})
    #-------------------- Save results
    #Make a character vector of the names of gene sets I tested
    setnames <- c(names(my.database))
    
    #Add a column Dataset indicating the respective gene set to all the dataframes generated from the individual datasets separate them for plotting 
    for(i in 1:length(setnames)){   
      ora_all_datasets[[i]]$Dataset <- rep(setnames[i], nrow(ora_all_datasets[[i]]))
    }
    
    #Make a big dataframe for all result along sets by binding the rows
    ora_all_datasets_df <-  rbindlist(ora_all_datasets)
    
    #Make terms lower case
    term_to_plot_ora <- sub(".*?_",'',tolower(ora_all_datasets_df$Description))
    #Substitute the underscore in all terms with gaps
    term_to_plot_ora  <- gsub( "_"," ", term_to_plot_ora )
    #Capitalise first letter
    term_to_plot_ora <- sapply(term_to_plot_ora ,function(i){CapStr(i)})
    
    #Add the "beautified" term names to my big dataframe 
    ora_all_datasets_df$term_to_plot_ora <- unlist(term_to_plot_ora)
    
    pathway.df <- rbind(pathway.df, as.data.table(ora_all_datasets_df) %>% .[, .(GeneSet=paste0(length(lst[[my.set]]), " genes in the ", my.set, " set"), Pathway=term_to_plot_ora, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Dataset)])
    # plot.ORA(ora_all_datasets, my.set)
  }
  return(pathway.df)
}

#-------------------- Run for hc genes
read_gene_sets <- function(set.name){
  tmp <- read.csv(paste0("/project_data/overlap_T2D_OA/pathway_analysis/", set.name, ".txt"), header=FALSE)
  tmp <- clusterProfiler::bitr(tmp[[1]], fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db) %>% setnames(., c("ENSEMBL", set.name))
  return(tmp[set.name])
}

gene.sets <- sapply(list("potential", "likely", "likely_knee", "likely_hip", "hc", "hc_knee", "hc_hip"), read_gene_sets)
# gene.sets <- sapply(list("hc", "hc_knee", "hc_hip"), read_gene_sets)
pathway.df <- my.ORA(30,2000, gene.sets, my.database=list(WikiPathways=m_t2g_C2_WIKIPATHWAYS, HPO=m_t2g_C5_HPO, REACTOME=m_t2g_C2_REACTOME, KEGG=m_t2g_C2_KEGG))
pathway.GO <- my.ORA(30,2000, gene.sets, my.database=list(GO=m_t2g_C5_BP), GOonly=TRUE, GOlevel=5)
pathway.df <- rbind(pathway.df, pathway.GO)
pathway.df <- pathway.df[!is.null(Pathway) & !is.na(geneID)]

fwrite(pathway.df[!is.na(pvalue)], "/project_data/overlap_T2D_OA/pathway_analysis/pathway_table.csv")

# potential, likely.knee (53): 150-500 / 300-1000
# potential, hc (16), hc.hip (7): 5-5000/ 30-5000 / 10-10000
# potential, hc (16), hc.knee (15): 200-5000 / 1000-5000 / 300-3000 / 200-10000
# hc, hc.knee: 1500-5000
# potential, hc: 35-5000 / 100-5000

for (mi in c(3, seq(10,40,10), seq(50,500,50), seq(600,1500,100))) {
  for (ma in c(100, 250, seq(500,10000,500))) {
    # assign("last.warning", NULL, envir = baseenv())
    if (mi<ma){
      tryCatch({
        dt.hc <- my.ORA(mi, ma, list(hc.hip.genes))
        # if (is.null(warnings())) res[[paste(mi, ma, sep="-")]] <- dt
        if(dt.hc[is.na(p.adjust),.N]==0) res[[paste(mi, ma, sep="-")]] <- dt.hc
        
        dt <- my.ORA(mi, ma, list(likely.hip.genes))
        if(dt[is.na(p.adjust),.N]==0) res[[paste(mi, ma, sep="-")]] <- dt
      }, error=function(e){})
    }
  }
}

my.ORA(5,5000, gene.sets, my.database=list(GO=m_t2g_C5_BP,WikiPathways=m_t2g_C2_WIKIPATHWAYS, HPO=m_t2g_C5_HPO, REACTOME=m_t2g_C2_REACTOME, KEGG=m_t2g_C2_KEGG))
my.ORA(4,1000, gene.sets, my.database=list(GO=m_t2g_C5_BP), GOonly=TRUE, GOlevel=5)

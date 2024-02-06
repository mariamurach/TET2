library(DESeq2)
library(tidyverse)
library(readxl)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(biomaRt)
library(pheatmap)

  font <- "Georgia"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.grid = element_line(colour = "grey95"),
      text = element_text(size = 8),
      plot.title = element_text(size = 8, margin = margin(rep(2,4))),
      axis.line = element_blank(),
      legend.key = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      complete = TRUE,
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white",
                                      colour = "grey20"), 
      legend.position="right"
    )
}
mor_normalization = function(data){
  library(dplyr)
  library(tibble)
  
  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}

safe_colorblind_palette <- c("#88CCEE",  "#DDCC77", "#CC6677","#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

####
mart <- useMart("ensembl","mmusculus_gene_ensembl")
####
plot_enrichments2 <- function(data_full, cell, padj_name = "p.adjust", w = 6, h = 6){
  data <- data_full %>% as_tibble()
  if(nrow(data) == 0) return()
  x <- ifelse(class(data_full) == "gseaResult", "enrichmentScore", "GeneRatio")
  if(x == "GeneRatio") {
    data$GeneRatio <- sapply(data$GeneRatio, function(x) eval(parse(text=x))) %>% unname
    Cnt = "Count"
    finame = paste0("geneEnrichment_", data_full@ontology)
  }else{
    Cnt = "setSize"
    finame = paste0("pathwayEnrichment_", data_full@setType)
  }
  data <- data %>% slice(1:10) 
  plot <- data %>%  dplyr::filter(p.adjust < 0.05) %>%
    head(40) %>% dplyr::mutate(Description = str_remove(Description,"- Mus(.*)")) %>% 
    ggplot(aes(get(x), reorder(Description, -get(padj_name)))) + 
    geom_point(colour="black", size = 5) +
    geom_point(aes(color = get(padj_name)), size = 4) +
    scale_color_gradient(high="blue", low="red") +
    scale_size(range = c(1, 4)) + theme_bw() %+replace% 
    theme(legend.text = element_text(size = 9), 
          legend.title = element_text(size = 9),
          legend.key.size = unit(1,"line")) + xlab(x) + labs(color = "p.adjust", y = "Description") +
    guides(size = "none") + ggtitle(paste0(cell)) 
  
  
  plot2 <- data %>%  dplyr::filter(p.adjust < 0.05) %>%
    head(40) %>% dplyr::mutate(Description = str_remove(Description,"- Mus(.*)")) %>% 
    ggplot(aes(get(x), reorder(Description, -get(padj_name)))) + 
    geom_bar(stat = "identity", aes(fill = get(padj_name))) + 
    scale_fill_gradient(high="blue", low="orange") +
    scale_size(range = c(1, 4)) + theme_bw() %+replace% 
    theme(legend.text = element_text(size = 9), 
          legend.title = element_text(size = 9),
          legend.key.size = unit(1,"line")) + xlab(x) + labs(fill = "p.adjust", y = "Description") + 
    guides(size = "none") + ggtitle(paste0(cell)) 
  
  finame_png <- paste0(res, cell,"/",finame,".png")
  finame_png1 <- paste0(res, cell,"/", "1_", finame,".png")
  
  finame_csv <- paste0(res, cell,"/",finame,".csv")
  
  ggsave(filename = finame_png1, plot = plot, width = w, height = h, units = "in")
  ggsave(filename = finame_png, plot = plot2, width = w, height = h, units = "in")
  write_csv(x = data, file = finame_csv)
  #print(finame)
  #return(plot)
  return(plot2)
}

wd <- "analysis"
setwd(wd)
res <- "analysis/results"
createResDir <- function(dir_name = ""){
  if (!dir.exists(paste0(res, dir_name))){
    dir.create(paste0(res, dir_name))
  }else{
    print("dir exists") 
  }}
createResDir()
#Read in counts and metadata
samples <- read_csv("data_table.csv") %>% dplyr::select(sample, library)
write_csv(samples, "samples.csv")

flist <-list.files(path = "counts",
          pattern = "*counts.txt$",
          full.names = T)

myfiles = lapply(flist, read_tsv, skip  = 1)
counts <- purrr::reduce(myfiles, left_join)
counts_og <- counts
counts$Geneid <- counts$Geneid %>% str_replace("\\.(.*)", "")
rownames(counts) <- counts$Geneid


#Setup column names as actual samples names
test <- names(counts[-(1:6)])
test <- test %>% str_replace("(.*?)\\/", "") %>% str_remove("Aligned.sortedByCoord.out.bam")
smpls <- data.frame(library = test)
test <- samples$sample[match(smpls$library, samples$library)]
test <- c(names(counts[1:6]), test)
names(counts) <- test
counts_og <- counts

#setup metadata
samples$condition <- map2(str_split(samples$sample, "_"),1, pluck) %>% unlist %>% factor
samples$cell_type <- map2(str_split(samples$sample, "_"),2, pluck) %>% unlist %>% factor
write_csv(samples, "samples.csv")

rownames(samples) <- samples$sample
cts <- counts[, rownames(samples)]
rownames(cts) <- rownames(counts)
coldata <- samples
coldata$condition <- factor(coldata$condition, levels = c("WT", "KO"))
all(coldata$sample %in% colnames(cts))
all(coldata$sample == colnames(cts))

### TODO: could be put in a for loop
B1a <- coldata %>% dplyr::filter(str_detect(cell_type, "B1a")) 
B1b <- coldata %>% dplyr::filter(str_detect(cell_type, "B1b")) 
B2 <- coldata %>% dplyr::filter(str_detect(cell_type, "B2")) 
## setup seperate tables for each cell type
B2_cts <- cts[B2$sample]
rownames(B2_cts) <- rownames(counts)
B1a_cts <- cts[B1a$sample]
rownames(B1a_cts) <- rownames(counts)
B1b_cts <- cts[B1b$sample]
rownames(B1b_cts) <- rownames(counts)

B2_cts <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(B2_cts)), 
                               colData = B2)
dds_B2 <- DESeqDataSet(B2_cts, design = model.matrix(~condition, B2))

B1a_cts <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(B1a_cts)), 
                               colData = B1a)
dds_B1a <- DESeqDataSet(B1a_cts, design = model.matrix(~condition, B1a))

B1b_cts <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(B1b_cts)), 
                               colData = B1b)
dds_B1b <- DESeqDataSet(B1b_cts, design = model.matrix(~condition, B1b))

### Differential analysis
for(cell in c("B1a", "B1b")) {
  res <- "db_check/"
  createResDir(cell)
  dt <- get(paste0("dds_", cell))
  dt$condition <- factor(dt$condition, levels = c("WT", "KO"))
  featureData <- data.frame(gene = rownames(dt))
  mcols(dt) <- DataFrame(mcols(dt), featureData)
  mcols(dt)
  keep <- rowSums(counts(dt)) >= 10
  dt <- dt[keep, ]
  data <- DESeq(dt, parallel = T)
  #With P adj < 0.05
  data_res <- results(data)
  resLFC <-
    lfcShrink(data,
              coef = resultsNames(data)[2],
              type = "apeglm",
              parallel = T)
  lfc <- as_tibble(resLFC)
  lfc$gene <- rownames(resLFC)
  
  ##### lfc shrinkage - trying to retrieve as much gene symbols connected to gene name as possible
  ensemble2gene <- getBM(attributes=c("external_gene_name", "ensembl_gene_id"),
                         filters = "ensembl_gene_id",
                         values = lfc$gene, 
                         mart = mart)
  lfc <- lfc %>% left_join(ensemble2gene, by = c("gene" = "ensembl_gene_id"))
  lfc$symbol2 <- lfc$external_gene_name
  last <- function(x){x[[length(x)]]}
  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                     keys = lfc$gene, 
                                     column = c("SYMBOL"), keytype = "ENSEMBL", multiVals=last)
  lfc$symbol <- symbols
  ##### PHEATMAP
  sig_index <- which(lfc$padj < 0.05)
  scaled_counts <- mor_normalization(as.data.frame(assay(data)))
  scaled_counts$gene <- rownames(scaled_counts)
  #scaled_counts <- scaled_counts %>% left_join(sig_symbols, by = c("gene" = "ENSEMBL"))

  sig_genes <- scaled_counts[sig_index,] %>% as_tibble()
  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                   keys = sig_genes$gene, 
                                   column = c("SYMBOL"),
                                   keytype = "ENSEMBL",
                                   multiVals=last)
  sig_genes$SYMBOL = symbols
  
  sig_genes <- sig_genes %>% drop_na(SYMBOL)
  sig_genes_df <- sig_genes %>% dplyr::select( 1:7 ) %>% as.data.frame()
  names(sig_genes_df) <- names(sig_genes_df) %>% str_sub(start = 1, end = 2)
  rownames(sig_genes_df) <- sig_genes$SYMBOL
  ph <- pheatmap(t(sig_genes_df), 
           cluster_rows = T, 
           cluster_cols = T, 
           scale = "column",
           show_rownames = T, 
           show_colnames = T, 
           fontsize_row = 9,
           fontsize_col = 8)
  if(cell == "B1b")
  {
    ggsave(plot = ph, filename = paste0(res, cell, "/", cell, "_heatmap.png") , w = 9.4, h = 4, dpi = 300)
  }
  ggsave(plot = ph, filename = paste0(res, cell, "/", cell, "_heatmap.png") , w =19, h = 4, dpi = 300)
  
  lfc <- lfc %>% arrange(desc(log2FoldChange))
  
  entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = lfc$symbol, columns = c("ENTREZID"), keytype = "SYMBOL")
  plotted_genes <- lfc %>% filter(padj < 0.05 & str_detect(symbol, "Ighv*") & abs(log2FoldChange) > 1 ) %>% 
    arrange(desc(log2FoldChange), desc(padj)) %>% pull(symbol)
   p <- EnhancedVolcano(lfc,
                  lab = lfc$symbol,
                  x = 'log2FoldChange',
                  y = 'padj', 
                  selectLab = plotted_genes,
                  axisLabSize = 16,
                  subtitleLabSize = 0,
                  captionLabSize = 0, 
                  titleLabSize = 16,
                  title = paste0(cell, ", ", "KO vs WT"), 
                  subtitle = "",
                  labSize = 7, 
                  xlab = bquote(~Log[2] ~ "(TET2KO/WT)"),
                  ylab = bquote(~-Log[10] ~ italic(p[adjust])),
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black'
                
  ) 
  p
  ggsave(plot = p, filename =  paste0(res, cell, "/", "volcano.pdf"), width = 9, height = 6, units = "in", dpi = 300)
  
  entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = lfc$symbol, columns = c("ENTREZID"), keytype = "SYMBOL")
  
  # gsea diff analysis
  lfc <- lfc %>% left_join(entrez, by = c("symbol" = "SYMBOL"))
  genes <- lfc %>% filter(padj < 0.05) 
  
  gns <- tibble(ENTREZ = genes$ENTREZID, 
                SYMBOL = genes$symbol, 
                log2FC = genes$log2FoldChange,
                padj = genes$padj)
  gns_write <- gns
  gns_write$log2FC <- round(gns_write$log2FC, 2)
  write_csv(gns_write, paste0(res, cell, "/", cell, "_diff_genes", ".csv"))
  #Hallmark
  .path.lfc <-  lfc %>% dplyr::select(log2FoldChange, symbol, symbol2) %>% unique
  ranks <- lfc$log2FoldChange
  
  names(ranks) <- lfc$ENTREZID
  head(ranks)
  
  for (path in c("C2", "C3", "C4", "C5", "C6", "C7", "H")) {
    pathways = msigdbr(species = "mouse", category = path)
    pathways = split(x = pathways$gene_symbol, f = pathways$gs_name)
    fgseaRes <-
      fgsea(
        pathways,
        ranks,
        minSize = 15,
        maxSize = 500,
        eps = 0
      )
    fgseaRes <- fgseaRes %>% arrange(padj)
    createResDir(cell)
    fname <-
      paste0(res, cell, "/genes_", cell, "_msigdb_", path, ".csv")
    write_csv(fgseaRes, fname)
  }

  ggo_BP <- enrichGO(
    gene     = gns$ENTREZ,
    OrgDb    = org.Mm.eg.db,
    ont      = "ALL",
    readable = TRUE
  )
  plot <- plot_enrichments2(ggo_BP, cell = cell, w = 7, h = 3 )
  ggo_BP <- ggo_BP %>% as_tibble
  write_csv(ggo_BP, paste0(res, cell, "/", cell, "_geneEnrichmentGO.csv"))
  
  
  ranks <- lfc$log2FoldChange
  
  names(ranks) <- lfc$ENTREZID
  head(ranks)
  ranks_og <- ranks
  ranks2 <- ranks_og
  
  #entrez as names
  tmp <- which(is.na(names(ranks)) | names(ranks) == "")
  ranks <- ranks[-tmp]
  ranks = ranks[!duplicated(names(ranks))]
  
  #symbol as names
  names(ranks2) <- lfc$symbol2
  tmp <- which(is.na(names(ranks2)) | names(ranks2) == "")
  ranks2 <- ranks2[-tmp]
  ranks2 = ranks2[!duplicated(names(ranks2))]
  ranks2 = ranks2[!duplicated(unname(ranks2))]
  
  input.save <- tibble(lfc = ranks)
  input.save$gene <- names(ranks)
  input.save <- input.save %>% arrange(desc(lfc))
  input.save <- input.save %>% mutate(rank = rank(lfc,  ties.method = "random"))
  
  write.csv(input.save, paste0(res, cell, "/", cell, "_gsea_input.csv"))
  
  
  input.save2 <- tibble(lfc = ranks2)
  input.save2$gene <- names(ranks2)
  input.save2 <- input.save2 %>% arrange(desc(lfc))
  input.save2 <- input.save2 %>% mutate(rank = rank(lfc,  ties.method = "random"))
  write.csv(input.save2, paste0(res, cell, "/", cell, "_gsea_input.csv"))
  gse_BP <- gseGO(gene     = sort(ranks, decreasing = T),
                  OrgDb    = org.Mm.eg.db,
                  ont      = "ALL", 
                  keyType = "ENTREZID")

  gse_BP2 <- gseGO(gene     = sort(ranks2, decreasing = T),
                  OrgDb    = org.Mm.eg.db,
                  ont      = "ALL", 
                  keyType = "SYMBOL")
  plot <- plot_enrichments2(gse_BP2, cell = cell, w = 7, h = 3 )
  gse_BP2 <- gse_BP2 %>% as_tibble()
  write_csv(gse_BP2, paste0(res, cell, "/", cell, "_pathwaysEnrichmentGO.csv"))
  
  
  .paths_plot <- gse_BP2 %>%  dplyr::filter(p.adjust < 0.05) %>%
    head(40) %>% 
    ggplot(aes(enrichmentScore, reorder(Description, -p.adjust))) + 
    geom_bar(stat = "identity", aes(fill = p.adjust)) + 
    scale_fill_gradient(high="blue", low="orange") +
    scale_size(range = c(1, 4)) + theme_prism(base_size = 9) + xlab("enrichmentScore") + 
    labs(fill = "p.adjust", y = "Description")
    
  
  theme_bw() %+replace% 
    theme(legend.text = element_text(size = 9), 
          legend.title = element_text(size = 9),
          legend.key.size = unit(1,"line"))
  ggsave(h = 5.5, w = 7.5, plot = .paths_plot, filename = paste0(res, cell, "/", cell, "_pathwaysEnrichmentGO.png"), dpi = 300)
  
  if(cell == "B1a"){
    B1a_igs<- gse_BP2[1, ]
    B1a_senso <- gse_BP2[4, ] 
    B1a_senso2 <- gse_BP2[6, ]
    B1a_plasma_membrane<- gse_BP2[2, ]
    B1a_immune_response<- gse_BP2[3, ]
    
    
    ### Ig production
    .pathway_genes <- str_split(B1a_igs$core_enrichment, "/") %>% unlist()
    .pathway_genes2 <- str_split(B1a_immune_response$core_enrichment, "/") %>% unlist()
    .pathway_genes <- c(.pathway_genes, .pathway_genes2) %>% unique
    .gene_set <- sig_genes_df[rownames(sig_genes_df) %in%.pathway_genes,]
    .ph <- pheatmap(.gene_set, 
                   cluster_rows = T, 
                   cluster_cols = F, 
                   scale = "row",
                   show_rownames = T, 
                   show_colnames = T, 
                   fontsize_row = 9,
                   fontsize_col =7,
                   clustering_method ="ward.D", 
                   angle_col = 0)
    ggsave(plot = .ph, filename = paste0(res, cell, "/", cell, "immunoglobulin_and_immune_response.png"), 
           w  = 3.5, h = 6, dpi = 300, unit = "in")
    #### plasma
    .pathway_genes <- str_split(B1a_plasma_membrane$core_enrichment, "/") %>% unlist()
    .gene_set <- sig_genes_df[rownames(sig_genes_df) %in%.pathway_genes,]
    .ph <- pheatmap(.gene_set, 
                    cluster_rows = T, 
                    cluster_cols = F, 
                    scale = "row",
                    show_rownames = T, 
                    show_colnames = T, 
                    fontsize_row = 9,
                    fontsize_col =7,
                    clustering_method ="ward.D", 
                    angle_col = 0)
    ggsave(plot = .ph, filename = paste0(res, cell, "/", cell, "plasma.png"), 
           h  = 6, w = 3.5, dpi = 300, unit = "in")
   
    # Sensory perception of smell
    .pathway_genes <- str_split(B1a_senso2$core_enrichment, "/") %>% unlist()
    .ensembl <- AnnotationDbi::select(org.Mm.eg.db, 
                                      keys = .pathway_genes, 
                                      column = c("ENSEMBL"), 
                                      keytype = "SYMBOL")
    .gene_set <- scaled_counts[scaled_counts$gene %in% .ensembl$ENSEMBL,]
    
    .gene_set <- .gene_set %>% left_join(.ensembl, by = c("gene" = "ENSEMBL"))
    rownames(.gene_set) = .gene_set$SYMBOL
    .gene_set <- .gene_set %>% dplyr::select( 1:7 ) %>% as.data.frame()
    names(.gene_set) <- names(.gene_set) %>% str_sub(start = 1, end = 2)
    .ph <- pheatmap(.gene_set, 
                    cluster_rows = T, 
                    cluster_cols = F, 
                    scale = "row",
                    show_rownames = T, 
                    show_colnames = T, 
                    fontsize_row = 9,
                    fontsize_col =7,
                    clustering_method ="ward.D", 
                    angle_col = 0)
    ggsave(plot = .ph, filename = paste0(res, cell, "/", cell, "_sensory6.png"), 
           h  = 2.2, w = 3.5, dpi = 300, unit = "in")
  ###### sensoru
    .pathway_genes <- str_split(B1a_senso$core_enrichment, "/") %>% unlist()
    .ensembl <- AnnotationDbi::select(org.Mm.eg.db, 
                                      keys = .pathway_genes, 
                                      column = c("ENSEMBL"), 
                                      keytype = "SYMBOL")
    .gene_set <- scaled_counts[scaled_counts$gene %in% .ensembl$ENSEMBL,]
    .gene_set <- .gene_set %>% left_join(.ensembl, by = c("gene" = "ENSEMBL"))
    
    rownames(.gene_set) = .gene_set$SYMBOL
    .gene_set <- .gene_set %>% dplyr::select( 1:7 ) %>% as.data.frame()
    
    names(.gene_set) <- names(.gene_set) %>% str_sub(start = 1, end = 2)
    .ph <- pheatmap(.gene_set, 
                    cluster_rows = T, 
                    cluster_cols = F, 
                    scale = "row",
                    show_rownames = T, 
                    show_colnames = T, 
                    fontsize_row = 9,
                    fontsize_col =7,
                    angle_col = 0,
                    legend = F)
    ggsave(plot = .ph, filename = paste0(res, cell, "/", cell, "_sensory3.png"), 
           w  = 3.5, h = 3, dpi = 300, unit = "in")
    }
  
  if(cell == "B1b"){
    B1b_sensor<- gse_BP[1, ]
    
    
    ### Ig production
    
    .pathway_genes <- str_split(B1b_sensor$core_enrichment, "/") %>% unlist()
    .ensembl <- AnnotationDbi::select(org.Mm.eg.db, 
                                      keys = .pathway_genes, 
                                      column = c("ENSEMBL"), 
                                      keytype = "SYMBOL")
    .gene_set <- scaled_counts[scaled_counts$gene %in% .ensembl$ENSEMBL,]
    .gene_set <- .gene_set %>% left_join(.ensembl, by = c("gene" = "ENSEMBL"))
    
    rownames(.gene_set) = .gene_set$SYMBOL.y
    
    .gene_set <- .gene_set %>% dplyr::select( 1:7 ) %>% as.data.frame()
    .gene_set <- .gene_set %>% dplyr::select(which(str_detect(names(.gene_set), "WT")), everything())
    names(.gene_set) <- names(.gene_set) %>% str_sub(start = 1, end = 2)

    .ph <- pheatmap(.gene_set, 
                    cluster_rows = F, 
                    cluster_cols = F, 
                    scale = "row",
                    show_rownames = T, 
                    show_colnames = T, 
                    fontsize_row = 9,
                    fontsize_col =9,
                    clustering_method ="ward.D", 
                    angle_col = 0)
    ggsave(plot = .ph, filename = paste0(res, cell, "/", cell, "olfactory.png"), 
           h  = 2.14, w = 4.5, dpi = 300, unit = "in")
  }
  
}



######
library(ggprism)
plt <- gse_BP2 %>% ggplot(aes(NES, reorder(Description, -p.adjust))) + 
  geom_bar(stat = "identity", aes(fill = p.adjust)) + 
  scale_fill_gradient(high="blue", low="orange") +
  scale_size(range = c(1, 4)) + xlab("enrichmentScore") + labs(fill = "p.adjust", y = "Description") + 
  guides(size = "none") + ggtitle(paste0(cell)) + theme_prism(base_size = 9)
theme_bw() %+replace% 
  theme(legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        legend.key.size = unit(1,"line"), 
        panel.grid = element_blank()) 
ggsave(h = 4, w = 6.5, plot = plt, filename = paste0(res, cell, "/", cell, "_pathwaysEnrichmentGO.png"), dpi = 300)


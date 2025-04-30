fgseaFromdds <- function(t, 
                         test, 
                         control, 
                         pathways, 
                         min,
                         max,
                         time, 
                         pathname, 
                         name) {

  
  df <- t %>% 
    results(contrast = c("group", test, control)) %>% 
    as.matrix()
  
  names <- rownames(df)
  rownames(df) <- NULL
  df <- cbind(names, df) %>% as.data.frame()
  df$baseMean <- as.numeric(df$baseMean)
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  df$lfcSE <- as.numeric(df$lfcSE)
  df$stat <- as.numeric(df$stat)
  df$pvalue <- as.numeric(df$pvalue)
  df$padj <- as.numeric(df$padj)
  
  ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                      key=df$names, 
                                      columns="SYMBOL",
                                      keytype="ENSEMBL")
  colnames(df)[1] <- "ENSEMBL"
  stat <- inner_join(df, ens2symbol) %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    dplyr::distinct() 
  
  ranks <- deframe(stat)
  
  fgseaRes <- fgsea::fgsea(pathways = pathways, 
                           stats = ranks, 
                           eps = 0,
                           minSize = min,
                           maxSize = max)
  
  #Gathering the top 10 up- and downregulated pathways, plotting and saving
  top_up <- fgseaRes[ES > 0][head(order(padj), n = 10), pathway]
  top_down <- fgseaRes[ES < 0][head(order(padj), n = 10), pathway]
  top_pathways <- c(top_up, rev(top_down))
  
  gsea_table <- plotGseaTable(pathways[top_pathways], 
                              ranks, 
                              fgseaRes, 
                              gseaParam = 0.5,
                              render = F)

  
  #Adding the leading edge genes to a separate column 
  t <- fgseaRes %>% as.data.frame() %>% mutate(leading = 0)
  
  
  
  #Creating a tidy table
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  table <- fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, log2err, ES, NES, size) %>% 
    arrange(padj)

  
  
  #Making a plot of the pathways and their normalized enrichment score, coloured by significant p values
  return(table)
}

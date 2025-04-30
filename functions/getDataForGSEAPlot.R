library(fgsea)
getDataForGSEAPlot <- function(t, test, control, pathways, pathwayname, name, time) {
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
  
  plotEnrichment(pathways[[pathwayname]],
                 ranks)
}

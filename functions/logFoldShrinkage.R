library(dplyr)
library(DESeq2)
#Annotation file available from ensembl: https://www.ensembl.org/Homo_sapiens/Info/Index
ref_annot <- read.table(here::here("seq_data/Homo_sapiens.GRCh38.95_nodupes_genid.txt"), header = TRUE, sep = "\t")

logFoldShrinkage <- function(t, c){
  r <- lfcShrink(t, coef=c, type="apeglm") #log-fold shrink
  
  p <- data.frame(Gene_id = rownames(r)) #Make a data frame
  genes <- as.data.frame(dplyr::select(ref_annot, Gene_name, Gene_id)) #A data frame with gene names and IDs
  out <- left_join(p, genes) #Join the data with the gene names
  r$geneName <- out$Gene_name 
  return(r)
  
  
  #Alternative: not using log-fold shrinkage
  #res <- results(dds, name=paste0("condition_", cond2,"_vs_",cond1))
  #res$geneName <- ref_annot[rownames(res)]$Gene_name
} 

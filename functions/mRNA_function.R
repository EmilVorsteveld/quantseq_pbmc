#Annotation file available from ensembl: https://www.ensembl.org/Homo_sapiens/Info/Index
ref_annot <- read.table(here::here("seq_data/Homo_sapiens.GRCh38.95_nodupes_genid.txt"), header = TRUE, sep = "\t")
#function to systematically exclude mtRNA, rRNA, lincRNA, snRNA,  snoRNA and miscRNA from ddsHTSeq objects
mRNA_function <- function(t){
  #Make into a DGEList
  dge <- DEFormats::as.DGEList(t)
  genes <- NULL
  
  #Get the counts per gene
  d <- dge$counts
  Gene_id <- rownames(d)
  rownames(d) <- NULL
  d <- as.data.frame(cbind(Gene_id, d))
  
  #Get the gene names and gene types
  annotate <- as.data.frame(ref_annot) %>% dplyr::select(Gene_id, Gene_type) %>% base::unique()
  
  #The genes from the data: join with the annotation file
  genes <- as.data.frame(rownames(dge$counts))
  colnames(genes) <- "Gene_id"
  genes <- left_join(genes, annotate, by = "Gene_id")
  
  #Add to the counts file
  
  d <- left_join(d, genes, by = "Gene_id")
  #d <- cbind(d, genes)
  
  #Filter all RNA types
  filtered <- d %>%
    filter(Gene_type != "Mt_rRNA" &
             Gene_type != "lincRNA" &
             Gene_type != "snRNA" &
             Gene_type != "Mt_tRNA" &
             Gene_type != "misc_RNA" &
             Gene_type != "snoRNA"
    ) %>% dplyr::select(- Gene_type)
  
  rownames(filtered) <- filtered$Gene_id
  filtered <- subset(filtered, select = -Gene_id)
  filtered <- as.matrix(filtered)
  #filtered <- filtered %>% dplyr::select(-Gene_id)

  class(filtered) <- "numeric"
  dge$counts <- filtered
  
  return(as.DESeqDataSet(dge))
}

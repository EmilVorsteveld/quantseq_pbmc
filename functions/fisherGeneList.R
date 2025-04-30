`%ni%` = Negate(`%in%`)
fisherGeneList <- function(core, other, geneset, conf) {
  non_core <- other[other %ni% core] %>% unique()
  
  t <- matrix(nrow = 2, ncol = 2)
  t[1,1] <- length(core[core %in% geneset])
  t[1,2] <- length(core[core %ni% geneset])
  
  t[2,1] <- length(non_core[non_core %in% geneset])
  t[2,2] <- length(non_core[non_core %ni% geneset])
  #print(t)
  return(fisher.test(t, conf.level = conf))
}

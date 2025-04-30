diffExpression <- function(t){
  t$group <- relevel(t$group, ref="RPMI")
  output <- DESeq(t)
  return(output)
}

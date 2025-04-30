fixConditionOrder <- function(t) {
  t$condition <- factor(t$condition, levels = c("LPS", "SAur", "PolIC", "CAlb"))
  return(t)
}

fixTimeOrder <- function(t) {
  t$`time-point` <- factor(t$`time-point`, levels = c("4hr", "24hr"))
  return(t)
}
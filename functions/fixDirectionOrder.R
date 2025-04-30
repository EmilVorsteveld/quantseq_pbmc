fixDirectionOrder <- function(t) {
  t$direction <- factor(t$direction, levels = c("up", "down"))
  return(t)
}

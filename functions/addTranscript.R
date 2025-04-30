addTranscript <- function(t) {
  names <- rownames(t) 
  t$Transcript <- names
  return(t)
}
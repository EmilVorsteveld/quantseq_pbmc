library(here)
opt.class <- c("PID", "imm.ass")

PID <- readRDS(file=here("data/PID.rda"))
imm.ass <- readRDS(file=here("data/imm.ass.rda"))

functClass <- function(file){
  file$funct.class <- ""
  
  #For each column, add the name of the class options if it is in the list
  for (i in 1:nrow(file)){
    for (class in opt.class){
      if (file$geneName[i] %in% get(class)$geneName){
        file$funct.class[i] <- c(as.character(class))
      }
    }
  }
  return(file)
}

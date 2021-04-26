dataframe_by_strain <- function (dataframe){
  x <- dataframe
  strain_abundance <- vector()
  strain_num <- vector()
  treatment <- factor()
  OD <- vector()
  for (i in 4:17){
    for (j in x[i]){
      strain_abundance<-c(strain_abundance,j)
      strain_num <- c(strain_num, rep(i-8,length(j)))
      treatment <- c(treatment,factor(x$treatment))
      OD <- c(OD,as.numeric(as.character(x$OD)))
    }
  }
  # add all to one data frame
  strain_analysis <- data.frame(strain_abundance, strain_num, treatment, OD)
  strain_analysis$strain_num <- as.factor(strain_analysis$strain_num)
  levels(strain_analysis$strain_num) <- c(1:7,11:17)
  strain_analysis$treatment <- as.factor(strain_analysis$treatment)
  levels(strain_analysis$treatment) <- levels(factor(x$treatment))
  return(strain_analysis)
}

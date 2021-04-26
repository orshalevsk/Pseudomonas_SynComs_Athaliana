#### converting to strain by strain data frame  ("strains analysis" variable)
dataframe_by_strain <- function (dataframe){
  x <- dataframe
  strain_abundance <- vector()
  strain_num <- vector()
  genotype <- factor()
  treatment <- factor()
  plate_position <- factor()
  experiment <- factor()
  tray <- factor()
  weight <- vector()
  for (i in 10:23){
    for (j in x[i]){
      strain_abundance<-c(strain_abundance,j)
      strain_num <- c(strain_num, rep(i-8,length(j)))
      genotype <- c(genotype,factor(x$genotype))
      treatment <- c(treatment,factor(x$treatment))
      plate_position <- c(plate_position,factor(x$plate_position))
      tray <- c(tray, as.character(x$tray))
      experiment <- c(experiment, as.character(x$exp))
      weight <- c(weight,as.numeric(as.character(x$weight)))
    }
  }
  # add all to one data frame
  strain_analysis <- data.frame(strain_abundance, strain_num, genotype, treatment, plate_position, tray)
  strain_analysis$strain_num <- as.factor(strain_analysis$strain_num)
  levels(strain_analysis$strain_num) <- c(1:7,11:17)
  strain_analysis$genotype <- as.factor(strain_analysis$genotype)
  levels(strain_analysis$genotype) <- levels(factor(x$genotype))
  strain_analysis$treatment <- as.factor(strain_analysis$treatment)
  levels(strain_analysis$treatment) <- levels(factor(x$treatment))
  strain_analysis$plate_position <- as.factor(strain_analysis$plate_position)
  levels(strain_analysis$plate_position) <- levels(factor(x$plate_position))
  strain_analysis$exp <- as.factor(experiment)
  strain_analysis$weight <- weight
  return(strain_analysis)
}

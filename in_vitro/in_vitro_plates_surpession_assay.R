###### summerizing the in vitro results i got in my SynCom project

library(gplots)
library(igraph)



##### starting from the plates suppression assay (from Vorholt's group)
plates_supression <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/Plates_suppression_assay_May_2019/ plates_supression_final_results_both_Exp_2020.csv")
plates_supression <- plates_supression[!is.na(plates_supression$P1),]
rownames(plates_supression) <- plates_supression$X
plates_supression <- plates_supression[,-1]
plates_supression$C3 <- as.numeric(as.character(plates_supression$C3))

heatmap.2(x = as.matrix(plates_supression),trace = "none")

library(ggplot2)
library(vegan)
library(lme4)
library(rstanarm)
library(pairwiseAdonis)
library(ellipse)
library(ggordiplots)
library(tsne)

#load Aug
setwd("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/weight+PAMP-triggered_August_2018/")
dataset_Aug <- read.csv2(file = "results/dataframe_barcodes_bacterialload_factorized_Cq_less_32_corrected_JAN_2019_single_threshold_50reads_filter.csv")
dataset_Aug$exp <- "Aug"

## load Oct
setwd("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/weight+HE_October_2018/results/")
dataset_Oct <- read.csv2("dataframe_barcodes_bacteriaload_factorized_Cq_less_32_corrected_Oct_2018_Exp_200reads_filter.csv")
levels(dataset_Oct$genotype)[levels(dataset_Oct$genotype)=="Col0_ics1"] <- "sid2-2"
dataset_Oct <- dataset_Oct[!is.na(dataset_Oct$X17),]
dataset_Oct$exp <- "Oct"
colnames_Aug <- colnames(dataset_Aug)[!colnames(dataset_Aug)=="hits_ratio"]

## load April exp
#setwd("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/April2018-July2018/Results/Miseq/")
#dataset_april <- read.csv2("dataframe_barcodes_bacterialload_factorized.csv")
#dataset_april$exp <- "April"
#dataset_april <- dataset_april[dataset_april$treatment %in% c("#1","#2","#3"),]

#colnames_all_data <- colnames(dataset_april)[colnames(dataset_april) %in% colnames(dataset_Oct)]

#merging Oct and Aug exp
full_dataset <- rbind(dataset_Oct[,colnames_Aug], dataset_Aug[,colnames_Aug])

#merging all three experiments
#full_dataset <- rbind(full_dataset_Aug_Oct[,colnames_all_data], dataset_april[,colnames_all_data])
#full_dataset <- full_dataset_Aug_Oct


#subsetting to local only
full_dataset_local <- full_dataset[full_dataset$genotype%in%c("Lu","Eyach","HE","Kus","Tu-Wal","Schl") & full_dataset$treatment%in%c("Control","#1","#2","#3"),] 
full_dataset_local$genotype <- droplevels.factor(full_dataset_local$genotype) # dropping unwanted levels (Col0, etc...)
levels(full_dataset_local$genotype)<- c("Ey15-2","HE-1","Kus3-1","Lu3-30","Schl-7", "Tu-Wal-2") #renaming host genotypes to original 1001 names
full_dataset_local$weight[full_dataset_local$weight<1] <- full_dataset_local$weight[full_dataset_local$weight<1]*1000


#removing outliers function - more than 2.5 sds from mean of chosen factor (weight / abundance / etc...) of genotype*treatment*experiment
outliers_removal <- function(dataset, factor){
  outliers_list <- vector()
  for (exp in unique(dataset$exp)){
    for (gen in unique(dataset$genotype)){
      for (treat in unique(dataset$treatment)){
        gen_by_treat_by_exp <- dataset$treatment==treat & dataset$genotype==gen & dataset$exp==exp
        mean_gen_by_treat_by_exp <- mean(dataset[gen_by_treat_by_exp,factor])
        sd_2.5 <-2.5*sd(dataset[gen_by_treat_by_exp,factor])
        outliers <- dataset[gen_by_treat_by_exp,factor] > (mean_gen_by_treat_by_exp+sd_2.5) | dataset[gen_by_treat_by_exp,factor] <(mean_gen_by_treat_by_exp-sd_2.5)
        if (sum(outliers)>0){
          outliers_list <- c(outliers_list,row.names(dataset[gen_by_treat_by_exp,][outliers,]))
        }
      }
    }
  }
  return(outliers_list)
}

#removing weight outliers
outliers_list <- outliers_removal(dataset = full_dataset_local, factor = "weight") # removing outliers - weight
length(outliers_list)/length(full_dataset_local$weight) # Sanity check: the % of outliers removed. about 1%, as expected (2.5 std are 99% of the normal distrubiton)
full_dataset_local_no_out <- full_dataset_local[-which(row.names(full_dataset_local) %in% outliers_list),]

hist(log(full_dataset_local_no_out$bacterial_load), breaks = 100) # Analysis - how does the abundance distributed? Normally. At least the log(abudnance)



#removing bacterial load outliers (log bacterial load). No outliers were found after all (list is empty)
full_dataset_local_no_out$log_bacterial_load <- log(full_dataset_local_no_out$bacterial_load)
outliers_list <- outliers_removal(dataset = full_dataset_local_no_out, factor = "log_bacterial_load") # removing outliers - log(bacterial_load). No outliers were found
if (length(outliers_list)>0){
  full_dataset_local_no_out <- full_dataset_local_no_out[-which(row.names(full_dataset_local_no_out) %in% outliers_list),]
}


#### adding patho and com under the mixed com
full_dataset_local_no_out$com3_patho[full_dataset_local_no_out$treatment=="#3"] <- rowSums(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",paste("X",c(1:7),sep = "")])
full_dataset_local_no_out$com3_commensal[full_dataset_local_no_out$treatment=="#3"] <- rowSums(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",paste("X",c(11:17),sep = "")])

'
#removing zeros. These cannot be log() transformed
for (i in colnames(full_dataset_local_no_out)){
  zero_row <- which(full_dataset_local_no_out[,i]==0)
  if (length(zero_row)>0){
    full_dataset_local_no_out <- full_dataset_local_no_out[-zero_row,]
  }
}
'

#####Fig XXXX - disimilarity between genotypes, i.e. the genotype effect - treatment-wise. Both PERMANOVA and ANOSIM (both are Bary-Curtis based)
#####Seems like a supplementary, and i should only report the % of similarity each 

#This part is to section by community. I will only include data from a certain community when analysing it (e.g. #1 = only columns X1-X7)
strains_col_names <- c("X1","X2","X3","X4","X5","X6","X7","X11","X12","X13","X14","X15","X16","X17")
pathogens_col_names <- c("X1","X2","X3","X4","X5","X6","X7")
commensals_col_names <- c("X11","X12","X13","X14","X15","X16","X17")


####Bray-Curtis based PERMANOVA test (adonis)
perm <- function(genotypes, formula) {
  samples_list_full_design_adonis <- genotypes[genotypes$treatment%in%c("#1","#2","#3"),]
  samples_list_full_design_adonis[samples_list_full_design_adonis$treatment=="#1",commensals_col_names] <- 0
  samples_list_full_design_adonis[samples_list_full_design_adonis$treatment=="#2",pathogens_col_names] <- 0
  permanova <- adonis(formula = formula ,data = samples_list_full_design_adonis,permutations = 2000, method = "bray")
  return(permanova)
}

# each experiment seperately
formula_treat1_perm <- as.formula("samples_list_full_design_adonis[pathogens_col_names] ~ genotype") #specify a formula. This is treatment-wise
formula_treat2_perm <- as.formula("samples_list_full_design_adonis[commensals_col_names] ~ genotype") #specify a formula. This is treatment-wise
formula_treat3_perm <- as.formula("samples_list_full_design_adonis[strains_col_names] ~ genotype") #specify a formula. This is treatment-wise

permanova_treat1_oct <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct",], formula = formula_treat1_perm) # specify the subset of treatments
permanova_treat2_oct <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",], formula = formula_treat2_perm) # specify the subset of treatments
permanova_treat3_oct <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",], formula = formula_treat3_perm) # specify the subset of treatments

permanova_treat1_Aug <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Aug",], formula = formula_treat1_perm) # specify the subset of treatments
permanova_treat2_Aug <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Aug",], formula = formula_treat2_perm) # specify the subset of treatments
permanova_treat3_Aug <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug",], formula = formula_treat3_perm) # specify the subset of treatments


permanova_treat1_April <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="April",], formula = formula_treat1_perm) # specify the subset of treatments
permanova_treat2_April <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="April",], formula = formula_treat2_perm) # specify the subset of treatments
permanova_treat3_April <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="April",], formula = formula_treat3_perm) # specify the subset of treatments



permanova_treat1_df_oct <- as.data.frame(permanova_treat1_oct$aov.tab[1:6])
permanova_treat1_df_oct$treatment <- "#1"
permanova_treat1_df_oct$exp <- "oct"
permanova_treat1_df_oct$factor <- rownames(permanova_treat1_df_oct)

permanova_treat2_df_oct <- as.data.frame(permanova_treat2_oct$aov.tab[1:6])
permanova_treat2_df_oct$treatment <- "#2"
permanova_treat2_df_oct$exp <- "oct"
permanova_treat2_df_oct$factor <- rownames(permanova_treat2_df_oct)

permanova_treat3_df_oct <- as.data.frame(permanova_treat3_oct$aov.tab[1:6])
permanova_treat3_df_oct$treatment <- "#3"
permanova_treat3_df_oct$exp <- "oct"
permanova_treat3_df_oct$factor <- rownames(permanova_treat3_df_oct)

permanova_treat1_df_Aug <- as.data.frame(permanova_treat1_Aug$aov.tab[1:6])
permanova_treat1_df_Aug$treatment <- "#1"
permanova_treat1_df_Aug$exp <- "Aug"
permanova_treat1_df_Aug$factor <- rownames(permanova_treat1_df_Aug)

permanova_treat2_df_Aug <- as.data.frame(permanova_treat2_Aug$aov.tab[1:6])
permanova_treat2_df_Aug$treatment <- "#2"
permanova_treat2_df_Aug$exp <- "Aug"
permanova_treat2_df_Aug$factor <- rownames(permanova_treat2_df_Aug)

permanova_treat3_df_Aug <- as.data.frame(permanova_treat3_Aug$aov.tab[1:6])
permanova_treat3_df_Aug$treatment <- "#3"
permanova_treat3_df_Aug$exp <- "Aug"
permanova_treat3_df_Aug$factor <- rownames(permanova_treat3_df_Aug)


permanova_by_treat_both_exp_seperatley <- rbind(permanova_treat1_df_oct,permanova_treat2_df_oct,permanova_treat3_df_oct,
                                                permanova_treat1_df_Aug,permanova_treat2_df_Aug,permanova_treat3_df_Aug)

# all experiments
formula_treat1_perm <- as.formula("samples_list_full_design_adonis[pathogens_col_names] ~ genotype*exp") #specify a formula. This is treatment-wise
formula_treat2_perm <- as.formula("samples_list_full_design_adonis[commensals_col_names] ~ genotype*exp") #specify a formula. This is treatment-wise
formula_treat3_perm <- as.formula("samples_list_full_design_adonis[strains_col_names] ~ genotype*exp") #specify a formula. This is treatment-wise

permanova_treat1 <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",], formula = formula_treat1_perm) # specify the subset of treatments
permanova_treat2 <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2",], formula = formula_treat2_perm) # specify the subset of treatments
permanova_treat3 <- perm(genotypes = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",], formula = formula_treat3_perm) # specify the subset of treatments


permanova_treat1_df <- as.data.frame(permanova_treat1$aov.tab[1:6])
permanova_treat1_df$treatment <- "#1"
permanova_treat1_df$factor <- rownames(permanova_treat1_df)

permanova_treat2_df <- as.data.frame(permanova_treat2$aov.tab[1:6])
permanova_treat2_df$treatment <- "#2"
permanova_treat2_df$factor <- rownames(permanova_treat2_df)

permanova_treat3_df <- as.data.frame(permanova_treat3$aov.tab[1:6])
permanova_treat3_df$treatment <- "#3"
permanova_treat3_df$factor <- rownames(permanova_treat3_df)

permanova_by_treat_both_exp <- rbind(permanova_treat1_df,permanova_treat2_df,permanova_treat3_df)
permanova_by_treat_both_exp$exp <- "all"

permanova_by_treat <- rbind(permanova_by_treat_both_exp_seperatley,permanova_by_treat_both_exp)
permanova_by_treat_genotype <- permanova_by_treat[permanova_by_treat$factor %in% c("genotype"),]

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/composition_PERMANOVA_genotype_in_treat_within_between_exp.pdf")
ggplot(aes(x = treatment, y = R2, color = exp), data = permanova_by_treat_genotype) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw()
#! dev.off()

#! write.csv2(x = permanova_by_treat, file = "~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/composition_PERMANOVA_genotype_in_treat_within_between_exp.csv")


####Bray-Curtis based ANOSIM. TODO: choose experiment or find a way to analyze both.
anosim_treat1_oct <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct",strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct"], permutations = 2000, distance = "bray")
anosim_treat2_oct <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct"], permutations = 2000, distance = "bray")
anosim_treat3_oct <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct"], permutations = 2000, distance = "bray")

anosim_treat1_Aug <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Aug",strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Aug"], permutations = 2000, distance = "bray")
anosim_treat2_Aug <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Aug",strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Aug"], permutations = 2000, distance = "bray")
anosim_treat3_Aug <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug",strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug"], permutations = 2000, distance = "bray")

anosim_treat1 <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" ,pathogens_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#1"], permutations = 2000, distance = "bray")
anosim_treat2 <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" ,commensals_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#2"], permutations = 2000, distance = "bray")
anosim_treat3 <- anosim( x = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" ,strains_col_names], grouping = full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#3"], permutations = 2000, distance = "bray")

anosim_all <- anosim( x = full_dataset_local_no_out[ ,strains_col_names], grouping = full_dataset_local_no_out$treatment, permutations = 200, distance = "bray")



anosim_summary <- data.frame("ANOSIM_statistic_R"=NA, "signif"=NA,"treatment"=NA,"exp"=NA)
anosim_summary <- rbind(anosim_summary, c(anosim_treat1_Aug$statistic,anosim_treat1_Aug$signif,"#1","Aug"), c(anosim_treat2_Aug$statistic,anosim_treat2_Aug$signif,"#2","Aug"), c(anosim_treat3_Aug$statistic,anosim_treat3_Aug$signif,"#3","Aug"),
      c(anosim_treat1_oct$statistic,anosim_treat1_oct$signif,"#1","Oct"),c(anosim_treat2_oct$statistic,anosim_treat2_oct$signif,"#2","Oct"),c(anosim_treat3_oct$statistic,anosim_treat3_oct$signif,"#3","Oct"))
anosim_summary <- anosim_summary[!is.na(anosim_summary$exp),]

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/composition_anosim_genotype_within_exp.pdf")
ggplot(aes(x = treatment, y = ANOSIM_statistic_R, color = exp), data = anosim_summary) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw()
#! dev.off()

#! write.csv2(x = anosim_summary, file = "~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/composition_anosim_genotype_within_exp.csv")



## pair-wise PERMANOVA - treatment-wise
pairwise_adonis_summary <- data.frame("treatment"=NA,"comparison"=NA,"Pval"=NA,"R2"=NA)

pariwaise_adonis_wild_subset_gen_com1 <- pairwise.adonis2(dist(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct",strains_col_names])~genotype,
                                                     data=full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct",])

pariwaise_adonis_wild_subset_gen_com2 <- pairwise.adonis2(dist(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",strains_col_names])~genotype,
                                                          data=full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",])

pariwaise_adonis_wild_subset_gen_com3 <- pairwise.adonis2(dist(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",strains_col_names])~genotype,
                                                          data=full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",])

for (comparison_index in 2:length(names(pariwaise_adonis_wild_subset_gen_com1))){
  R2_current_comparison <- pariwaise_adonis_wild_subset_gen_com1[names(pariwaise_adonis_wild_subset_gen_com1)[comparison_index]][1][[1]][[5]][1]                                                                        
  Pval_cccurrent_comparison <- pariwaise_adonis_wild_subset_gen_com1[names(pariwaise_adonis_wild_subset_gen_com1)[comparison_index]][1][[1]][[6]][1]                                                                        
  current_comparison <- names(pariwaise_adonis_wild_subset_gen_com1)[comparison_index]
  treatment <- "#1"
  pairwise_adonis_summary <- rbind(pairwise_adonis_summary,c(treatment, current_comparison,Pval_cccurrent_comparison, R2_current_comparison))
}

for (comparison_index in 2:length(names(pariwaise_adonis_wild_subset_gen_com2))){
  R2_current_comparison <- pariwaise_adonis_wild_subset_gen_com2[names(pariwaise_adonis_wild_subset_gen_com2)[comparison_index]][1][[1]][[5]][1]                                                                        
  Pval_cccurrent_comparison <- pariwaise_adonis_wild_subset_gen_com2[names(pariwaise_adonis_wild_subset_gen_com2)[comparison_index]][1][[1]][[6]][1]                                                                        
  current_comparison <- names(pariwaise_adonis_wild_subset_gen_com2)[comparison_index]
  treatment <- "#2"
  pairwise_adonis_summary <- rbind(pairwise_adonis_summary,c(treatment, current_comparison,Pval_cccurrent_comparison, R2_current_comparison))
}

for (comparison_index in 2:length(names(pariwaise_adonis_wild_subset_gen_com3))){
  R2_current_comparison <- pariwaise_adonis_wild_subset_gen_com3[names(pariwaise_adonis_wild_subset_gen_com3)[comparison_index]][1][[1]][[5]][1]                                                                        
  Pval_cccurrent_comparison <- pariwaise_adonis_wild_subset_gen_com3[names(pariwaise_adonis_wild_subset_gen_com3)[comparison_index]][1][[1]][[6]][1]                                                                        
  current_comparison <- names(pariwaise_adonis_wild_subset_gen_com3)[comparison_index]
  treatment <- "#3"
  pairwise_adonis_summary <- rbind(pairwise_adonis_summary,c(treatment, current_comparison,Pval_cccurrent_comparison, R2_current_comparison))
}

significant_pairwise_adonis_summary <- pairwise_adonis_summary[pairwise_adonis_summary$Pval < 0.05, ]
significant_pairwise_adonis_summary <- significant_pairwise_adonis_summary[!is.na(significant_pairwise_adonis_summary$treatment),]

for (index in 1:length(significant_pairwise_adonis_summary$treatment)){
  significant_pairwise_adonis_summary$gen1[index] <- strsplit(x = significant_pairwise_adonis_summary$comparison, split = "_vs_")[[index]][1]
  significant_pairwise_adonis_summary$gen2[index] <- strsplit(x = significant_pairwise_adonis_summary$comparison, split = "_vs_")[[index]][2]
}

#! write.csv2(x = significant_pairwise_adonis_summary,file = "~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/pairwise_gen_by_treat_anosim.csv")


#### NMDS analysis and plotting

# NMDS of both experiments, all treatments. This is to examine differences between experiments
NMDS_all_com=metaMDS(comm = full_dataset_local_no_out[,strains_col_names], distance = "bray", k=3,trymax = 40 ,autotransform = T,plot = T) # k = the number of final reduced dimensions


# NMDS only of october experiment (less stochstic, abundance-wise). All treatments. This is to examine differences between treatments.
NMDS_all_com_oct=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct" ,strains_col_names], distance = "bray", k=3,trymax = 200 ,autotransform = T,plot = T) # k = the number of final reduced dimensions


# NMDS only of october experiment (less stochstic, abundance-wise). Pathogens and commensals only. This is to examine differences between the same strains in different treatments.
NMDS_all_com_oct_patho =metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#1","#3") & full_dataset_local_no_out$exp=="Oct" ,pathogens_col_names], distance = "bray", k=2,trymax = 300 ,autotransform = T,plot = T) # k = the number of final reduced dimensions
NMDS_all_com_oct_commensal =metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#2","#3") & full_dataset_local_no_out$exp=="Oct" ,commensals_col_names], distance = "bray", k=2,trymax = 300 ,autotransform = T,plot = T) # k = the number of final reduced dimensions


# NMDS only of october experiment (less stochstic, abundance-wise). Treatments-wise. This is to examine differences between genotypes, within a treatment.
NMDS_com1_oct=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct" ,pathogens_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T,plot = T) # k = the number of final reduced dimensions
NMDS_com2_oct=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",commensals_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T, plot = T) # k = the number of final reduced dimensions
NMDS_com3_oct=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",strains_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T, plot = T) # k = the number of final reduced dimensions

#now getting the centroids of each genotype in each dataset
x <-full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct" ,pathogens_col_names]
envfit_com1_oct <- envfit(x ~ full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct"] , full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct",])

x <-full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct" ,commensals_col_names]
envfit_com2_oct <- envfit(x ~ full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct"] , full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",])

x <-full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct" ,strains_col_names]
envfit_com3_oct <- envfit(x ~ full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct"] , full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",])


# NMDS only of August experiment (less stochstic, abundance-wise). Treatments-wise. This is to examine differences between genotypes, within a treatment.
NMDS_com1_Aug=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Aug" ,pathogens_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T,plot = T) # k = the number of final reduced dimensions
NMDS_com2_Aug=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Aug",commensals_col_names], distance = "bray", k=4,trymax = 200, autotransform = T, plot = T) # k = the number of final reduced dimensions
NMDS_com3_Aug=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug",strains_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T, plot = T) # k = the number of final reduced dimensions



# all communities
NMDS_com1=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" ,pathogens_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T,plot = T) # k = the number of final reduced dimensions
NMDS_com2=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" ,commensals_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T, plot = T) # k = the number of final reduced dimensions
NMDS_com3=metaMDS(comm = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" ,strains_col_names], distance = "bray", k=4,trymax = 200 ,autotransform = T, plot = T) # k = the number of final reduced dimensions


weight_by_NMDS <- as.data.frame(cbind(NMDS_com3_oct$points[,4],full_dataset_local_no_out$weight[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#3"]))
weight_by_NMDS <- as.data.frame(cbind(NMDS_com1_oct$points[,1],full_dataset_local_no_out$weight[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#1"]))
weight_by_NMDS <- as.data.frame(cbind(NMDS_com3_oct$points[,4],full_dataset_local_no_out$weight[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#3"]))
weight_by_NMDS <- as.data.frame(cbind(NMDS_com3_oct$points[,4],full_dataset_local_no_out$weight[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#3"]))

cor.test(weight_by_NMDS$V1,weight_by_NMDS$V2)

#### correlation of load and NMDS1
load_by_NMDS_com1 <- as.data.frame(cbind(NMDS_com1_oct$points[,1],full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#1"]))
load_by_NMDS_com2 <- as.data.frame(cbind(NMDS_com2_oct$points[,1],full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#2"]))
load_by_NMDS_com3 <- as.data.frame(cbind(NMDS_com3_oct$points[,1],full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment=="#3"]))

cor.test(load_by_NMDS_com1$V1,load_by_NMDS_com1$V2)
cor.test(load_by_NMDS_com2$V1,load_by_NMDS_com2$V2)
cor.test(load_by_NMDS_com3$V1,load_by_NMDS_com3$V2)



ggplot(data = load_by_NMDS_com1, aes(x = V1, y=V2))+
  geom_point() +
  geom_smooth(method = "lm")


#! stressplot(NMDS_all_com_oct) # check for the stress - some kind of goodness of fit/model. Read more about it if i want to understand
### plotting both exp, with confidence ellipse.
NMDS = data.frame(MDS1 = NMDS_all_com$points[,1], MDS2 = 
                    NMDS_all_com$points[,2],MDS3=NMDS_all_com$points[,3],group=full_dataset_local_no_out$treatment, exp=full_dataset_local_no_out$exp)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$group),mean)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in unique(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))}

cent <- aggregate(cbind(NMDS$MDS1, NMDS$MDS2) ~ treatment, data = full_dataset_local_no_out, FUN = mean)

NMDS_all_com_plot_ellipse<-ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(aes(color = group, shape = exp), stroke =0.2,size=2.5,alpha=0.7) +
  geom_point(data = cent,aes(x=V1,y=V2), size = 10, shape="*") +  
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")
  )+
  theme(legend.text = element_text(size = 14)
  )+
  theme(legend.title = element_text(size = 14)
  )+
  geom_polygon(data=df_ell, aes(x=MDS1, y=MDS2, color = group, fill=group), linetype="dashed", size =1, alpha=0.2) +
  scale_colour_manual(name = "Treatment", values = c("#d95f02", "#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  geom_text(x=-1, y=-0.8, label="stress=0.08017") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_both_exp_with_centroids.pdf", useDingbats=F)
print(NMDS_all_com_plot_ellipse)
#dev.off()


##now plotting heatmap for both exp, all samples, by treatment (using Bray Curtis distances)
full_dataset_local_no_out$color_by_treat[full_dataset_local_no_out$treatment=="#1"] <- "#d95f02"
full_dataset_local_no_out$color_by_treat[full_dataset_local_no_out$treatment=="#2"] <- "#1b9e77"
full_dataset_local_no_out$color_by_treat[full_dataset_local_no_out$treatment=="#3"] <- "#7570b3"

full_dataset_local_no_out$color_by_exp[full_dataset_local_no_out$exp=="Aug"] <- "grey30"
full_dataset_local_no_out$color_by_exp[full_dataset_local_no_out$exp=="Oct"] <- "grey80"

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/heatmap_all_treat_both_exp.pdf", useDingbats = F)  
#png("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/heatmap_all_treat_both_exp.png",width = 5*600,height = 5*600,res = 600, pointsize = 8) # 5 x 300 pixels, 300 pixels per inch
gplots::heatmap.2(as.matrix(vegdist(full_dataset_local_no_out[,strains_col_names], method="bray")),trace = 'none',
                  RowSideColors = full_dataset_local_no_out$color_by_treat, ColSideColors = full_dataset_local_no_out$color_by_exp,labRow = FALSE, labCol = FALSE)
#dev.off()


##now plotting heatmap per treatment, by genotype (using Bray Curtis distances)
full_data_com1<- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",]
full_data_com1$col_by_gen[full_data_com1$genotype=="Ey15-2"] <- "#d73027"
full_data_com1$col_by_gen[full_data_com1$genotype=="HE-1"] <- "#fc8d59"
full_data_com1$col_by_gen[full_data_com1$genotype=="Kus3-1"] <- "#fee090"
full_data_com1$col_by_gen[full_data_com1$genotype=="Lu3-30"] <- "#e0f3f8"
full_data_com1$col_by_gen[full_data_com1$genotype=="Schl-7"] <- "#91bfdb"
full_data_com1$col_by_gen[full_data_com1$genotype=="Tu-Wal-2"] <- "#4575b4"

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/heatmap_all_treat_.pdf", useDingbats = F)  
#png("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/heatmap_all_treat_.png",width = 5*600,height = 5*600,res = 600, pointsize = 8) # 5 x 300 pixels, 300 pixels per inch
gplots::heatmap.2(as.matrix(vegdist(full_data_com1[,strains_col_names], method="bray")),trace = 'none',
                  RowSideColors = full_data_com1$col_by_gen,labRow = FALSE, labCol = FALSE)
#dev.off()




# plotting the difference between experiments (colors)
ggplot(aes(x=NMDS_all_com$points[,1], y=NMDS_all_com$points[,2], color = exp, shape = treatment), 
       data = full_dataset_local_no_out)+
  geom_point(size = 3) +
  theme_bw()

# plotting the difference between treatments (colors) in both experiments (shape)
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_both_exp.pdf", useDingbats=F)
ggplot(aes(x=NMDS_com2$points[,3], y=NMDS_com2$points[,4], color = genotype, shape = exp), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2",])+
  geom_point(size = 3) +
  theme_bw()
#dev.off()



# plotting the difference between treatments - only in oct
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_with_load.pdf" , useDingbats=F)
ggplot(aes(x=NMDS_all_com_oct$points[,1], y=NMDS_all_com_oct$points[,2], color = log(bacterial_load), shape=treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",])+
  scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()
#dev.off()


# plotting the difference between treatments - only in oct
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_with_load_dim2_3.pdf" , useDingbats=F)
ggplot(aes(x=NMDS_all_com_oct$points[,2], y=NMDS_all_com_oct$points[,3], color = log(bacterial_load), shape=treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",])+
  scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()
#dev.off()



# plotting correlation betweeb NMDS dimension 2 to weight
ggplot(aes(x=NMDS_all_com_oct$points[,2], y=weight), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",])+
  #scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()

## plotting the difference between treatments - only in oct, and subsetting to OTU5 and non-OTU5
# pathogens
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_pathogens_only.pdf")
ggplot(aes(x=NMDS_all_com_oct_patho$points[,1], y=NMDS_all_com_oct_patho$points[,2], color = treatment, shape=treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment%in%c("#1","#3"),])+
  #scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_pathogens_only_with_load.pdf")
ggplot(aes(x=NMDS_all_com_oct_patho$points[,1], y=NMDS_all_com_oct_patho$points[,2], color = log(bacterial_load), shape=treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment%in%c("#1","#3"),])+
  scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()
#! dev.off()


# commensals
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_commensals_only.pdf")
ggplot(aes(x=NMDS_all_com_oct_commensal$points[,1], y=NMDS_all_com_oct_commensal$points[,2], color = treatment, shape=treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment%in%c("#2","#3"),])+
  #scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_commensals_only_with_load.pdf")
ggplot(aes(x=NMDS_all_com_oct_commensal$points[,1], y=NMDS_all_com_oct_commensal$points[,2], color = log(bacterial_load), shape=treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct" & full_dataset_local_no_out$treatment%in%c("#2","#3"),])+
  scale_color_gradient(low = "blue",high = "red") +
  geom_point(size = 3) +
  theme_bw()
#! dev.off()




#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_all_treat_oct_exp_by_treat.pdf", useDingbats=F)
ggplot(aes(x=NMDS_all_com_oct$points[,1], y=NMDS_all_com_oct$points[,2], color = treatment), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",])+
  geom_point(size = 3) +
  theme_bw()
#dev.off()


min_load <- min(log(full_dataset_local_no_out$bacterial_load))
max_load<- max(log(full_dataset_local_no_out$bacterial_load))

NMDS_com3_site<- as.data.frame(scores(NMDS_com3_oct, "species"))
NMDS_com3_site$strain <- rownames(NMDS_com3_site)  

# plotting the difference between genotypes, within treatments - only in oct
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_MixedCom_oct_exp_by_gen_load_species.pdf", useDingbats = F)
ggplot()+
  geom_point(size = 3,aes(x=NMDS_com3_oct$points[,1], y=NMDS_com3_oct$points[,2], color = log(bacterial_load), shape=genotype), # try also color="log(bacterial_load)", shape="genotype"
             data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",]) +
  scale_colour_gradient(low = "blue", high = "red") +
  geom_text(data =NMDS_com3_site, aes(x=NMDS1,y=NMDS2,label=strain),size=6,vjust=0, alpha=1) +
  annotate(geom="text", x=-0, y=-0.7, label=as.character(NMDS_com3_oct$stress),
           color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #coord_cartesian(xlim = c(-1.5,2), ylim = c(-1,1)) +
  #stat_ellipse(mapping = NULL, data = NULL, geom = "path",
  #             position = "identity", type = "t", level = 0.95,
  #             segments = 51, na.rm = FALSE, show.legend = NA,
  #             inherit.aes = TRUE) +
#! dev.off()


NMDS_com2_site<- as.data.frame(scores(NMDS_com2_oct, "species"))
NMDS_com2_site$strain <- rownames(NMDS_com2_site)  

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_CommenCom_oct_exp_by_gen_load_species.pdf", useDingbats = F)
ggplot()+
  geom_point(size = 3,aes(x=NMDS_com2_oct$points[,1], y=NMDS_com2_oct$points[,2], color = log(bacterial_load), shape=genotype),  # try also color="log(bacterial_load)", shape="genotype"
             data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",]) +
  geom_text(data =NMDS_com2_site, aes(x=NMDS1,y=NMDS2,label=strain),size=6,vjust=0,alpha=1) +
  scale_colour_gradient(low = "blue", high = "red") +
  annotate(geom="text", x=0.6, y=0.6, label=as.character(NMDS_com2_oct$stress),
           color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #stat_ellipse(mapping = NULL, data = NULL, geom = "path",
  #             position = "identity", type = "norm", level = 0.95,
  #             segments = 51, na.rm = FALSE, show.legend = NA,
  #             inherit.aes = TRUE) +
  #coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-0.6,0.6)) +
#! dev.off()




NMDS_com1_site<- as.data.frame(scores(NMDS_com1_oct, "species"))
NMDS_com1_site$strain <- rownames(NMDS_com1_site)  

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_PathoCom_oct_exp_by_gen_load_species.pdf", useDingbats = F)
ggplot()+
  geom_point(size = 3,aes(x=NMDS_com1_oct$points[,1], y=NMDS_com1_oct$points[,2], color=log(bacterial_load), shape=genotype),  # try also color="log(bacterial_load)", shape="genotype"
             data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Oct",]) +
  geom_text(data =NMDS_com1_site, aes(x=NMDS1,y=NMDS2,label=strain),size=6,vjust=0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  annotate(geom="text", x=0.6, y=0.6, label=as.character(NMDS_com1_oct$stress),
           color="black") +
  scale_colour_gradient(low = "blue", high = "red") 
  #coord_cartesian(xlim = c(-2,2.5), ylim = c(-1.2,1.2)) +
  #stat_ellipse(mapping = NULL, data = NULL, geom = "path",
  #             position = "identity", type = "t", level = 0.95,
  #             segments = 51, na.rm = FALSE, show.legend = NA,
  #             inherit.aes = TRUE) +
#! dev.off()




# plotting the difference between genotypes, within treatments - only in oct
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_MixedCom_oct_exp_by_gen.pdf", useDingbats = F)
ggplot(aes(x=NMDS_com3_oct$points[,1], y=NMDS_com3_oct$points[,2], color = genotype), # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",])+
  geom_point(size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#coord_cartesian(xlim = c(-1.5,2), ylim = c(-1,1)) +
#stat_ellipse(mapping = NULL, data = NULL, geom = "path",
#             position = "identity", type = "t", level = 0.95,
#             segments = 51, na.rm = FALSE, show.legend = NA,
#             inherit.aes = TRUE) +
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_CommenCom_oct_exp_by_gen.pdf", useDingbats = F)
ggplot(aes(x=NMDS_com2_oct$points[,1], y=NMDS_com2_oct$points[,2], color = genotype),  # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Oct",])+
  geom_point(size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#stat_ellipse(mapping = NULL, data = NULL, geom = "path",
#             position = "identity", type = "norm", level = 0.95,
#             segments = 51, na.rm = FALSE, show.legend = NA,
#             inherit.aes = TRUE) +
#coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-0.6,0.6)) +
#! dev.off()


#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_PathoCom_oct_exp_by_gen.pdf", useDingbats = F)
ggplot()+
  geom_point(data = NMDS_com1_results,size = 3,aes(x=NMDS1, y=NMDS2, color=genotype)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#coord_cartesian(xlim = c(-2,2.5), ylim = c(-1.2,1.2)) +
#stat_ellipse(mapping = NULL, data = NULL, geom = "path",
#             position = "identity", type = "t", level = 0.95,
#             segments = 51, na.rm = FALSE, show.legend = NA,
#             inherit.aes = TRUE) +
#! dev.off()


  







# plotting the difference between genotypes, within treatments - only in Aug
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_MixedCom_Aug_exp_by_gen_load.pdf")
ggplot(aes(x=NMDS_com3_Aug$points[,1], y=NMDS_com3_Aug$points[,2], color = log(bacterial_load), shape=genotype), # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug",])+
  geom_point(size = 3) +
  scale_colour_gradient(low = "blue", high = "red") +
  #coord_cartesian(xlim = c(-1.5,2), ylim = c(-1,1)) +
  #stat_ellipse(mapping = NULL, data = NULL, geom = "path",
  #             position = "identity", type = "t", level = 0.95,
  #             segments = 51, na.rm = FALSE, show.legend = NA,
  #             inherit.aes = TRUE) +
  theme_bw()
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_CommenCom_Aug_exp_by_gen_load.pdf")
ggplot(aes(x=NMDS_com2_Aug$points[,1], y=NMDS_com2_Aug$points[,2], color = log(bacterial_load), shape=genotype),  # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$exp=="Aug",])+
  geom_point(size = 3) +
  scale_colour_gradient(low = "blue", high = "red") +
  #stat_ellipse(mapping = NULL, data = NULL, geom = "path",
  #             position = "identity", type = "norm", level = 0.95,
  #             segments = 51, na.rm = FALSE, show.legend = NA,
  #             inherit.aes = TRUE) +
  #coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-0.6,0.6)) +
  theme_bw()
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_PathoCom_Aug_exp_by_gen_load.pdf")
ggplot(aes(x=NMDS_com1_Aug$points[,1], y=NMDS_com1_Aug$points[,2], color = log(bacterial_load), shape=genotype),  # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1" & full_dataset_local_no_out$exp=="Aug",])+
  geom_point(size = 3) +
  scale_colour_gradient(low = "blue", high = "red") +
  #coord_cartesian(xlim = c(-2,2.5), ylim = c(-1.2,1.2)) +
  #stat_ellipse(mapping = NULL, data = NULL, geom = "path",
  #             position = "identity", type = "t", level = 0.95,
  #             segments = 51, na.rm = FALSE, show.legend = NA,
  #             inherit.aes = TRUE) +
  theme_bw()
#! dev.off()


min(significant_pairwise_adonis_summary$R2)
max(significant_pairwise_adonis_summary$R2)
hist(as.numeric(significant_pairwise_adonis_summary$R2))
median(as.numeric(significant_pairwise_adonis_summary$R2))
## now the same NMDS with relative abundance
full_dataset_local_no_out_relative <- full_dataset_local_no_out
full_dataset_local_no_out_relative[,strains_col_names] <- full_dataset_local_no_out_relative[,strains_col_names]/rowSums(full_dataset_local_no_out_relative[,strains_col_names])

# NMDS only of october experiment (less stochstic, abundance-wise). Treatments-wise. This is to examine differences between genotypes, within a treatment.
NMDS_com1_oct_relative=metaMDS(comm = full_dataset_local_no_out_relative[full_dataset_local_no_out_relative$treatment=="#1" & full_dataset_local_no_out_relative$exp=="Oct" ,pathogens_col_names], distance = "bray", k=2,trymax = 2000 ,autotransform = T,plot = T) # k = the number of final reduced dimensions
NMDS_com2_oct_relative=metaMDS(comm = full_dataset_local_no_out_relative[full_dataset_local_no_out_relative$treatment=="#2" & full_dataset_local_no_out_relative$exp=="Oct",commensals_col_names], distance = "bray", k=2,trymax = 2000 ,autotransform = T, plot = T) # k = the number of final reduced dimensions
NMDS_com3_oct_relative=metaMDS(comm = full_dataset_local_no_out_relative[full_dataset_local_no_out_relative$treatment=="#3" & full_dataset_local_no_out_relative$exp=="Oct",strains_col_names], distance = "bray", k=2,trymax = 2000 ,autotransform = T, plot = T) # k = the number of final reduced dimensions

# plotting the difference between genotypes, within treatments - only in oct
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_MixedCom_oct_exp_by_gen_load.pdf")
ggplot(aes(x=NMDS_com1_oct_relative$points[,1], y=NMDS_com1_oct_relative$points[,2], color = genotype, shape=genotype), # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out_relative[full_dataset_local_no_out_relative$treatment=="#1" & full_dataset_local_no_out_relative$exp=="Oct",])+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#! dev.off()


#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_MixedCom_oct_exp_by_gen_load.pdf")
ggplot(aes(x=NMDS_com2_oct_relative$points[,1], y=NMDS_com2_oct_relative$points[,2], color = genotype, shape=genotype), # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out_relative[full_dataset_local_no_out_relative$treatment=="#2" & full_dataset_local_no_out_relative$exp=="Oct",])+
  geom_point(size = 3) +
  theme_bw()
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/composition_by_genotype/NMDS_MixedCom_oct_exp_by_gen_load.pdf")
ggplot(aes(x=NMDS_com3_oct_relative$points[,1], y=NMDS_com3_oct_relative$points[,2], color = genotype, shape=genotype), # try also color="log(bacterial_load)", shape="genotype"
       data = full_dataset_local_no_out_relative[full_dataset_local_no_out_relative$treatment=="#3" & full_dataset_local_no_out_relative$exp=="Oct",])+
  geom_point(size = 3) +
  theme_bw()
#! dev.off()




#### tsne analysis and plots as well
tsne_com1 <- tsne::tsne(X = vegdist(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct",pathogens_col_names], method="bray"),k = 4)
tsne_com1_df <- as.data.frame(tsne_com2)

tsne_com1_df$genotype <- full_dataset_local_no_out$genotype[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Oct"]
ggplot(data = tsne_com2_df, aes(x = V1, y=V2, color=genotype))+
  geom_point(size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))







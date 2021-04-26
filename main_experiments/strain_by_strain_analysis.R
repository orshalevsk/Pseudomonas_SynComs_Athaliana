library(ggplot2)
library(vegan)
library(lme4)
library(rstanarm)
library(pairwiseAdonis)
library(Hmisc)
library(dplyr)
library(RColorBrewer)
library(loo)

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


#merging Oct and Aug exp
full_dataset <- rbind(dataset_Oct[,colnames(dataset_Oct)], dataset_Aug[,colnames_Aug])


#subsetting to local only
full_dataset_local <- full_dataset[full_dataset$genotype%in%c("Lu","Eyach","HE","Kus","Tu-Wal","Schl") & full_dataset$treatment%in%c("Control","#1","#2","#3"),] 
full_dataset_local$genotype <- droplevels.factor(full_dataset_local$genotype) # dropping unwanted levels (Col0, etc...)
levels(full_dataset_local$genotype)<- c("Ey15-2","HE-1","Kus3-1","Lu3-30","Schl-7", "Tu-Wal-2") #renaming host genotypes to original 1001 names


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


#### adding patho and com under the mixed com
full_dataset_local_no_out$com3_patho[full_dataset_local_no_out$treatment=="#3"] <- rowSums(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",10:16])
full_dataset_local_no_out$com3_commensal[full_dataset_local_no_out$treatment=="#3"] <- rowSums(full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",17:23])

full_dataset_local_no_out$com3_patho[full_dataset_local_no_out$treatment%in%c("#1","#2")] <- full_dataset_local_no_out$bacterial_load[full_dataset_local_no_out$treatment%in%c("#1","#2")]
full_dataset_local_no_out$com3_commensal[full_dataset_local_no_out$treatment%in%c("#1","#2")] <- full_dataset_local_no_out$bacterial_load[full_dataset_local_no_out$treatment%in%c("#1","#2")]


#reordering the genotype levels - everything will be compared to Kus3-1 (random decision)
full_dataset_local_no_out$genotype <- factor(full_dataset_local_no_out$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))


## now building strain-by-strain dataset

source("~/ownCloud/My papers/Syncoms_paper/scripts/dataframe_to_strains_df_convertor.R")
dataset_by_strains <- dataframe_by_strain(dataframe = full_dataset_local_no_out)
dataset_by_strains <- dataset_by_strains[dataset_by_strains$strain_num%in%c(1:7) & dataset_by_strains$treatment%in%c("#1","#3") | dataset_by_strains$strain_num%in%c(11:17) & dataset_by_strains$treatment%in%c("#2","#3"),]
# now removing zeros. These cannot be log() transformed
for (i in colnames(dataset_by_strains)){
  zero_row <- which(dataset_by_strains[,i]==0)
  if (length(zero_row)>0){
    dataset_by_strains <- dataset_by_strains[-zero_row,]
  }
}

#reordering the genotype levels - everything will be compared to Kus3-1 (random decision)
dataset_by_strains$genotype <- factor(dataset_by_strains$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))


## Since both commensals and pathogens in community 3 behaved the same (for weight prediction) - plotting the correlation between commensals and pathogens in com#3.
## Result - perfectly correlated.

#Regardless of genotypes.
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/MixedCom_OTU5_non-OTU5_total_load_corr.pdf")
ggplot(aes(x=log(com3_commensal),y=log(com3_patho)), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",])+
  geom_point()+
  geom_smooth(method='lm') +
  theme_bw()
#! dev.off()

#simple peason corr test
cor.test(log(full_dataset_local_no_out$com3_commensal),y=log(full_dataset_local_no_out$com3_patho)) # R^2 = 0.912 , pval - super significance

#Genotype-wise.
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/MixedCom_OTU5_non-OTU5_total_load_corr_by_gen.pdf")
ggplot(aes(x=log(com3_commensal),y=log(com3_patho), color = genotype, group = genotype), 
       data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",])+
  geom_point()+
  geom_smooth(method='lm') +
  theme_bw() 
#! dev.off()


strain_type_by_genotype_stan_mod <- stan_glm(com3_patho ~ genotype*com3_commensal + exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",])
strain_type_by_genotype_stan_mod$stan_summary # doesn't seem that Ey15-2 has a different slope than the others by the model...

## Now let's examine correlation of all vs all strains, treatment-wise.
full_dataset_local_no_out_strains_only <- full_dataset_local_no_out[,c("exp","treatment","genotype",paste("X",c(1:7,11:17),sep = ""))]

# removing zeros. These cannot be log() transformed
for (i in colnames(full_dataset_local_no_out_strains_only)){
  zero_row <- which(full_dataset_local_no_out_strains_only[,i]==0)
  if (length(zero_row)>0){
    full_dataset_local_no_out_strains_only <- full_dataset_local_no_out_strains_only[-zero_row,]
  }
}

full_dataset_local_no_out_strains_com1 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#1",paste("X",c(1:7),sep = "")])
full_dataset_local_no_out_strains_com2 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#2",paste("X",c(11:17),sep = "")])
full_dataset_local_no_out_strains_com3 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#3",paste("X",c(1:7,11:17),sep = "")])



# calculate all vs. all correlations (with pval). Seems like a few strains have higher correlation than others, e.g. 14 and 15 (in both treatment 2 and 3)
library(reshape2)
corr_matriix_com1 <- melt(rcorr(as.matrix(full_dataset_local_no_out_strains_com1))[[1]])
corr_matriix_com2 <- melt(rcorr(as.matrix(full_dataset_local_no_out_strains_com2))[[1]])
corr_matriix_com3 <- melt(rcorr(as.matrix(full_dataset_local_no_out_strains_com3))[[1]])

#write.csv2(x = corr_matriix_com1,file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_correlations_table_com1_cytoscope.csv")
#write.csv2(x = corr_matriix_com2,file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_correlations_table_com2_cytoscope.csv")
#write.csv2(x = corr_matriix_com3,file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_correlations_table_com3_cytoscope.csv")

#now plot all, and color by experiment. CAN ALSO TRY BY GENOTYPE
cols_1 <- character(nrow(full_dataset_local_no_out_strains_com1))
cols_1[] <- "black"
cols_1[full_dataset_local_no_out_strains_only$exp[full_dataset_local_no_out_strains_only$treatment=="#1"]== "Aug"] <- "blue"
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_PathoCom.pdf")
pairs(full_dataset_local_no_out_strains_com1)
#! dev.off()

cols_2 <- character(nrow(full_dataset_local_no_out_strains_com2))
cols_2[] <- "black"
cols_2[full_dataset_local_no_out_strains_only$exp[full_dataset_local_no_out_strains_only$treatment=="#2"]== "Aug"] <- "blue"
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_CommensalCom.pdf")
pairs(full_dataset_local_no_out_strains_com2)
#! dev.off()

full_dataset_local_no_out_strains_com3 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#3" & full_dataset_local_no_out_strains_only$genotype=="Kus3-1",paste("X",c(1:7,11:17),sep = "")])
cols_3 <- character(nrow(full_dataset_local_no_out_strains_com3))
cols_3[] <- "black"
cols_3[full_dataset_local_no_out_strains_only$exp[full_dataset_local_no_out_strains_only$treatment=="#3" & full_dataset_local_no_out_strains_only$genotype=="Kus3-1"]== "Aug"] <- "blue"
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_MixedCom.pdf",width = 14,height = 14)
pairs(full_dataset_local_no_out_strains_com3, col = cols_3)
#! dev.off()


## now correlate abundance of each strain with total abudnance
library(tidyr)
library(ggpubr)
dataset_by_strains_and_abs_abundance <- full_dataset_local_no_out %>% gather(strain_num, strain_abundance, X1:X17)


dataset_by_strains_and_abs_abundance$strain_num <- factor(dataset_by_strains_and_abs_abundance$strain_num, levels = paste("X",c(1:7,11:17),sep = ""))
levels(dataset_by_strains_and_abs_abundance$strain_num) <- c("P1","P2","P3","P4","P5","P6","P7",
                                                             "C1","C2","C3","C4","C5","C6","C7")

#! pdf("~/ownCloud/My papers/Syncoms_paper/Final_edited_figures_first_review/Figures_sent_20200723/Figures/Supplementary_figures/FigureX_unedited.pdf")
ggplot(data = dataset_by_strains_and_abs_abundance[dataset_by_strains_and_abs_abundance$treatment=="#3",], aes(x = log(strain_abundance), y = log(bacterial_load)))+
  geom_point(alpha=0.2, stroke = 0.2) +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = -7) +
  stat_cor(p.accuracy = 0.0001, label.y = 4) +
  facet_wrap(. ~ strain_num) +
  ylab("log10 (cumulative isolate load)") +
  xlab("log10 (isolate abundance)") +
  theme_bw()
#! ßdev.off()

library(gplots)
colors = c(seq(-1,0,length=100),seq(0,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow"))(n = 199)

# distance & hierarchical clustering
distance= dist(as.matrix(cor(full_dataset_local_no_out_strains_com3), method ="euclidean"))    
hcluster = hclust(distance, method ="ward.D")


heatmap.2(as.matrix(cor(full_dataset_local_no_out_strains_com3)), col=my_palette, 
           density.info="none", trace="none",dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="none")


##### now the same strain-strain pairwise correlations, only with relative abundanc

full_dataset_local_no_out_strains_com1_relative <- full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#1",paste("X",c(1:7),sep = "")]/rowSums(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#1",paste("X",c(1:7),sep = "")])
full_dataset_local_no_out_strains_com2_relative <- full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#2",paste("X",c(11:17),sep = "")]/rowSums(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#2",paste("X",c(11:17),sep = "")])
full_dataset_local_no_out_strains_com3_relative <- full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#3",paste("X",c(1:7,11:17),sep = "")]/rowSums(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#3",paste("X",c(1:7,11:17),sep = "")])



# calculate all vs. all correlations (with pval). Seems like a few strains have higher correlation than others, e.g. 14 and 15 (in both treatment 2 and 3)
library(reshape2)
corr_matriix_com1_relative  <- melt(rcorr(as.matrix(full_dataset_local_no_out_strains_com1_relative ))[[1]])
corr_matriix_com2_relative  <- melt(rcorr(as.matrix(full_dataset_local_no_out_strains_com2_relative ))[[1]])
corr_matriix_com3_relative  <- melt(rcorr(as.matrix(full_dataset_local_no_out_strains_com3_relative ))[[1]])

#! write.csv2(x = corr_matriix_com1_relative ,file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_correlations_table_com1_relative_cytoscope.csv")
#! write.csv2(x = corr_matriix_com2_relative ,file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_correlations_table_com2_relative_cytoscope.csv")
#! write.csv2(x = corr_matriix_com3_relative ,file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_correlations_table_com3_relative_cytoscope.csv")



#now plot all, and color by experiment. CAN ALSO TRY BY GENOTYPE
cols_1 <- character(nrow(full_dataset_local_no_out_strains_com1_relative))
cols_1[] <- "black"
#cols_1[full_dataset_local_no_out_strains_only$exp[full_dataset_local_no_out_strains_only$treatment=="#1"]== "Aug"] <- "blue"
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_PathoCom_relative_adun.pdf")
pairs(full_dataset_local_no_out_strains_com1_relative)
#! dev.off()

cols_2 <- character(nrow(full_dataset_local_no_out_strains_com2_relative))
cols_2[] <- "black"
#cols_2[full_dataset_local_no_out_strains_only$exp[full_dataset_local_no_out_strains_only$treatment=="#2"]== "Aug"] <- "blue"
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_CommensalCom_relative_adun.pdf")
pairs(full_dataset_local_no_out_strains_com2_relative)
#! dev.off()

cols_3 <- character(nrow(full_dataset_local_no_out_strains_com3_relative))
cols_3 <- character(nrow(full_dataset_local_no_out_strains_com3_relative))
cols_3[] <- "black"
#cols_3[full_dataset_local_no_out_strains_only$exp[full_dataset_local_no_out_strains_only$treatment=="#3" & full_dataset_local_no_out_strains_only$genotype=="Kus3-1"]== "Aug"] <- "blue"
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_MixedCom_relative_adun.pdf",width = 14,height = 14)
pairs(full_dataset_local_no_out_strains_com3_relative, col = cols_3)
#! dev.off()

# now exploring the genotype*treatment effect for abundance change in MixedCom in comparison to exclusive communities.
fit_commensals_without_gen <- stan_glm(log(strain_abundance) ~ treatment+genotype + exp , data = dataset_by_strains[dataset_by_strains$treatment%in% c("#2","#3") & dataset_by_strains$strain_num %in% c(11:17),])
fit_commensals_with_gen <- stan_glm(log(strain_abundance) ~ treatment*genotype + exp , data = dataset_by_strains[dataset_by_strains$treatment%in% c("#2","#3") & dataset_by_strains$strain_num %in% c(11:17),])

fit_pathogens_without_gen <- stan_glm(log(strain_abundance) ~ treatment+genotype + exp , data = dataset_by_strains[dataset_by_strains$treatment%in% c("#1","#3") & dataset_by_strains$strain_num %in% c(1:7),])
fit_pathogens_with_gen <- stan_glm(log(strain_abundance) ~ treatment*genotype + exp , data = dataset_by_strains[dataset_by_strains$treatment%in% c("#1","#3") & dataset_by_strains$strain_num %in% c(1:7),])


loo_results_summary <- data.frame("strain_num"=NA,"model"=NA,"elpd_diff"=NA, "se_diff"=NA)
for (strain_num in c(1:7,11:17)){
  if (strain_num %in% c(1:7)){
    fit_without_gen <- stan_glm(log(strain_abundance) ~ treatment+genotype + exp, data = dataset_by_strains[dataset_by_strains$strain_num==strain_num & dataset_by_strains$treatment %in% c("#1","#3"),])
    fit_with_gen <- stan_glm(log(strain_abundance) ~ treatment*genotype + exp, data = dataset_by_strains[dataset_by_strains$strain_num==strain_num & dataset_by_strains$treatment %in% c("#1","#3"),])
  } else {
    fit_without_gen <- stan_glm(log(strain_abundance) ~ treatment+genotype + exp, data = dataset_by_strains[dataset_by_strains$strain_num==strain_num & dataset_by_strains$treatment %in% c("#2","#3"),])
    fit_with_gen <- stan_glm(log(strain_abundance) ~ treatment*genotype + exp, data = dataset_by_strains[dataset_by_strains$strain_num==strain_num & dataset_by_strains$treatment %in% c("#2","#3"),])
  }
  
  loo_output <- loo_compare(loo(fit_with_gen),loo(fit_without_gen))
  loo_df <- as.data.frame(loo_output)
  loo_df_temp <- data.frame("strain_num"=strain_num,"model"=rownames(loo_df),"elpd_diff"=loo_df$elpd_diff, "se_diff"=loo_df$se_diff)
  loo_results_summary <- rbind(loo_results_summary,loo_df_temp)
}
loo_results_summary <- loo_results_summary[!is.na(loo_results_summary),]
#! write.csv(x = loo_results_summary, file = "~/ownCloud/My papers/Syncoms_paper/Final_edited_figures/Supplementary/table_S2_model_per_strain_by_genotype.csv")

# now i will run Baesiyan model to understand the suppression of each strain (ignoring the genotype)
stan_summary_load_by_treat <- data.frame("strain"=NA, "mod_2.5" = NA,"mod_50" = NA, "mod_97.5" = NA)

chosen_dataset <- dataset_by_strains #choose data set to work with (subset to one exp?)
# first for pathogenic strains
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment*exp , data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#3") & chosen_dataset$strain_num %in% c(strain),])
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary["treatment#3",c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary)))
}
# now for commensal strains
for (strain in c(11:17)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment*exp , data = chosen_dataset[chosen_dataset$treatment%in% c("#2","#3") & chosen_dataset$strain_num %in% c(strain),])
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary["treatment#3",c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary)))
}

stan_summary_load_by_treat$strain <- factor(stan_summary_load_by_treat$strain)
stan_summary_load_by_treat <- stan_summary_load_by_treat[!is.na(stan_summary_load_by_treat$strain),]
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(1:7)] <- "OTU5"
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(11:17)] <- "non_OTU5"
stan_summary_load_by_treat$strain_type <- factor(stan_summary_load_by_treat$strain_type)

# now plotting the stan model results, taking into account both experiments
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/all_strains_change_com3.pdf", width = 7, height = 10, useDingbats = F)
ggplot(aes(x=strain,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) 
  #theme(axis.text.x = element_text(angle = 45)) 
#dev.off()


# now i will run Baesiyan model to understand the suppression of each strain, per genotype
stan_summary_load_by_treat <- data.frame("strain"=NULL, "mod_2.5" = NULL,"mod_50" = NULL, "mod_97.5" = NULL, "genotype"=NULL)

dataset_by_strains$genotype <- factor(dataset_by_strains$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))
chosen_dataset <- dataset_by_strains #choose data set to work with (subset to one exp?)

# first for pathogenic strains
for (gen in levels(chosen_dataset$genotype)){
  for (strain in c(1:7)){
    strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + exp , data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#3") & chosen_dataset$strain_num %in% c(strain) & chosen_dataset$genotype==gen ,])
    strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary["treatment#3",c("2.5%","50%","97.5%")]
    strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
    colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
    stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"genotype"=gen))
  }
}

# now for commensal strains
for (gen in levels(chosen_dataset$genotype)){
  for (strain in c(11:17)){
    strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + exp , data = chosen_dataset[chosen_dataset$treatment%in% c("#2","#3") & chosen_dataset$strain_num %in% c(strain) & chosen_dataset$genotype==gen,])
    strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary["treatment#3",c("2.5%","50%","97.5%")]
    strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
    colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
    stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"genotype"=gen))
  }
}



stan_summary_load_by_treat$strain <- factor(stan_summary_load_by_treat$strain)
stan_summary_load_by_treat <- stan_summary_load_by_treat[!is.na(stan_summary_load_by_treat$strain),]
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(1:7)] <- "OTU5"
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(11:17)] <- "non_OTU5"
stan_summary_load_by_treat$strain_type <- factor(stan_summary_load_by_treat$strain_type)

# now plotting the stan model results, taking into account both experiments
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/all_strains_change_com3_by_gen_dim.pdf", useDingbats = F)
ggplot(aes(x=strain,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1)+
  facet_grid(genotype ~ .)
#! dev.off()




# now i will run Baesiyan model to understand the suppression of each strain per genotype (interactions tern - genotype*treatment)
stan_summary_load_by_treat <- data.frame("strain"=NA, "mod_2.5" = NA,"mod_50" = NA, "mod_97.5" = NA, "treatment"= NA)

dataset_by_strains$genotype <- factor(dataset_by_strains$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))
chosen_dataset <- dataset_by_strains #choose data set to work with (subset to one exp?)
# first for pathogenic strains
coef_tested <- c("treatment#3:genotypeHE-1","treatment#3:genotypeEy15-2","treatment#3:genotypeLu3-30", "treatment#3:genotypeSchl-7", "treatment#3:genotypeTu-Wal-2")
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ exp*treatment*genotype, data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#3") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(strain_abundance_summary)
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=rownames(strain_abundance_summary)))
}
# now for commensal strains
for (strain in c(11:17)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ exp*treatment*genotype, data = chosen_dataset[chosen_dataset$treatment%in% c("#2","#3") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(strain_abundance_summary)
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=rownames(strain_abundance_summary)))
}

stan_summary_load_by_treat$strain <- factor(stan_summary_load_by_treat$strain)
stan_summary_load_by_treat <- stan_summary_load_by_treat[!is.na(stan_summary_load_by_treat$strain),]
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(1:7)] <- "OTU5"
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(11:17)] <- "non_OTU5"
stan_summary_load_by_treat$strain_type <- factor(stan_summary_load_by_treat$strain_type)



gen_split <- strsplit(stan_summary_load_by_treat$treatment, ":genotype")
even <- function(x) x%%2 == 0
stan_summary_load_by_treat$genotype <- unlist(gen_split)[even(1:length(gen_split))]


# now plotting the stan model results, taking into account both experiments and this time also the genotype interactions with MixedCom
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/all_strains_change_com3_by_genotype_updated_dim.pdf", useDingbats = F)
ggplot(aes(x=strain,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  scale_color_manual(values = c("#a50026","#4393c3")) +
  facet_grid(genotype ~ .) 
#! dev.off()

# now plotting the stan model results, taking into account both experiments and this time also the genotype - only pathogens
stan_summary_load_by_treat_patho <- stan_summary_load_by_treat[stan_summary_load_by_treat$strain%in%c(1:7),]
pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/patho_strains_change_com3_by_genotype.pdf",width = 7, height = 9)
ggplot(aes(x=treatment,y=stan_summary_load_by_treat_patho$mod_50, ymin = stan_summary_load_by_treat_patho$mod_2.5, ymax = stan_summary_load_by_treat_patho$mod_97.5),
       data = stan_summary_load_by_treat_patho)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(strain ~ .)
dev.off()

# now plotting the stan model results, taking into account both experiments and this time also the genotype - only commensals
stan_summary_load_by_treat_comm <- stan_summary_load_by_treat[stan_summary_load_by_treat$strain%in%c(11:17),]
pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/comm_strains_change_com3_by_genotype.pdf",width = 7, height = 9)
ggplot(aes(x=treatment,y=stan_summary_load_by_treat_comm$mod_50, ymin = stan_summary_load_by_treat_comm$mod_2.5, ymax = stan_summary_load_by_treat_comm$mod_97.5),
       data = stan_summary_load_by_treat_comm)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(strain ~ .)
dev.off()


#lastly, i will validate microbe-microbe correlation (all pair-wise) and will account for the genotype, WITHIN commuities
strains_col_names <- c("X1","X2","X3","X4","X5","X6","X7","X11","X12","X13","X14","X15","X16","X17")
pathogens_col_names <- c("X1","X2","X3","X4","X5","X6","X7")
commensals_col_names <- c("X11","X12","X13","X14","X15","X16","X17")


chosen_dataset <- full_dataset_local_no_out_strains_only #choose data set to work with (subset to one exp?)

# first for PathoCom
stan_summary_pariwise_strains_by_gen_PathoCom <- data.frame(matrix(NA, nrow = 7, ncol = 7))
colnames(stan_summary_pariwise_strains_by_gen_PathoCom) <- pathogens_col_names

for (strain_focal in pathogens_col_names){
  stan_summary_pariwise_strains_by_gen <- data.frame("strain"=NA, "mod_50" = NA)
  current_patho_subset <- pathogens_col_names[!pathogens_col_names %in% strain_focal]
  for (strain_formula in current_patho_subset){
    current_formula <- paste(strain_focal,paste(strain_formula, collapse = "+"),sep = "~", collapse = "")
    current_formula <- as.formula(paste(current_formula,"+ genotype*exp"))
    strain_abundance_by_treat_mod <- stan_glm(current_formula , data = chosen_dataset[chosen_dataset$treatment%in% c("#1"),])
    temp_df <- as.data.frame(t(c(strain_formula,strain_abundance_by_treat_mod$stan_summary[strain_formula,c("2.5%")])))
    colnames(temp_df) <- c("strain", "mod_50")
    stan_summary_pariwise_strains_by_gen <- rbind(stan_summary_pariwise_strains_by_gen,temp_df)
  }
  stan_summary_pariwise_strains_by_gen <- stan_summary_pariwise_strains_by_gen[!is.na(stan_summary_pariwise_strains_by_gen$strain),c("strain","mod_50")]
  own_strain <- data.frame("strain"=strain_focal, "mod_50"=1)
  stan_summary_pariwise_strains_by_gen <- rbind(stan_summary_pariwise_strains_by_gen,own_strain)
  stan_summary_pariwise_strains_by_gen_PathoCom[,strain_focal] <- stan_summary_pariwise_strains_by_gen$mod_50[match(pathogens_col_names, stan_summary_pariwise_strains_by_gen$strain)]
}

# now for CommenCom
stan_summary_pariwise_strains_by_gen_CommenCom <- data.frame(matrix(NA, nrow = 7, ncol = 7))
colnames(stan_summary_pariwise_strains_by_gen_CommenCom) <- commensals_col_names

for (strain_focal in commensals_col_names){
  stan_summary_pariwise_strains_by_gen <- data.frame("strain"=NA, "mod_50" = NA)
  current_Commen_subset <- commensals_col_names[!commensals_col_names %in% strain_focal]
  for (strain_formula in current_Commen_subset){
    current_formula <- paste(strain_focal,paste(strain_formula, collapse = "+"),sep = "~", collapse = "")
    current_formula <- as.formula(paste(current_formula,"+ genotype*exp"))
    strain_abundance_by_treat_mod <- stan_glm(current_formula , data = chosen_dataset[chosen_dataset$treatment%in% c("#2"),])
    temp_df <- as.data.frame(t(c(strain_formula,strain_abundance_by_treat_mod$stan_summary[strain_formula,c("2.5%")])))
    colnames(temp_df) <- c("strain", "mod_50")
    stan_summary_pariwise_strains_by_gen <- rbind(stan_summary_pariwise_strains_by_gen,temp_df)
  }
  stan_summary_pariwise_strains_by_gen <- stan_summary_pariwise_strains_by_gen[!is.na(stan_summary_pariwise_strains_by_gen$strain),c("strain","mod_50")]
  own_strain <- data.frame("strain"=strain_focal, "mod_50"=1)
  stan_summary_pariwise_strains_by_gen <- rbind(stan_summary_pariwise_strains_by_gen,own_strain)
  stan_summary_pariwise_strains_by_gen_CommenCom[,strain_focal] <- stan_summary_pariwise_strains_by_gen$mod_50[match(commensals_col_names, stan_summary_pariwise_strains_by_gen$strain)]
}



# Lastly for MixedCom
stan_summary_pariwise_strains_by_gen_MixedCom <- data.frame(matrix(NA, nrow = 14, ncol = 14))
colnames(stan_summary_pariwise_strains_by_gen_MixedCom) <- strains_col_names

for (strain_focal in strains_col_names){
  stan_summary_pariwise_strains_by_gen <- data.frame("strain"=NA, "mod_50" = NA)
  current_patho_subset <- strains_col_names[!strains_col_names %in% strain_focal]
  for (strain_formula in current_patho_subset){
    current_formula <- paste(strain_focal,paste(strain_formula, collapse = "+"),sep = "~", collapse = "")
    current_formula <- as.formula(paste(current_formula,"+ genotype*exp"))
    strain_abundance_by_treat_mod <- stan_glm(current_formula , data = chosen_dataset[chosen_dataset$treatment%in% c("#3"),])
    temp_df <- as.data.frame(t(c(strain_formula,strain_abundance_by_treat_mod$stan_summary[strain_formula,c("2.5%")])))
    colnames(temp_df) <- c("strain", "mod_50")
    stan_summary_pariwise_strains_by_gen <- rbind(stan_summary_pariwise_strains_by_gen,temp_df)
  }
  stan_summary_pariwise_strains_by_gen <- stan_summary_pariwise_strains_by_gen[!is.na(stan_summary_pariwise_strains_by_gen$strain),c("strain","mod_50")]
  own_strain <- data.frame("strain"=strain_focal, "mod_50"=1)
  stan_summary_pariwise_strains_by_gen <- rbind(stan_summary_pariwise_strains_by_gen,own_strain)
  stan_summary_pariwise_strains_by_gen_MixedCom[,strain_focal] <- stan_summary_pariwise_strains_by_gen$mod_50[match(strains_col_names, stan_summary_pariwise_strains_by_gen$strain)]
}


# now plotting the stan model results, taking into account both experiments
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/all_strains_change_com3.pdf", width = 7, height = 10)
ggplot(aes(x=strain,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  theme(axis.text.x = element_text(angle = 45)) 
#! dev.off()

#! write.csv2(x = stan_summary_pariwise_strains_by_gen_PathoCom, file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_PathoCom_STAN_genotype*exp_cdl50%.csv")
#! write.csv2(x = stan_summary_pariwise_strains_by_gen_CommenCom, file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_CommenCom_STAN_genotype*exp_cdl50%.csv")
#! write.csv2(x = stan_summary_pariwise_strains_by_gen_MixedCom, file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_MixedCom_STAN_genotype*exp_cdl50%.csv")
stan_summary_pariwise_strains_by_gen_PathoCom <- sapply(stan_summary_pariwise_strains_by_gen_PathoCom, as.numeric)
heatmap.2(as.matrix(stan_summary_pariwise_strains_by_gen_PathoCom), trace='none')




library(reshape2)
corr_matriix_com1_stan_mod  <- melt(as.matrix(stan_summary_pariwise_strains_by_gen_PathoCom))
corr_matriix_com2_stan_mod  <- melt(as.matrix(stan_summary_pariwise_strains_by_gen_CommenCom))
corr_matriix_com3_stan_mod  <- melt(as.matrix(stan_summary_pariwise_strains_by_gen_MixedCom))

#! write.csv2(x = corr_matriix_com1_stan_mod, file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_PathoCom_STAN_genotype*exp_cdl50%_cytoscope.csv")
#! write.csv2(x = corr_matriix_com2_stan_mod, file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_CommenCom_STAN_genotype*exp_cdl50%_cytoscope.csv")
#! write.csv2(x = corr_matriix_com3_stan_mod, file = "~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/all_pairwise_MixedCom_STAN_genotype*exp_cdl50%_cytoscope.csv")





full_dataset_local_no_out_strains_com1 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#1" & full_dataset_local_no_out_strains_only$genotype=="Ey15-2",paste("X",c(1:7),sep = "")])
full_dataset_local_no_out_strains_com2 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#2",paste("X",c(11:17),sep = "")])
full_dataset_local_no_out_strains_com3 <- log(full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#3" & full_dataset_local_no_out_strains_only$genotype=="Ey15-2" & full_dataset_local_no_out_strains_only$exp=="Aug",paste("X",c(1:7,11:17),sep = "")])



# calculate all vs. all correlations (with pval). Seems like a few strains have higher correlation than others, e.g. 14 and 15 (in both treatment 2 and 3)
rcorr(as.matrix(full_dataset_local_no_out_strains_com1))
rcorr(as.matrix(full_dataset_local_no_out_strains_com2))
rcorr(as.matrix(full_dataset_local_no_out_strains_com3))



################ now i will analyze strain-strain interactions in HE, individual commensals + PathoCom experiment #######
individual_commen_pathogens_dataset <- full_dataset[full_dataset$treatment %in% c("#1","#3","#1+s13", "#1+s14", "#1+s15", "#1+s17") & full_dataset$exp=="Oct",]
individual_commen_pathogens_dataset <- individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$genotype%in%c("HE"),] 
individual_commen_pathogens_dataset$genotype <- droplevels.factor(individual_commen_pathogens_dataset$genotype) # dropping unwanted levels (Col0, etc...)
levels(individual_commen_pathogens_dataset$genotype)<- c("HE-1") #renaming host genotypes to original 1001 names

# now removing zeros. These cannot be log() transformed
#for (i in colnames(individual_commen_pathogens_dataset)){
#  zero_row <- which(individual_commen_pathogens_dataset[,i]==0)
#  if (length(zero_row)>0){
#    individual_commen_pathogens_dataset <- individual_commen_pathogens_dataset[-zero_row,]
#  }
#}

for (i in (paste("X",c(1:7,11:17),sep = ""))){
  zero_row <- individual_commen_pathogens_dataset[,i]==0
  zero_row <- zero_row[!is.na(zero_row)]
  if (length(zero_row)>0){
    individual_commen_pathogens_dataset[zero_row,i] <- 0.0000000000001
  }
}



individual_commen_pathogens_dataset_com1 <- log(individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$treatment=="#1" ,paste("X",c(1:7),sep = "")])
individual_commen_pathogens_dataset_com3 <- log(individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$treatment=="#3",paste("X",c(1:7,13,14,15,17),sep = "")])
individual_commen_pathogens_dataset_com1_s13 <- log(individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$treatment=="#1+s13",paste("X",c(1:7,13),sep = "")])
individual_commen_pathogens_dataset_com1_s14 <- log(individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$treatment=="#1+s14",paste("X",c(1:7,14),sep = "")])
individual_commen_pathogens_dataset_com1_s15 <- log(individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$treatment=="#1+s15",paste("X",c(1:7,15),sep = "")])
individual_commen_pathogens_dataset_com1_s17 <- log(individual_commen_pathogens_dataset[individual_commen_pathogens_dataset$treatment=="#1+s17",paste("X",c(1:7,17),sep = "")])

colnames(individual_commen_pathogens_dataset_com1) <- c(paste("P",c(1:7),sep = ""))
colnames(individual_commen_pathogens_dataset_com3) <- c(paste("P",c(1:7),sep = ""),paste("C",c(3,4,5,7),sep = ""))
colnames(individual_commen_pathogens_dataset_com1_s13) <- c(paste("P",c(1:7),sep = ""),paste("C",c(3),sep = ""))
colnames(individual_commen_pathogens_dataset_com1_s14) <- c(paste("P",c(1:7),sep = ""),paste("C",c(4),sep = ""))
colnames(individual_commen_pathogens_dataset_com1_s15) <- c(paste("P",c(1:7),sep = ""),paste("C",c(5),sep = ""))
colnames(individual_commen_pathogens_dataset_com1_s17) <- c(paste("P",c(1:7),sep = ""),paste("C",c(7),sep = ""))


# calculate all vs. all correlations (with pval). 
r <- rcorr(as.matrix(individual_commen_pathogens_dataset_com1))

write.csv2(x = melt(rcorr(as.matrix(individual_commen_pathogens_dataset_com1))[[1]]),"~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/individual_commensals_exp_com1_cytoscope.csv")
write.csv2(x = melt(rcorr(as.matrix(individual_commen_pathogens_dataset_com3))[[1]]),"~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/individual_commensals_exp_com3_cytoscope.csv")
write.csv2(x = melt(rcorr(as.matrix(individual_commen_pathogens_dataset_com1_s13))[[1]]),"~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/individual_commensals_exp_com1-C3_cytoscope.csv")
write.csv2(x = melt(rcorr(as.matrix(individual_commen_pathogens_dataset_com1_s14))[[1]]),"~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/individual_commensals_exp_com1-C4_cytoscope.csv")
write.csv2(x = melt(rcorr(as.matrix(individual_commen_pathogens_dataset_com1_s15))[[1]]),"~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/individual_commensals_exp_com1-C5_cytoscope.csv")
write.csv2(x = melt(rcorr(as.matrix(individual_commen_pathogens_dataset_com1_s17))[[1]]),"~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Within_treatments/individual_commensals_exp_com1-C7_cytoscope.csv")


pairs(as.matrix(individual_commen_pathogens_dataset_com1))
pairs(as.matrix(individual_commen_pathogens_dataset_com3))
pairs(as.matrix(individual_commen_pathogens_dataset_com1_s13))
pairs(as.matrix(individual_commen_pathogens_dataset_com1_s14))
pairs(as.matrix(individual_commen_pathogens_dataset_com1_s15))
pairs(as.matrix(individual_commen_pathogens_dataset_com1_s17))


# now building strain-by-strain dataset
source("~/ownCloud/My papers/Syncoms_paper/scripts/dataframe_to_strains_df_convertor_individual_commensals.R")
dataset_by_strains_individual_commen <- dataframe_by_strain(dataframe = individual_commen_pathogens_dataset)
dataset_by_strains_individual_commen <- dataset_by_strains_individual_commen[dataset_by_strains_individual_commen$strain_num%in%c(1:7) & dataset_by_strains_individual_commen$treatment%in%c("#1","#3","#1+s13", "#1+s14", "#1+s15", "#1+s17") | 
                                                                               dataset_by_strains_individual_commen$strain_num%in%c(11:17) & dataset_by_strains_individual_commen$treatment%in%c("#3") | 
                                                                               dataset_by_strains_individual_commen$strain_num==13 & dataset_by_strains_individual_commen$treatment=="#1+s13" |
                                                                               dataset_by_strains_individual_commen$strain_num==14 & dataset_by_strains_individual_commen$treatment=="#1+s14" |
                                                                               dataset_by_strains_individual_commen$strain_num==15 & dataset_by_strains_individual_commen$treatment=="#1+s15" |
                                                                               dataset_by_strains_individual_commen$strain_num==17 & dataset_by_strains_individual_commen$treatment=="#1+s17",]

# now i will run Baesiyan model to understand the suppression of each strain per genotype
stan_summary_load_by_treat <- data.frame("strain"=NA, "mod_2.5" = NA,"mod_50" = NA, "mod_97.5" = NA, "treatment"= NA)

chosen_dataset <- dataset_by_strains_individual_commen #choose data set to work with (subset to one exp?)
# first for pathogenic strains
coef_tested <- "treatment#3"
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + tray , data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#3") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=coef_tested))
}
# now for pathogens +s13 
coef_tested <- "treatment#1+s13"
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + tray, data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#1+s13") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=coef_tested))
}

# now for pathogens +s14 
coef_tested <- "treatment#1+s14"
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + tray, data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#1+s14") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=coef_tested))
}


# now for pathogens +s15 
coef_tested <- "treatment#1+s15"
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + tray , data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#1+s15") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=coef_tested))
}


# now for pathogens +s17 
coef_tested <- "treatment#1+s17"
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(log(strain_abundance) ~ treatment + tray , data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#1+s17") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary),"treatment"=coef_tested))
}


stan_summary_load_by_treat$strain <- factor(stan_summary_load_by_treat$strain)
stan_summary_load_by_treat <- stan_summary_load_by_treat[!is.na(stan_summary_load_by_treat$strain),]
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(1:7)] <- "OTU5"
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(11:17)] <- "non_OTU5"
stan_summary_load_by_treat$strain_type <- factor(stan_summary_load_by_treat$strain_type)


# now plotting the stan model results, taking into account both experiments and this time also the genotype
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/pathogens_change_individual_commensals.pdf",width = 7, height = 14)
ggplot(aes(x=treatment,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(strain ~ .)
#! dev.off()

ggplot(aes(x = strain_num, y= log(strain_abundance), fill=treatment),data = dataset_by_strains_individual_commen[dataset_by_strains_individual_commen$strain_num %in% c(1:7),])+
  geom_boxplot()
  
# now plotting the stan model results, per strain
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/pathogens_change_individual_commensals.pdf",width = 7, height = 14)
ggplot(aes(x=treatment,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(strain ~ .)
#! dev.off()

# now the same, but not strain by strain, but OTU5 in general
individual_commen_pathogens_dataset$patho_load <- rowSums(individual_commen_pathogens_dataset[,paste("X",c(1:7),sep = "")])

strain_abundance_by_treat_mod <- stan_glm(log(patho_load) ~ treatment + tray , data = individual_commen_pathogens_dataset, seed = 12345)
strain_abundance_summary <- as.data.frame(strain_abundance_summary)
strain_abundance_summary$treatment <- rownames(strain_abundance_summary)

# now plotting the stan model results, all OTU5
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/pathogens_change_individual_commensals.pdf",width = 7, height = 14)
ggplot(aes(x=treatment,y=strain_abundance_summary$`50%`, ymin = strain_abundance_summary$`2.5%`, ymax = strain_abundance_summary$`97.5%`),
       data = strain_abundance_summary)+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()


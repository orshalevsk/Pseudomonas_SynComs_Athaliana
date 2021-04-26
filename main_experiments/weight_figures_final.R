'''
This script produces all the original figures for Shalev et al. 2020 paper (Syncoms and host genotypes research)
'''

#loading libraries
library(ggplot2)
library(dabestr)
library(lme4)
library(rstan)
library(rstanarm)
library(shinystan)
library(bayestestR)

#loading and configuring Aug exp
setwd("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/weight+PAMP-triggered_August_2018/results/")
dataset_Aug <- read.csv2("weight_full_dataset.csv")
dataset_Aug$weight <- as.numeric(as.character(dataset_Aug$weight))
dataset_Aug <- dataset_Aug[!is.na(dataset_Aug$weight),]
dataset_Aug$treatment <- factor(dataset_Aug$treatment, levels = c("Control","#1","#1_boil","#2","#3"))
dataset_Aug$exp <- "Aug"
dataset_Aug$exp <- as.factor(dataset_Aug$exp)
dataset_Aug$weight <- as.numeric(as.character(dataset_Aug$weight))
dataset_Aug$genotype <- factor(dataset_Aug$genotype, levels = c("Eyach","Kus","Schl","Lu","HE","Tu-Wal","Tu-Wal_ics1","Tu-Wal_acd6","Col0","Col0_ics1","Col0_acd6"))
levels(dataset_Aug$genotype)[levels(dataset_Aug$genotype)=="Col0_ics1"] <- "sid2-2"
dataset_Aug$tray <- as.factor(dataset_Aug$tray)
dataset_Aug$pot <- as.factor(dataset_Aug$pot)

#loading and configuring Oct exp
default_dir <- "~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/weight+HE_October_2018/"
setwd(default_dir)
dataset_oct <- read.csv2("results/dataset_before_weight.csv")
dataset_oct$exp <- "Oct"
weight_oct <- read.csv2("results/weight_20181128.csv")
dataset_oct <- merge(dataset_oct, weight_oct, by = "position")
dataset_oct$weight <- as.numeric(as.character(dataset_oct$weight))
dataset_oct$exp <- as.factor(dataset_oct$exp)
dataset_oct$genotype <- factor(dataset_oct$genotype, levels = c("Eyach","Kus","Schl","Lu","HE","Tu-Wal","Tu-Wal_ics1","Tu-Wal_acd6","Col0","sid2-2","Col0_acd6"))
dataset_oct$treatment <- factor(dataset_oct$treatment, levels = c("Control", "#1", "#2","#3","#1+s13","#1+s14","#1+s15","#1+s17"))
dataset_oct$position <- paste("O",dataset_oct$position,sep = "")
colnames(dataset_oct)[colnames(dataset_oct)=="weight.only"] <- "weight.only."
dataset_oct$tray <- as.factor(paste(dataset_oct$tray.x,"oct"))
dataset_oct$pot <- as.factor(dataset_oct$pot.x)

#merging Oct and Aug exp
full_dataset <- rbind(dataset_oct[,c("treatment","genotype","exp","weight","tray","pot")], dataset_Aug[,c("treatment","genotype","exp","weight","tray","pot")])

#subsetting to local only
full_dataset_local <- full_dataset[full_dataset$genotype%in%c("Lu","Eyach","HE","Kus","Tu-Wal","Schl") & full_dataset$treatment%in%c("Control","#1","#2","#3"),] 
full_dataset_local$genotype <- droplevels.factor(full_dataset_local$genotype) # dropping unwanted levels (Col0, etc...)
levels(full_dataset_local$genotype)<- c("Ey15-2","Kus3-1","Schl-7","Lu3-30","HE-1", "Tu-Wal-2") #renaming host genotypes to original 1001 names

#removing outliers function - more than 2.5 sds from mean weight of genotype*treatment*experiment
outliers_removal <- function(dataset){
  outliers_list <- vector()
  for (exp in unique(dataset$exp)){
    for (gen in unique(dataset$genotype)){
      for (treat in unique(dataset$treatment)){
        gen_by_treat_by_exp <- dataset$treatment==treat & dataset$genotype==gen & dataset$exp==exp
        mean_gen_by_treat_by_exp <- mean(dataset$weight[gen_by_treat_by_exp])
        sd_2.5 <-2.5*sd(dataset$weight[gen_by_treat_by_exp])
        outliers <- dataset$weight[gen_by_treat_by_exp] > (mean_gen_by_treat_by_exp+sd_2.5) | dataset$weight[gen_by_treat_by_exp] <(mean_gen_by_treat_by_exp-sd_2.5)
        if (sum(outliers)>0){
          outliers_list <- c(outliers_list,row.names(dataset[gen_by_treat_by_exp,][outliers,]))
        }
      }
    }
  }
  return(outliers_list)
}

outliers_list <- outliers_removal(dataset = full_dataset_local) # removing outliers
length(outliers_list)/length(full_dataset_local$weight) # Sanity check: the % of outliers removed. about 1%, as expected (2.5 std are 99% of the normal distrubiton)
full_dataset_local_no_out <- full_dataset_local[-which(row.names(full_dataset_local) %in% outliers_list),]
############################## Figure 1 + Supplementary 1 - Weight reduction by pathogens and protection by commensals are host genotype-dependent #############################################


#####Fig 1A - representative exp for weight of host*syncom (Aug results) 
full_dataset_local_rep_exp$gen_by_treat <- paste(full_dataset_local_rep_exp$genotype,full_dataset_local_rep_exp$treatment) # with outliers
full_dataset_local_rep_exp_no_out$gen_by_treat <- paste(full_dataset_local_rep_exp_no_out$genotype,full_dataset_local_rep_exp_no_out$treatment) #without outliers

#ploting using dabestr - Baesyian based approach

index_dabest <- list() # creating an index for the dabest function (what is control, and to seperate by genotype)
i = 1 # helping index for the loop that will come
sorted_genotypes <- c("Lu3-30","Kus3-1","Schl-7","HE-1", "Tu-Wal-2","Ey15-2") # sorting the genotypes for the plot. I used Tu-Wal before the end, and Eyach at the end to mark its difference 

for (gen in sorted_genotypes){
  current_genotype <- paste(gen, c("Control","#1","#2","#3")) # hard coding - if i change the treatments name - i must change it here also. Did it to sort in the right way
  index_dabest[[i]] <- current_genotype
  i=i+1
}

unpaired_mean_diff <- dabest(full_dataset_local_rep_exp_no_out, gen_by_treat, weight, paired = F, idx = index_dabest, seed = 12345, ci = 95)
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_genotype_by_treatments/figure_weight_all_treatments_by_genotypes.pdf", width = 14, height = 6)
plot(unpaired_mean_diff, color.column = treatment)
#! dev.off()
results_print <- as.data.frame(unpaired_mean_diff$result[,c("test_group","difference","control_size","bca_ci_low","bca_ci_high")])
write.csv2(x = results_print,file = "~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_genotype_by_treatments/dabest_output.csv")




temp_mod <- BayesFactor::lmBF(formula = weight ~ genotype, data = full_dataset_local_rep_exp[full_dataset_local_rep_exp_no_out$treatment=="#1",], posterior = T, iter=10000)
summary(temp_mod)[[2]]

temp_mod <- BayesFactor::lmBF(formula = weight ~ genotype + exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",], posterior = T, iter=10000)



#####Fig 1B - mixled linear model using all data (both exp), and accounting for tray, pot and exp as random effects to find genotype*treatment interactions
##i.e. proving that eyach is an outlier for 'Mixed' com, using all meatdata and after correcting for it

full_dataset_local_no_out$genotype <- factor(full_dataset_local_no_out$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))
# frequentist models first (before Baesiyan). This part is to find the best model without too much computation time (models comparison part)

weight_mod <- lmer(data = full_dataset_local_no_out, formula = weight ~ genotype + treatment + exp + (1|pot) + (1|tray),REML = F) 
weight_mod1 <- lmer(data = full_dataset_local_no_out, formula = weight ~ genotype*treatment + exp + (1|pot) + (1|tray),REML = F) # with g*t interactions
anova(weight_mod, weight_mod1) # comparing the models - proving that genotype*treatment interactions improves the model, hence neccesary

weight_mod2 <- lmer(data = full_dataset_local_no_out, formula = weight ~ genotype*treatment + exp + exp*genotype + (1|pot) + (1|tray),REML = F) # addind exp*genotype interactions
anova(weight_mod1, weight_mod2) # comparing the models - proving that exp*genotype interactions improves the model, hence neccesary

weight_mod3 <- lmer(data = full_dataset_local_no_out, formula = weight ~ genotype*treatment + exp + exp*genotype + exp*treatment + (1|pot) + (1|tray),REML = F) # addind exp*genotype interactions
anova(weight_mod2, weight_mod3) # comparing the models - proving that exp*treatment interactions improves the model, hence neccesary

bayesfactor_models(weight_mod1 , weight_mod2, weight_mod3 , denominator = weight_mod)
bayesfactor_models(weight_mod1 , weight_mod2, weight_mod3 , denominator = weight_mod)

# Decided to keep on with model #3. It corrects for the genotypes and treamtments differences due to the experiment.
# Not g*t due to exp, but the latter will be too convoluted to interpret.

#Baesiyan part

#This is quite a slow process (15-20 minutes) - Baesiyan lmer model with many covariables. Run manually only when needed
#! weight_mod_stan <- stan_lmer(data = full_dataset_local_no_out, formula = weight ~ genotype*treatment + exp + exp*genotype + exp*treatment + (1|pot) + (1|tray) ,seed=12345, control = list(adapt_delta = 0.9999),iter=5000)

#This model is faster but less accurate. Good only when trying out things, and not for the final figure! 
#! weight_mod_stan <- stan_lmer(data = full_dataset_local_no_out, formula = weight ~ genotype*treatment + exp + exp*genotype + exp*treatment + (1|pot) + (1|tray) ,seed=12345) 
#! launch_shinystan(weight_mod_stan)
weight_mod_stan$stan_summary # sanity check. Check for Rhat. should be less than 1.01 (which is the case :) )
prior_summary(weight_mod_stan)

coeffcients_interactions <- c("genotypeEy15-2:treatment#1",	"genotypeSchl-7:treatment#1",	"genotypeLu3-30:treatment#1",	"genotypeHE-1:treatment#1",	"genotypeTu-Wal-2:treatment#1",
  "genotypeEy15-2:treatment#2",	"genotypeSchl-7:treatment#2",	"genotypeLu3-30:treatment#2",	"genotypeHE-1:treatment#2",	"genotypeTu-Wal-2:treatment#2",	
  "genotypeEy15-2:treatment#3",	"genotypeSchl-7:treatment#3",	"genotypeLu3-30:treatment#3",	"genotypeHE-1:treatment#3",	"genotypeTu-Wal-2:treatment#3") # all names are here rownames(weight_mod_stan$stan_summary)

coeffcients_interactions_summary <- weight_mod_stan$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")]
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)  

gen_by_treat_name <- strsplit(coeffcients_interactions_summary$condition,":")
for (i in 1:length(gen_by_treat_name)){
  gen <- gen_by_treat_name[i][[1]][1]
  treat <- gen_by_treat_name[i][[1]][2]
  coeffcients_interactions_summary$genotype[i] <- gen
  coeffcients_interactions_summary$treatment[i] <- treat
}



# Plotting the effect of each genotype * treatment, all in comparison to genotype Kus. This is a bit tricky to interpret since Kus will not be there. 
# At the end i am producing 95% CI for all effects. Focusing on Eyach - #3. Again - this is accross both exp, and accounting for exp, pot and tray as random effects

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_genotype_by_treatments/supplementary_weight_all_treatments_by_genotypes_stan_model.pdf", width = 7,height = 10, useDingbats = F)
ggplot(data = coeffcients_interactions_summary, aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(aes(group=treatment, color=treatment),position=position_dodge(width=1)) +
  theme_bw() +
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_grid(treatment ~ .) 
#! dev.off()


weight_mod_stan <- stan_glm(data = full_dataset_local_no_out, formula = weight ~ treatment + exp ,seed=12345) 

coeffcients_treament <- c("treatment#1",	"treatment#2",	"treatment#3")
coeffcients_treament <- weight_mod_stan$stan_summary[coeffcients_treament,c("2.5%","50%","97.5%")]
coeffcients_treament_summary <- as.data.frame(coeffcients_treament)
coeffcients_treament_summary$condition <- rownames(coeffcients_treament_summary)  


# Plotting the effect of each treatment, all in comparison to Control.
# At the end i am producing 95% CI for all effects.

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_genotype_by_treatments/supplementary_weight_by_treatment_no_gen.pdf", width = 7,height = 10, useDingbats = F)
ggplot(data = coeffcients_treament_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(aes(group=condition, color=condition),position=position_dodge(width=1)) +
  theme_bw() +
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  geom_hline(yintercept=0, linetype="dashed")
#! dev.off()



treatment_by_gen_df <- data.frame()
for (treat in c("#1","#2","#3")){
  weight_mod_stan <- stan_glm(data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("Control",treat),], formula = weight ~ exp + treatment*genotype ,seed=12345) 
  gen_names <- paste(paste("treatment",treat,":", sep = ""),c("genotypeEy15-2","genotypeSchl-7","genotypeLu3-30","genotypeHE-1","genotypeTu-Wal-2"), sep = "")
  temp_df <- as.data.frame(weight_mod_stan$stan_summary[gen_names,c("2.5%","50%","97.5%")])
  temp_df$genotype <- rownames(temp_df)  
  temp_df$treatment <- treat
  treatment_by_gen_df <- rbind(treatment_by_gen_df, temp_df)
}

gen_by_treat_name <- strsplit(treatment_by_gen_df$genotype,"genotype")
for (i in 1:length(gen_by_treat_name)){
  gen <- gen_by_treat_name[i][[1]][2]
  treatment_by_gen_df$genotype_updated[i] <- gen
}


# Plotting the effect of each treatment, all in comparison to Control.
# At the end i am producing 95% CI for all effects.
ggplot(data = treatment_by_gen_df, aes(x=genotype_updated,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(aes(group=genotype_updated, color=genotype_updated),position=position_dodge(width=1)) +
  theme_bw() +
  facet_grid(treatment ~ .)






### another model to understand treatment*genotype, using bayes factor
weight_mod_stan <- BayesFactor::lmBF(data = full_dataset_local_no_out, formula = weight ~ genotype*treatment*exp + tray,posterior = T, iter=10000) 
weight_mod_stan_df <- as.data.frame(summary(weight_mod_stan)[[2]])
weight_mod_stan_df$condition <- rownames(weight_mod_stan_df)

gen_by_treat_name <- strsplit(weight_mod_stan_df$condition,".&.")
for (i in 1:length(gen_by_treat_name)){
  gen <- gen_by_treat_name[i][[1]][1]
  gen <- strsplit(gen, "treatment-")[[1]][2]
  treat <- gen_by_treat_name[i][[1]][2]
  weight_mod_stan_df$genotype[i] <- gen
  weight_mod_stan_df$treatment[i] <- treat
}
weight_mod_stan_df <- weight_mod_stan_df[!is.na(weight_mod_stan_df$treatment),]
weight_mod_stan_df <- weight_mod_stan_df[!is.na(weight_mod_stan_df$genotype),]

weight_mod_stan_df <- weight_mod_stan_df[weight_mod_stan_df$treatment!="Control",]

# Plotting the effect of each treatment, all in comparison to Control.
# At the end i am producing 95% CI for all effects.

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_genotype_by_treatments/supplementary_weight_all_treatments_by_genotypes_bayesian_factor_model.pdf", width = 7,height = 10,useDingbats = F)
ggplot(data = weight_mod_stan_df, aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(aes(group=treatment, color=treatment),position=position_dodge(width=1)) +
  theme_bw() +
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  facet_grid(treatment ~ .) 
#! dev.off()





#####Fig S1A - representative exp for weight of boiled pathogens VS. pathogens (Aug results) 
#subsetting to local only and control, #1 and boiled_#1 - here i cleaned outliers, and have no version without outliers cleaning
dataset_Aug_local <- dataset_Aug[dataset_Aug$genotype%in%c("Lu","Eyach","HE","Kus","Tu-Wal","Schl") & dataset_Aug$treatment%in%c("Control","#1","#1_boil"),] 
dataset_Aug_local$treatment <- droplevels.factor(dataset_Aug_local$treatment) # dropping unwanted levels (#2, #3, etc...)
dataset_Aug_local$genotype <- droplevels.factor(dataset_Aug_local$genotype) # dropping unwanted levels (Col0, etc...)
levels(dataset_Aug_local$genotype)<- c("Ey15-2","Kus3-1","Schl-7","Lu3-30","HE-1", "Tu-Wal-2") #renaming host genotypes to original 1001 names

outliers_list_Aug <- outliers_removal(dataset = dataset_Aug_local) # removing outliers
length(outliers_list_Aug)/length(dataset_Aug_local$weight) # Sanity check: the % of outliers removed. about 1%, as expected (2.5 std are 99% of the normal distrubiton)
dataset_Aug_local_no_out <- dataset_Aug_local[-which(row.names(dataset_Aug_local) %in% outliers_list_Aug),]

###TODO: print the plot (dabest) for boiled #1 com as supplementary

#ploting using dabestr - Baesyian based approach

index_dabest <- list() # creating an index for the dabest function (what is control, and to seperate by genotype)
i = 1 # helping index for the loop that will come
sorted_genotypes <- c("Lu3-30","Kus3-1","Schl-7","HE-1", "Tu-Wal-2","Ey15-2") # sorting the genotypes for the plot. I used Tu-Wal before the end, and Eyach at the end to mark its difference 

for (gen in sorted_genotypes){
  current_genotype <- paste(gen, c("Control","#1","#1_boil")) # hard coding - if i change the treatments name - i must change it here also. Did it to sort in the right way
  index_dabest[[i]] <- current_genotype
  i=i+1
}

dataset_Aug_local_no_out$gen_by_treat <- paste(dataset_Aug_local_no_out$genotype,dataset_Aug_local_no_out$treatment) #without outliers

unpaired_mean_diff <- dabest(dataset_Aug_local_no_out, gen_by_treat, weight, paired = F, idx = index_dabest, seed = 12345, ci = 95)
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/figure_S1/figureS1_A_weight.pdf", width = 14, height = 6)
plot(unpaired_mean_diff, color.column = treatment)
#! dev.off()



unpaired_mean_diff_only_treatment <- dabest(dataset_Aug_local_no_out, treatment, weight, paired = F, idx = c("Control","#1","#1_boil"), seed = 12345, ci = 95)


boiled_fit_anova <- anova(lm(formula = weight ~ treatment + genotype,data = dataset_Aug_local_no_out))
boiled_fit_anova$`Sum Sq`

summary(lm(formula = weight ~ treatment,data = dataset_Aug_local_no_out))

boiled_fit <- stan_glm(formula = weight ~ treatment, data = dataset_Aug_local_no_out)
summary(boiled_fit )
boiled_fit$coefficients

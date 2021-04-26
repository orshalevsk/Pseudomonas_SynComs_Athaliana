library(ggplot2)
library(dabestr)
library(vegan)
library(rstan)
library(rstanarm)

### March exp
weight_results_mixed_Eyach_March <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/Eyach_mixed-6_vs_mixed_March_2019/weight_results_20190409.csv")

weight_results_mixed_Eyach_March$weight <- as.numeric(as.character(weight_results_mixed_Eyach_March$weight))
weight_results_mixed_Eyach_March <- weight_results_mixed_Eyach_March[weight_results_mixed_Eyach_March$treatment%in%c("Control","#3","#3-S6"),]
levels(weight_results_mixed_Eyach_March$treatment) <- c("Control","Mixed dP6","Mixed")
weight_results_mixed_Eyach_March$exp <- "March"

### July exp
weight_results_mixed_Eyach_July <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/Eyach_mixed-6_vs_mixed_July_2019/results/Eyach_mixed_strain6_July_2019_weight_results.csv")
weight_results_mixed_Eyach_July<- weight_results_mixed_Eyach_July[!is.na(weight_results_mixed_Eyach_July$Tray),]
weight_results_mixed_Eyach_July$weight<- as.numeric(as.character(weight_results_mixed_Eyach_July$weight))
levels(weight_results_mixed_Eyach_July$Treatment) <- c("", "Mixed", "Mixed_unknown", "Mixed dP1", "Mixed dP6", "Control")
colnames(weight_results_mixed_Eyach_July)[colnames(weight_results_mixed_Eyach_July)=="Treatment"] <- "treatment"
weight_results_mixed_Eyach_July$exp <- "July"
weight_results_mixed_Eyach_July <- weight_results_mixed_Eyach_July[,c("weight","position","treatment","exp")]


weight_results_mixed_dropout <- rbind(weight_results_mixed_Eyach_July,weight_results_mixed_Eyach_March)
#subset only to factor levels that are shared between both experiments
weight_results_mixed_dropout <- weight_results_mixed_dropout[weight_results_mixed_dropout$treatment%in% c("Control","Mixed","Mixed dP6"),]

weight_results_mixed_dropout_no_out <- weight_results_mixed_dropout[-which(row.names(weight_results_mixed_dropout) %in% outliers_list),]




unpaired_mean_diff <- dabest(weight_results_mixed_dropout_no_out[weight_results_mixed_dropout_no_out$exp=="March",], treatment, weight,
                             idx = list(c("Control","Mixed","Mixed dP6")),
                             paired = FALSE)
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_mixed_dropout_exp/weight_Eyach_mixed_dropout_March_exp.pdf")
plot(unpaired_mean_diff,color.column = treatment)
#! dev.off()

weight_results_mixed_dropout_no_out$treatment <- factor(weight_results_mixed_dropout_no_out$treatment, levels = c("Control","Mixed","Mixed dP6"))
weight_mixed_dropout_stan_mod <- stan_glm(formula = weight ~ exp*treatment, data = weight_results_mixed_dropout_no_out)
weight_mixed_dropout_stan_mod$stan_summary


# Producing 95% CI for all effects.

coeffcients_interactions <- c("treatmentMixed",	"treatmentMixed dP6",	"expMarch:treatmentMixed",	"expMarch:treatmentMixed dP6")
coeffcients_interactions_summary <- weight_mixed_dropout_stan_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")]
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)  

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_mixed_dropout_exp/weight_Eyach_mixed_dropout_March_exp_stan_validation.pdf")
ggplot(data = coeffcients_interactions_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()

# Here the PathoCom and PathoCom d6 (April 2021)
weight_results_pathocom_d6_Eyach <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/RE_2021_RNAseq_and_dP6/PathoCom_d6_ranodmized_results.csv")
weight_results_pathocom_d6_Eyach$weight <- as.numeric(as.character(weight_results_pathocom_d6_Eyach$weight))

weight_results_mixed_dropout_no_out <- weight_results_mixed_dropout[-which(row.names(weight_results_mixed_dropout) %in% outliers_list),]

levels(weight_results_pathocom_d6_Eyach$Treatment) <- gsub(pattern = " ", replacement = "", x = levels(weight_results_pathocom_d6_Eyach$Treatment))


unpaired_mean_diff <- dabestr::dabest(weight_results_pathocom_d6_Eyach, Treatment, weight,
                             idx = list(c("control","PathoCom","PathoCom_d6")),
                             paired = FALSE)
difference <- mean_diff(unpaired_mean_diff)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_mixed_dropout_exp/weight_Eyach_PathoCom_p6_dropout_March_exp.pdf", useDingbats = F)
plot(difference,color.column = Treatment)
#! dev.off()

weight_results_mixed_dropout_no_out$treatment <- factor(weight_results_mixed_dropout_no_out$treatment, levels = c("Control","Mixed","Mixed dP6"))
weight_mixed_dropout_stan_mod <- stan_glm(formula = weight ~ exp*treatment, data = weight_results_mixed_dropout_no_out)
weight_mixed_dropout_stan_mod$stan_summary


# Producing 95% CI for all effects. Focusing on Eyach - #3. Again - this is accross both exp, and accounting for exp, pot and tray as random effects

coeffcients_interactions <- c("treatmentMixed",	"treatmentMixed dP6",	"expMarch:treatmentMixed",	"expMarch:treatmentMixed dP6")
coeffcients_interactions_summary <- weight_mixed_dropout_stan_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")]
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)  

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/weight/weight_mixed_dropout_exp/weight_Eyach_mixed_dropout_March_exp_stan_validation.pdf")
ggplot(data = coeffcients_interactions_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()




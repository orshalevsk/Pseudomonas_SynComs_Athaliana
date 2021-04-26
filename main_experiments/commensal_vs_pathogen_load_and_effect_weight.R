library(ggplot2)
library(vegan)
library(lme4)
library(rstanarm)
library(pairwiseAdonis)
library(Hmisc)
library(dplyr)
library(RColorBrewer)
library(randomForest)
library(ggpubr)
library(loo)
library(BayesFactor)

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
full_dataset_local$genotype <- factor(full_dataset_local$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))


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


#preparing strain-by-strain dataframe
source("~/ownCloud/My papers/Syncoms_paper/scripts/dataframe_to_strains_df_convertor.R")
dataset_by_strains <- dataframe_by_strain(dataframe = full_dataset_local_no_out)
dataset_by_strains <- dataset_by_strains[dataset_by_strains$strain_num%in%c(1:7) & dataset_by_strains$treatment%in%c("#1","#3") | dataset_by_strains$strain_num%in%c(11:17) & dataset_by_strains$treatment%in%c("#2","#3"),]



# removing zeros. These cannot be log() transformed
for (i in colnames(dataset_by_strains)){
  zero_row <- which(dataset_by_strains[,i]==0)
  if (length(zero_row)>0){
    dataset_by_strains <- dataset_by_strains[-zero_row,]
  }
}
#### Modeling weight by bacterial load and interactions with treatment - to understand wether the slope of bacterial load is different and depednds on the treatment type.
#### Seems like it doesn't. Though the slope is the same  
#### the intercept is different. I.e. for the same load, different treatments show different average weight. For pathogens, the same load shows bigger reduction in weight.
#### The conclusion is that for the same amount of cells of pathogens and commensals, pathogens will reduce the weight more.

# First i will start with a simple lm. This is the simplest model, and is following the proporionist approach (not Baesiyan yet)
weight_by_load_mod_lm <- lm(formula = weight ~ log(bacterial_load)*treatment + genotype + exp, 
                            data = full_dataset_local_no_out) # examining interactions of load with treatment

summary(weight_by_load_mod_lm) # no significancy for interactions. Bacterial load is a good predictor of weight. Model is ok (R^2=62.6)
#? The intersect of different treatments is different ? (treatment is significant here)
plot(weight_by_load_mod_lm) # sanity check is the model ok?
hist(residuals(weight_by_load_mod_lm)) # sanity check is the model ok? - normal dist of residuals? Approximatley yes


weight_by_load_by_gen_mod_lm <- lm(formula = weight ~ log(bacterial_load)*genotype + treatment + exp,  data = full_dataset_local_no_out) # examining interactions of load with genotype


summary(weight_by_load_by_gen_mod_lm) #no significancy for interactions. Bacterial load is a good predictor of weight.Model is ok (R^2=62.7)
                                        #? The intersect of different treatments is different ? (treatment is significant here)
plot(weight_by_load_by_gen_mod_lm) # sanity check is the model ok?
hist(residuals(weight_by_load_by_gen_mod_lm)) # sanity check is the model ok? - normal dist of residuals? Approximatley yes


weight_by_load_by_gen_mod_lm_com3 <- lm(formula = weight ~ log(bacterial_load)*genotype + exp, 
                                   data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2",]) # examining interactions of load with genotype

weight_by_load_by_gen_mod_lm_com3 <- lm(formula = weight ~ log(bacterial_load)*genotype + exp, 
                                        data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",]) # examining interactions of load with genotype


summary(weight_by_load_by_gen_mod_lm_com3) #no significancy for interactions. Bacterial load is a good predictor of weight.Model is ok (R^2=62.7)
#? The intersect of different treatments is different ? (treatment is significant here)

# Here are attempts to try and fit an lm to the data

#Now modeling the same, but taking the Baesiyan approach


weight_by_load_mod <- lmBF(formula = weight ~ log_bacterial_load*treatment*genotype, data = full_dataset_local_no_out,posterior = T,iter=10000)
chains_sum <- as.data.frame(summary(weight_by_load_mod)[[2]])

#weight_by_load_mod <- stan_glm(formula = weight ~ log(bacterial_load)*treatment*genotype + exp, data = full_dataset_local_no_out)


weight_by_load <- aov(formula = weight ~ log_bacterial_load*treatment  , data = full_dataset_local_no_out)
weight_by_load <- aov(formula = log_bacterial_load ~ treatment*exp  , data = full_dataset_local_no_out[!full_dataset_local_no_out$treatment=="#1",])

ggplot(data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",], aes(x=treatment, y=log_bacterial_load))+
  geom_violin()
summary(weight_by_load)

plot(weight_by_load_mod)

coeffcients_interactions <- c("log(bacterial_load):treatment#2", "log(bacterial_load):treatment#3", "log(bacterial_load)", "treatment#2", "treatment#3")
coeffcients_interactions_summary <- weight_by_load_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_model_stan_validation.pdf", width = 7,height = 7)
ggplot(data = coeffcients_interactions_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()


###NOW LETS PLOT THE RAW DATA THAT THE MODEL SUPPORTED
##Plotting the regression of bacterial load to weight, to examplify the models and the relationship of load-weight.
#both experiments
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_both_exp.pdf", useDingbats = F)
ggplot(aes(x=log(bacterial_load),y=weight, color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#1","#2"),])+
  geom_point(aes(color = treatment))+
  #geom_smooth(method='glm') +
  geom_smooth(method = 'loess') +
  #xlim(c(-7.5,2)) +
  theme_bw()
  #facet_grid(genotype ~ .)
#dev.off()


weight_by_load_with_treatment <- stan_glm(formula = weight ~ log_bacterial_load*treatment + genotype + exp, data = full_dataset_local_no_out)
weight_by_load_without_treatment <- stan_glm(formula = weight ~ log_bacterial_load + genotype + exp, data = full_dataset_local_no_out)

weight_by_load_with_treatment_lm <- lm(formula = weight ~ log_bacterial_load*treatment + genotype + exp, data = full_dataset_local_no_out)
weight_by_load_without_treatment_lm <- lm(formula = weight ~ log_bacterial_load + genotype + exp, data = full_dataset_local_no_out)




loo_compare(loo(weight_by_load_with_treatment),loo(weight_by_load_without_treatment))

library(rstan)
stan_plot(weight_by_load_mod, point_est = "mean", show_density = TRUE, fill_color = "maroon")
plot(weight_by_load_mod)
weight_by_load_mod$stanfit

full_dataset_local_no_out$treatment <- factor(full_dataset_local_no_out$treatment)

weight_by_load_mod <- brms::brm(formula = weight ~ log_bacterial_load*treatment + genotype + exp, data = full_dataset_local_no_out)
library("modelr")
library("tidybayes")

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_both_exp_with_Baesiyan_confidence.pdf", useDingbats = F)
full_dataset_local_no_out %>%
  group_by(treatment) %>%
  #data_grid(log_bacterial_load = seq_range(log_bacterial_load, n = 51)) %>%
  add_fitted_draws( n = 100,weight_by_load_mod) %>%
  ggplot(aes(x = log_bacterial_load, y = weight, color = ordered(treatment))) +
  stat_lineribbon(aes(y = .value)) +
  geom_point() +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  theme_bw()
#! dev.off()



new_x<-data.frame(x1=new_X[,2],x2=rep(c("Min","Mean","Max"),each=20))
new_y<-extract(m_norm,pars="y_pred")
pred<-apply(new_y[[1]],2,quantile,probs=c(0.025,0.5,0.975)) #the median line with 95% credible intervals
#plot
plot(dat$x1,y_norm,pch=16)
lines(new_x$x1[1:20],pred[2,1:20],col="red",lwd=3)
lines(new_x$x1[1:20],pred[2,21:40],col="orange",lwd=3)
lines(new_x$x1[1:20],pred[2,41:60],col="blue",lwd=3)
lines(new_x$x1[1:20],pred[1,1:20],col="red",lwd=1,lty=2)
lines(new_x$x1[1:20],pred[1,21:40],col="orange",lwd=1,lty=2)
lines(new_x$x1[1:20],pred[1,41:60],col="blue",lwd=1,lty=2)
lines(new_x$x1[1:20],pred[3,1:20],col="red",lwd=1,lty=2)
lines(new_x$x1[1:20],pred[3,21:40],col="orange",lwd=1,lty=2)
lines(new_x$x1[1:20],pred[3,41:60],col="blue",lwd=1,lty=2)
legend("topright",legend=c("Min","Mean","Max"),lty=1,col=c("red","orange","blue"),bty = "n",title = "Effect of x2 value on\nthe regression")

#Now the same trick, only treating pathogens and commensals seperately (load of com3_patho #3 vs #1, and com3_comm #3 s #2)
weight_by_load_mod_patho <- stan_glm(formula = weight ~ log(com3_patho)*treatment + genotype + exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#1","#3"),])
coeffcients_interactions_patho <- c("log(com3_patho):treatment#3", "log(com3_patho)", "treatment#3")
coeffcients_interactions_summary_patho <- weight_by_load_mod_patho$stan_summary[coeffcients_interactions_patho,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_interactions_summary_patho <- as.data.frame(coeffcients_interactions_summary_patho)
coeffcients_interactions_summary_patho$condition <- rownames(coeffcients_interactions_summary_patho)

weight_by_load_mod_commensal<- stan_glm(formula = weight ~ log(com3_commensal)*treatment + genotype + exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#2","#3"),])
coeffcients_interactions_commensal <- c("log(com3_commensal):treatment#3", "log(com3_commensal)", "treatment#3")
coeffcients_interactions_summary_commensal <- weight_by_load_mod_commensal$stan_summary[coeffcients_interactions_commensal,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_interactions_summary_commensal <- as.data.frame(coeffcients_interactions_summary_commensal)
coeffcients_interactions_summary_commensal$condition <- rownames(coeffcients_interactions_summary_commensal)


coeffcients_interactions_summary_both_commensal_patho <- rbind(coeffcients_interactions_summary_patho,coeffcients_interactions_summary_commensal)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_model_stan_validation.pdf", width = 7,height = 7)
ggplot(data = coeffcients_interactions_summary_both_commensal_patho, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()


###NOW LETS PLOT THE RAW DATA THAT THE MODEL SUPPORTED
###Plotting the regression of load per type (pathogen vs. commensal) to weight, to examplify the models and the relationship of load-weight per type of strain.
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_only_patho_both_exp.pdf")
ggplot(aes(x=log(com3_patho),y=weight, color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#1","#3"),])+
  geom_point(aes(color = treatment))+
  geom_smooth(method='lm') +
  theme_bw() 
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_only_commensal_both_exp.pdf")
ggplot(aes(x=log(com3_commensal),y=weight, color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in%c("#2","#3"),])+
  geom_point(aes(color = treatment))+
  geom_smooth(method='lm') +
  theme_bw() 
#! dev.off()

######################





for (gen in factor(full_dataset_local$genotype)){
  med_com2 <- median(full_dataset_local_no_out$weight[full_dataset_local_no_out$genotype==gen & full_dataset_local_no_out$treatment=="#2"])
  full_dataset_local_no_out$weight_div_by_med_com2[full_dataset_local_no_out$genotype==gen]<- full_dataset_local_no_out$weight[full_dataset_local_no_out$genotype==gen]/med_com2
}
# Here i moved to the segmented approach

#full_dataset_local_no_out$weight<- full_dataset_local_no_out$weight_div_by_med_com2
library(segmented)


hist(full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out_original$treatment=="#2"],breaks = 100)
full_dataset_local_no_out_original <- full_dataset_local_no_out




#full_dataset_local_no_out <- full_dataset_local_no_out[full_dataset_local_no_out$log_bacterial_load>-5,]
#filter_outlier_com2 <- !(full_dataset_local_no_out$treatment=="#2" & full_dataset_local_no_out$log_bacterial_load>0.5)
#full_dataset_local_no_out<-full_dataset_local_no_out[filter_outlier_com2,]
#full_dataset_local_no_out <- full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",]
full_dataset_local_no_out <- full_dataset_local_no_out[full_dataset_local_no_out$genotype!="Tu-Wal-2",]


my.lm<- glm(formula = weight ~ log_bacterial_load, data = full_dataset_local_no_out)

my.lm_com1<- glm(formula = weight ~ log_bacterial_load , data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",])
my.lm_com2<- glm(formula = weight ~ log_bacterial_load , data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2",])
my.lm_com3<- glm(formula = weight ~ log_bacterial_load , data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",])



# have to provide estimates for breakpoints.
# after looking a the data, 
my.seg_com1 <- segmented(my.lm_com1, 
                         seg.Z = ~ log_bacterial_load, 
                         psi = list(log_bacterial_load = c(-3,-1.5)))

my.seg_com2 <- segmented(my.lm_com2, 
                         seg.Z = ~ log_bacterial_load, 
                         psi = list(log_bacterial_load = c(-3,-1.5)))

my.seg_com3 <- segmented(my.lm_com3, 
                         seg.Z = ~ log_bacterial_load, 
                         psi = list(log_bacterial_load = c(-3,-1.5)))

davies.test(obj = my.lm,seg.Z = ~ log_bacterial_load,k = 3,type = "wald")
pscore.test(obj = my.lm,seg.Z = ~ log_bacterial_load,k = 3)

my.seg_com1$psi
slope(my.seg_com1)

my.fitted_com1 <- fitted(my.seg_com1)
my.fitted_com2 <- fitted(my.seg_com2)
my.fitted_com3 <- fitted(my.seg_com3)

my.model_com1 <- data.frame(bacterial_load = full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out$treatment=="#1"], weight = my.fitted_com1)
my.model_com2 <- data.frame(bacterial_load = full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out$treatment=="#2"], weight = my.fitted_com2)
my.model_com3 <- data.frame(bacterial_load = full_dataset_local_no_out$log_bacterial_load[full_dataset_local_no_out$treatment=="#3"], weight = my.fitted_com3)

my.model_com1$treatment <- "#1"
my.model_com2$treatment <- "#2"
my.model_com3$treatment <- "#3"

my.model <- rbind(my.model_com1,my.model_com2,my.model_com3)

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_both_exp_segmented_version.pdf",useDingbats = F)
ggplot(my.model, aes(x = bacterial_load, y = weight, color=treatment)) +
  geom_line() +
  #geom_smooth(data = full_dataset_local_no_out, stat='smooth', method = "lm", se=TRUE, alpha=0.07, 
  #            aes(color = NULL,fill=treatment, colour = factor(treatment),x=log_bacterial_load)) +
  theme_bw() +
  geom_point(data = full_dataset_local_no_out, aes(x=log_bacterial_load),alpha = 0.5,stroke=0.2 ,size=2)+
  geom_ribbon(data = full_dataset_local_no_out, stat='smooth', method = "loess", se=TRUE, alpha=0.07, 
              aes(color = NULL,fill=treatment, group = factor(treatment),x=log_bacterial_load)) +
  #geom_line(stat='smooth', method = "loess", alpha=0.5,linetype = "dashed") +
  scale_fill_manual(values = c("#d95f02",	"#1b9e77",	"#7570b3")) +
  scale_color_manual(values = c("#d95f02",	"#1b9e77",	"#7570b3")) +
  ylim(c(0,310))
#dev.off()




full_dataset_local_no_out <- full_dataset_local_no_out_original


#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_both_exp_no_Tu-Wal.pdf",useDingbats = F)
ggplot(full_dataset_local_no_out[full_dataset_local_no_out$genotype!="Lu3-30",], aes(x = log_bacterial_load, y = weight, color=treatment)) +
  #geom_line() +
  geom_smooth(method = "lm",alpha=0.07, aes(fill=treatment, group = factor(treatment)))+
  theme_bw() +
  geom_point(alpha = 0.5,stroke=0.2 ,size=2)+
  #geom_ribbon(data = full_dataset_local_no_out, stat='smooth', method = "lm", se=TRUE, alpha=0.07, 
  #            aes(color = NULL,fill=treatment, group = factor(treatment),x=log_bacterial_load)) +
  scale_fill_manual(values = c("#d95f02",	"#1b9e77",	"#7570b3")) +
  scale_color_manual(values = c("#d95f02",	"#1b9e77",	"#7570b3"))
  #facet_grid(genotype ~ .,scales = 'free')
#dev.off()


weight_by_load_mod <- lmBF(formula = weight ~ log_bacterial_load*treatment + genotype, data = full_dataset_local_no_out,posterior = T,iter=10000)
chains_sum <- as.data.frame(summary(weight_by_load_mod)[[2]])

weight_by_load_mod <- stan_glm(formula = weight ~ log_bacterial_load*treatment + genotype, data = full_dataset_local_no_out)
weight_by_load_mod$stan_summary


coeffcients_interactions <- c("log_bacterial_load:treatment#2", "log_bacterial_load:treatment#3")
coeffcients_interactions_summary <- weight_by_load_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)

#coeffcients_interactions_summary[,c("2.5%","50%","97.5%")] <- coeffcients_interactions_summary[,c("2.5%","50%","97.5%")]+2 # adding a factor due to the model
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_with_genotype_model_stan_validation_MODIFIED.pdf", useDingbats = F)
ggplot(data = coeffcients_interactions_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color=condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  #geom_text(y=10, x="log_bacterial_load:treatment#3" ,label="weight ~ log_bacterial_load*treatment + genotype")+
  ylim(c(-10,23)) +
  scale_color_manual(values = c("#1b9e77" ,"#7570b3"))+
  theme_bw()
#dev.off()


### now i will plot the intercept difference between the treatment
coeffcients_interactions <- c("treatment#2", "treatment#3")
coeffcients_interactions_summary <- weight_by_load_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_with_genotype_model_stan_validation_intercept.pdf", useDingbats = F)
ggplot(data = coeffcients_interactions_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_text(y=10, x="treatment#3" ,label="weight ~ log_bacterial_load*treatment + genotype")+
  theme_bw()
#dev.off()



#fit <- glm(formula = weight ~ log_bacterial_load*treatment, data = full_dataset_local_no_out[full_dataset_local_no_out$genotype!="Tu-Wal-2",])
#library(interactions)
#library(ggstance)

#sim_slopes_bac_load <- sim_slopes(fit, pred = log_bacterial_load, modx = treatment, johnson_neyman = F)
#plot(sim_slopes_bac_load)
#sim_slopes_bac_load$slopes
##### now i will check if every strain behaves the same in different treamtnets - hence it is only a matter of load
stan_summary_load_by_treat <- data.frame("strain"=NA, "mod_2.5" = NA,"mod_50" = NA, "mod_97.5" = NA)

dataset_by_strains$genotype <- factor(dataset_by_strains$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2"))

# removing zeros. These cannot be log() transformed
for (i in colnames(dataset_by_strains)){
  zero_row <- which(dataset_by_strains[,i]==0)
  if (length(zero_row)>0){
    dataset_by_strains <- dataset_by_strains[-zero_row,]
  }
}

chosen_dataset <- dataset_by_strains[dataset_by_strains$exp=="Oct",] #choose data set to work with (subset to one exp?)
# first for pathogenic strains
coef_tested <- c("treatment#3")
for (strain in c(1:7)){
  strain_abundance_by_treat_mod <- stan_glm(weight ~ treatment*log(strain_abundance), data = chosen_dataset[chosen_dataset$treatment%in% c("#1","#3") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary)))
}
# now for commensal strains
for (strain in c(11:17)){
  strain_abundance_by_treat_mod <- stan_glm(weight ~ treatment*log(strain_abundance), data = chosen_dataset[chosen_dataset$treatment%in% c("#2","#3") & chosen_dataset$strain_num %in% c(strain),], seed = 12345)
  strain_abundance_summary <- strain_abundance_by_treat_mod$stan_summary[coef_tested,c("2.5%","50%","97.5%")]
  strain_abundance_summary <- as.data.frame(t(strain_abundance_summary))
  colnames(strain_abundance_summary) <- c("mod_2.5", "mod_50", "mod_97.5")
  stan_summary_load_by_treat <- rbind(stan_summary_load_by_treat,cbind("strain"=strain,rbind(strain_abundance_summary)))
}

stan_summary_load_by_treat$strain <- factor(stan_summary_load_by_treat$strain)
stan_summary_load_by_treat <- stan_summary_load_by_treat[!is.na(stan_summary_load_by_treat$strain),]
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(1:7)] <- "OTU5"
stan_summary_load_by_treat$strain_type[stan_summary_load_by_treat$strain %in% c(11:17)] <- "non_OTU5"
stan_summary_load_by_treat$strain_type <- factor(stan_summary_load_by_treat$strain_type)


# now plotting the stan model results, taking into account both experiments and this time also the genotype
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/all_strains_change_com3_by_genotype.pdf",width = 7, height = 14)
ggplot(aes(x=strain,y=stan_summary_load_by_treat$mod_50, ymin = stan_summary_load_by_treat$mod_2.5, ymax = stan_summary_load_by_treat$mod_97.5),
       data = stan_summary_load_by_treat)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()

# now plotting the stan model results, taking into account both experiments and this time also the genotype - only pathogens
stan_summary_load_by_treat_patho <- stan_summary_load_by_treat[stan_summary_load_by_treat$strain%in%c(1:7),]
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/microbe-microbe_interactions/Between_treatments/patho_strains_change_com3_by_genotype.pdf",width = 7, height = 9)
ggplot(aes(x=treatment,y=stan_summary_load_by_treat_patho$mod_50, ymin = stan_summary_load_by_treat_patho$mod_2.5, ymax = stan_summary_load_by_treat_patho$mod_97.5),
       data = stan_summary_load_by_treat_patho)+
  geom_pointrange(aes(color = strain_type),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(strain ~ .)
#! dev.off()





##Plotting the regression of bacterial load to weight, to examplify the models and the relationship of load-weight, genotype-wise!
ggplot(aes(x=log(bacterial_load),y=weight, color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$treatment %in% c("#1"),])+
  geom_point(aes(color =genotype))+
  geom_smooth(method='lm') +
  theme_bw()

##Plotting the regression of bacterial load to weight, to examplify the models and the relationship of load-weight, genotype-wise!
ggplot(aes(x=log(X5),y=weight, color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$treatment %in% c("#1", "#3"),])+
  geom_point(aes(color =treatment))+
  geom_smooth(method='lm') +
  theme_bw()



#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_treatment_association_both_exp.pdf")
ggplot(aes(x=log(bacterial_load),y=weight, color = treatment), data = full_dataset_local_no_out)+
  geom_point(aes(color = treatment))+
  geom_smooth(method='lm') +
  theme_bw() +
  facet_grid(genotype ~ .)
#! dev.off()

## strain-wise
ggplot(aes(x=log(strain_abundance),y=weight, color = strain_num), data = dataset_by_strains[dataset_by_strains$treatment!="#2" & dataset_by_strains$strain_num%in% c(1:7),])+
  geom_point(aes(color = strain_num))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(strain_num ~ .)



'''
Decided to reomve this part for now... But essentially it looks for host*microbe-microbe interactions

## But it may be that the genotype explains most of the correlations (different genotypes have different loads - as will be detailed later in this script)
## Hence, does different host genotypes show different microbe-microbe interations?

weight_by_load_mod_lm <- lm(formula = cbind(log(X15)) ~ log(X17)*genotype, 
                            data = full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$exp=="Aug",]) # examining interactions of load with treatment


weight_by_load_mod_lm <- lm(formula = cbind(log(X17)) ~ log(X1)*genotype+ log(X2)*genotype+ log(X3)*genotype+ log(X4)*genotype + log(X5)*genotype + log (X6)*genotype + log(X7)*genotype + log(X11)*genotype + log(X12)*genotype + log(X13)*genotype + log(X14)*genotype + 
                              log(X15)*genotype + log(X16)*genotype + exp, 
                            data = full_dataset_local_no_out_strains_only)

summary(weight_by_load_mod_lm)


correlation_fun <- function(dataset)
{
  dataset_mat <- as.matrix(dataset[,paste("X",c(1:7,11:17),sep = "")])
  return(rcorr(log(dataset_mat)))
}

corr_dataset <- full_dataset_local_no_out_strains_only[full_dataset_local_no_out_strains_only$treatment=="#3",c("exp","genotype",paste("X",c(1:7,11:17),sep = ""))]
correlation_fun(dataset = corr_dataset[corr_dataset$exp=="Oct",])
correlation_fun(dataset = corr_dataset[corr_dataset$exp=="Aug",])
'''


## Now i will try to understand the difference in load between different genotypes, by treatment or by pathogens / commensals in com3
# first by Baesiyan model
load_by_treatment_mod <- stan_glm(formula = log(bacterial_load) ~ treatment, data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug",]) # regardless of genotype

coeffcients_interactions <- c("treatment#2", "treatment#3")
coeffcients_interactions_summary <- load_by_treatment_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")]
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$genotype <- rownames(coeffcients_interactions_summary)  

load_by_treatment_by_gen_mod <- stan_glm(formula = log(bacterial_load) ~ treatment*genotype, data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug",]) # genotype-wise

coeffcients_interactions_gen <- c("treatment#2:genotypeHE-1",	"treatment#3:genotypeHE-1",	"treatment#2:genotypeEy15-2",
                                  "treatment#3:genotypeEy15-2",	"treatment#2:genotypeLu3-30",	"treatment#3:genotypeLu3-30",
                                  "treatment#2:genotypeSchl-7",	"treatment#3:genotypeSchl-7",	"treatment#2:genotypeTu-Wal-2",
                                  "treatment#3:genotypeTu-Wal-2",	"treatment#3:genotypeTu-Wal-2")
coeffcients_interactions_gen_summary <- load_by_treatment_by_gen_mod$stan_summary[coeffcients_interactions_gen,c("2.5%","50%","97.5%")]
coeffcients_interactions_gen_summary <- as.data.frame(coeffcients_interactions_gen_summary)
coeffcients_interactions_gen_summary$genotype <- rownames(coeffcients_interactions_gen_summary)  

# load differences between genotypes infected with PathoCom
load_PathoCom_by_gen_mod <- stan_glm(formula = log(bacterial_load) ~ genotype*exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",]) # genotype-wise

load_PathoCom_by_gen_mod_glm <- glm(formula = log(X5) ~ genotype*exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",]) # genotype-wise
summary(load_PathoCom_by_gen_mod_glm)

coeffcients_interactions_gen <- c("genotypeEy15-2",	"genotypeSchl-7",	"genotypeLu3-30",
                                  "genotypeHE-1",	"genotypeTu-Wal-2")
coeffcients_interactions_gen_summary <- load_PathoCom_by_gen_mod $stan_summary[coeffcients_interactions_gen,c("2.5%","50%","97.5%")]
coeffcients_interactions_gen_summary <- as.data.frame(coeffcients_interactions_gen_summary)
coeffcients_interactions_gen_summary$genotype <- rownames(coeffcients_interactions_gen_summary)  

# plot the PathoCom load by genotype, stan model results

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_genotype_PathoCom_stan_validation.pdf.pdf",width = 7, height = 9)
ggplot(aes(x=genotype,y=coeffcients_interactions_gen_summary$`50%`, ymin = coeffcients_interactions_gen_summary$`2.5%`, ymax = coeffcients_interactions_gen_summary$`97.5%`),
       data = coeffcients_interactions_gen_summary)+
  geom_pointrange(aes(color = genotype),position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()



#regardless of the genotype - density plot
ggplot(aes(x=log(bacterial_load), color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Oct",])+
  geom_density() +
  theme_bw()+
  facet_grid(genotype ~ .) # add for the genotype

#genotype-wise by treatment - density plot
ggplot(aes(x=log(bacterial_load), color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug",])+
  geom_density() +
  theme_bw() +
  facet_grid(treatment ~ .) # add for the genotype


#regardless of the genotype - Baesiyan lm summary
ggplot(aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`), data = coeffcients_interactions_summary)+
  geom_pointrange(aes(group=genotype, color=genotype),position=position_dodge(width=1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  theme_bw()

# the pathogens and commensal proportion within community #3. Need to be done, and redone more systemically
hist(full_dataset_local_no_out$com3_patho[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug"]/full_dataset_local_no_out$bacterial_load[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug"])

## #3 vs #1 loads - only pathogens
#regardless of the genotype - density plot only for pathogens in com #3 vs #1
ggplot(aes(x=log(com3_patho), color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#1","#3"),])+
  geom_density() +
  theme_bw() +
  facet_grid(genotype ~ .) # add for the genotype

#regardless of the genotype - violin plot only for pathogens in com #3 vs #1
ggplot(aes(x = treatment, y=log(com3_patho), fill = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#1","#3"),])+
  geom_violin() +
  theme_bw() +
  facet_grid(genotype ~ .) # add for the genotype

#genotype-wise by treatment - density plot
ggplot(aes(x=log(com3_patho), color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#1","#3"),])+
  geom_density() +
  theme_bw() +
  facet_grid(treatment ~ .) # add for the genotype

## #3 vs #2 loads - only commensals
#regardless of the genotype - density plot only for commensals in com #3 vs #2
ggplot(aes(x=log(com3_commensal), color = treatment), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#2","#3"),])+
  geom_density() +
  theme_bw() +
  facet_grid(genotype ~ .) # add for the genotype

#genotype-wise by treatment - density plot
ggplot(aes(x=log(com3_commensal), color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#2","#3"),])+
  geom_density() +
  theme_bw() +
  facet_grid(treatment ~ .) # add for the genotype


# genotype-wise load of bacterial load in PathoCom, density plot
ggplot(aes(x=log(bacterial_load), color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#1"),])+
  geom_density() +
  theme_bw() + 
  facet_grid(genotype ~ .)

# genotype-wise load of bacterial load in PathoCom, density plot
ggplot(aes(x=log(X6), color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#1"),])+
  geom_density() +
  theme_bw() + 
  facet_grid(genotype ~ .)


# genotype-wise load of bacterial load in PathoCom
ggplot(aes(y=log(X6), x=genotype ,color = genotype), data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment %in% c("#1"),])+
  geom_violin() +
  geom_point() +
  theme_bw() 

######################## now i will bin pathogens and commensals from mixed and exclusive (bin by OTU5 / non-OTU5, and regardless of the treatment)
####################### This follows an analysis i already made, but now i will not take each treatmet, but the load by strain-type.
#build a data frame for pathogens vs commensals in 'Mixed' SynCom, vs the exclusive syncoms
com3_patho_dataframe <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3", c("com3_patho","genotype","treatment","exp","weight")]
com3_patho_dataframe$load_type <- "pathogens_mixed"
names(com3_patho_dataframe)[1] <- "load"

com3_comm_dataframe <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3", c("com3_commensal","genotype","treatment","exp","weight")]
com3_comm_dataframe$load_type <- "commensals_mixed"
names(com3_comm_dataframe)[1] <- "load"

com1_patho_dataframe <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1", c("bacterial_load","genotype","treatment","exp","weight")]
com1_patho_dataframe$load_type <- "pathogens_patho"
names(com1_patho_dataframe)[1] <- "load"

com2_comm_dataframe <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2", c("bacterial_load","genotype","treatment","exp","weight")]
com2_comm_dataframe$load_type <- "commensals_comm"
names(com2_comm_dataframe)[1] <- "load"

com3_load_dataframe <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3", c("bacterial_load","genotype","treatment","exp","weight")]
com3_load_dataframe$load_type <- "total_load_mixed"
names(com3_load_dataframe)[1] <- "load"


patho_vs_commensal <- rbind(com3_comm_dataframe,com3_patho_dataframe,com2_comm_dataframe,com1_patho_dataframe,com3_load_dataframe)
patho_vs_commensal$load_type <- as.factor(patho_vs_commensal$load_type)
patho_vs_commensal$strain_type[patho_vs_commensal$load_type%in% c("commensals_comm","commensals_mixed")] <- "non-OTU5"
patho_vs_commensal$strain_type[patho_vs_commensal$load_type%in% c("pathogens_patho","pathogens_mixed")] <- "OTU5"

#plotting pathogens vs commensals in exclusive and mixed.
ggplot(aes(x=log(load),y=weight, color = load_type), data = patho_vs_commensal[!patho_vs_commensal$load_type=="total_load_mixed",])+
  geom_point(aes(color = load_type))+
  geom_smooth(method='lm') +
  theme_bw()


#plotting pathogens vs commensals correlting with weight, grouped by strain type (pathogen or commensal)
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_strain_type_correlation.pdf")
ggplot(aes(x=log(load),y=weight, color = strain_type), data = patho_vs_commensal[!patho_vs_commensal$load_type=="total_load_mixed",])+
  geom_point(aes(color = strain_type))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(genotype ~ .) # OPTIONAL
#! dev.off()

#plotting pathogens vs commensals correlting with weight, grouped by strain type (pathogen or commensal), Excluding MixedCom
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_strain_type_correlation.pdf")
ggplot(aes(x=log(load),y=weight, color = strain_type), data = patho_vs_commensal[!patho_vs_commensal$load_type%in%c("total_load_mixed","pathogens_mixed","commensals_mixed"),])+
  geom_point(aes(color = strain_type))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(genotype ~ .) # OPTIONAL
#! dev.off()


#plotting pathogens vs commensals correlting with weight, grouped by strain type (pathogen or commensal), only in MixedCom
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_strain_type_correlation.pdf")
ggplot(aes(x=log(load),y=weight, color = strain_type), data = patho_vs_commensal[!patho_vs_commensal$load_type%in%c("total_load_mixed","pathogens_patho","commensals_comm"),])+
  geom_point(aes(color = strain_type))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(genotype ~ .) # OPTIONAL
#! dev.off()


#Now validating the former plor by stan_glm function, and plotting the coeeficients results
# first removing zeros. These cannot be log() transformed
for (i in colnames(dataset_by_strains)){
  zero_row <- which(dataset_by_strains[,i]==0)
  if (length(zero_row)>0){
    dataset_by_strains <- dataset_by_strains[-zero_row,]
  }
}


#patho_vs_commensal$load_min <- log(patho_vs_commensal$load) - min(log(patho_vs_commensal$load)) # should i force the intercept at minimum load??

#weight_by_strain_type_mod <- lmer(formula = weight ~ log(load)*treatment + (1| genotype) + exp , data = patho_vs_commensal)
#weight_by_strain_type_mod1 <- lmer(formula = weight ~ log(load) + (1| genotype) + exp , data = patho_vs_commensal)
#anova(weight_by_strain_type_mod,weight_by_strain_type_mod1)

#### now i will make stran model (lmer) and compare them to understand if the treatment and taxonomy (OTU5 or not), effect the regression slope of weight
weight_by_treatment_mod <- stan_lmer(formula = weight ~ log(load)*treatment + (1| genotype) + exp , data = patho_vs_commensal)
weight_by_treatment_mod1 <- stan_lmer(formula = weight ~ log(load)+treatment + (1| genotype) + exp , data = patho_vs_commensal)
a_treat<- loo(weight_by_treatment_mod) # epld_loo is a metric like AIC and BIC
b_treat<- loo(weight_by_treatment_mod1) # epld_loo is a metric like AIC and BIC
loo::compare(a_treat,b_treat) # epld_loo comparison. The lower elpd_loo is the better model 


weight_by_strain_type_mod <- stan_lmer(formula = weight ~ log(load)*strain_type + (1| genotype) + exp , data = patho_vs_commensal[patho_vs_commensal$treatment=="#3",])
weight_by_strain_type_mod1 <- stan_lmer(formula = weight ~ log(load)+strain_type + (1| genotype) + exp , data = patho_vs_commensal[patho_vs_commensal$treatment=="#3",])
a_strain<- loo(weight_by_strain_type_mod) # epld_loo is a metric like AIC and BIC
b_strain<- loo(weight_by_strain_type_mod1) # epld_loo is a metric like AIC and BIC
loo::compare(a_strain,b_strain) # epld_loo comparison. The lower elpd_loo is the better model 



weight_by_strain_type_mod <- stan_glm(formula = weight ~ log(load)*strain_type + genotype + exp , data = patho_vs_commensal, seed=1234)

coeffcients_weight_by_strain_type <- c("log(load):strain_typeOTU5", "strain_typeOTU5")
coeffcients_weight_by_strain_type_summary <-weight_by_strain_type_mod$stan_summary[coeffcients_weight_by_strain_type,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_weight_by_strain_type_summary <- as.data.frame(coeffcients_weight_by_strain_type_summary)
coeffcients_weight_by_strain_type_summary$condition <- rownames(coeffcients_weight_by_strain_type_summary)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_strain_type_model_stan_validation.pdf", width = 7,height = 7)
ggplot(data = coeffcients_weight_by_strain_type_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()

#plotting genotype-wise, load correlting with weight, grouped by strain type (pathogen or commensal). This is to show genotype tolerance / resistance.
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_strain_type_correlation.pdf")
ggplot(aes(x=log(load),y=weight, color = genotype), data = patho_vs_commensal[!patho_vs_commensal$load_type=="total_load_mixed",])+
  geom_point(aes(color = genotype))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(strain_type ~ .) # OPTIONAL
#! dev.off()

# stan model for regression of load per genotype, within a treatment (demonstrating tolerance)
weight_by_genotype_by_treatment_com1 <- stan_glm(formula = weight ~ log(bacterial_load)*genotype + exp , data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",], seed=1234)
weight_by_genotype_by_treatment_com2 <- stan_glm(formula = weight ~ log(bacterial_load)*genotype + exp , data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2",], seed=1234)
weight_by_genotype_by_treatment_com3 <- stan_glm(formula = weight ~ log(bacterial_load)*genotype + exp , data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",], seed=1234)

#now plotting results
col_names_mod <- c("log(bacterial_load):genotypeEy15-2", "log(bacterial_load):genotypeSchl-7", "log(bacterial_load):genotypeLu3-30", "log(bacterial_load):genotypeHE-1", "log(bacterial_load):genotypeTu-Wal-2")
weight_by_genotype_by_treatment_com1_summary <-weight_by_genotype_by_treatment_com1$stan_summary[col_names_mod,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
weight_by_genotype_by_treatment_com1_summary <- as.data.frame(weight_by_genotype_by_treatment_com1_summary)
weight_by_genotype_by_treatment_com1_summary$treatment <- "PathoCom"

weight_by_genotype_by_treatment_com2_summary <-weight_by_genotype_by_treatment_com2$stan_summary[col_names_mod,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
weight_by_genotype_by_treatment_com2_summary <- as.data.frame(weight_by_genotype_by_treatment_com2_summary)
weight_by_genotype_by_treatment_com2_summary$treatment <- "CommenCom"

weight_by_genotype_by_treatment_com3_summary <-weight_by_genotype_by_treatment_com3$stan_summary[col_names_mod,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
weight_by_genotype_by_treatment_com3_summary <- as.data.frame(weight_by_genotype_by_treatment_com3_summary)
weight_by_genotype_by_treatment_com3_summary$treatment <- "MixedCom"

weight_by_genotype_by_treatment_summary <- rbind(weight_by_genotype_by_treatment_com1_summary, weight_by_genotype_by_treatment_com2_summary, weight_by_genotype_by_treatment_com3_summary)
weight_by_genotype_by_treatment_summary$condition <- row.names(weight_by_genotype_by_treatment_summary)

strsplit(x = weight_by_genotype_by_treatment_summary$condition, split = "genotype")
weight_by_genotype_by_treatment_summary$genotype <- rep(c("Ey15-2","Schl-7","Lu3-30","HE-1","Tu-Wal-2"),3)
weight_by_genotype_by_treatment_summary$treatment <- factor(weight_by_genotype_by_treatment_summary$treatment,levels = c("PathoCom","CommenCom","MixedCom"))

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_treatment_model_stan_validation.pdf", width = 7,height = 7, useDingbats = F)
ggplot(data = weight_by_genotype_by_treatment_summary, aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`,color=treatment))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3")) +
  facet_grid(treatment ~ .)
#! dev.off()




library(RColorBrewer)
display.brewer.all()
#plotting genotype-wise, load correlting with weight, grouped by treatment. This is to show genotype tolerance / resistance.
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_treatment_correlation.pdf", useDingbats = F)
ggplot(aes(x=log(bacterial_load),y=weight, color = genotype), data = full_dataset_local_no_out)+
  geom_point(aes(color = genotype))+
  geom_smooth(method='lm', alpha=0.1,aes(fill = genotype)) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2")+
  facet_grid(treatment ~ .) # OPTIONAL
#! dev.off()


##plotting strain-wise, within exclusive communities
# first PathoCom
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_strain_type_correlation.pdf")
ggplot(aes(x=log(strain_abundance),y=weight, color = strain_num), data = dataset_by_strains[dataset_by_strains$treatment=="#1",])+
  geom_point()+
  geom_smooth(method='lm') +
  theme_bw() 
#! dev.off()

# now CommenCom
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_strain_type_correlation.pdf")
ggplot(aes(x=log(strain_abundance),y=weight, color = strain_num), data = dataset_by_strains[dataset_by_strains$treatment=="#2",])+
  geom_point()+
  geom_smooth(method='lm') +
  theme_bw() 
#! dev.off()

#Now validating the genotype by strain type plot,  by stan_glm function, and plotting the coeeficients results
patho_vs_commensal$genotype<- factor(x = patho_vs_commensal$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2")) #reordering the levels

weight_by_strain_type_by_gen_mod_patho <- stan_glm(formula = weight ~ log(load)*genotype + exp, data = patho_vs_commensal[patho_vs_commensal$strain_type=="OTU5",], seed=1234)
coeffcients_weight_by_strain_type_by_gen_patho <- c("log(load):genotypeEy15-2", "log(load):genotypeSchl-7", "log(load):genotypeLu3-30", "log(load):genotypeHE-1", "log(load):genotypeTu-Wal-2")
coeffcients_weight_by_strain_type_by_gen_sum_patho <-weight_by_strain_type_by_gen_mod_patho$stan_summary[coeffcients_weight_by_strain_type_by_gen_patho,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_weight_by_strain_type_by_gen_sum_patho <- as.data.frame(coeffcients_weight_by_strain_type_by_gen_sum_patho)
coeffcients_weight_by_strain_type_by_gen_sum_patho$condition <- rownames(coeffcients_weight_by_strain_type_by_gen_sum_patho)

weight_by_strain_type_by_gen_mod_commensal <- stan_glm(formula = weight ~ log(load)*genotype + exp, data = patho_vs_commensal[patho_vs_commensal$strain_type=="non-OTU5",], seed=1234)
coeffcients_weight_by_strain_type_by_gen_commensal <- c("log(load):genotypeEy15-2", "log(load):genotypeSchl-7", "log(load):genotypeLu3-30", "log(load):genotypeHE-1", "log(load):genotypeTu-Wal-2")
coeffcients_weight_by_strain_type_by_gen_sum_commensal <-weight_by_strain_type_by_gen_mod_commensal$stan_summary[coeffcients_weight_by_strain_type_by_gen_commensal,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_weight_by_strain_type_by_gen_sum_commensal <- as.data.frame(coeffcients_weight_by_strain_type_by_gen_sum_commensal)
coeffcients_weight_by_strain_type_by_gen_sum_commensal$condition <- rownames(coeffcients_weight_by_strain_type_by_gen_sum_commensal)

coeffcients_weight_by_strain_type_by_gen_sum_commensal$strain_type <- "non-OTU5"
coeffcients_weight_by_strain_type_by_gen_sum_patho$strain_type <- "OTU5"

coeffcients_weight_by_strain_type_by_gen_sum <- rbind(coeffcients_weight_by_strain_type_by_gen_sum_commensal,coeffcients_weight_by_strain_type_by_gen_sum_patho)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_strain_by_genotype_type_model_stan_validation.pdf", width = 7,height = 7)
ggplot(data = coeffcients_weight_by_strain_type_by_gen_sum, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color = strain_type))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
#! dev.off()


#############plotting genotype-wise, load correlting with weight, strain_by_strain (pathogen or commensal). This is to show genotype tolerance / resistance.
dataset_by_strains$strain_type[dataset_by_strains$strain_num%in%c(1:7)] <- "OTU5"
dataset_by_strains$strain_type[dataset_by_strains$strain_num%in%c(11:17)] <- "non-OTU5"

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_by_strain_number_correlation.pdf")
ggplot(aes(x=log(strain_abundance),y=weight, color = strain_num), data = dataset_by_strains)+
  geom_point(aes(color = strain_num))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(strain_type ~ .)
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_strain_number_correlation.pdf")
ggplot(aes(x=log(strain_abundance),y=weight, color = strain_num), data = dataset_by_strains)+
  geom_point(aes(color = strain_num))+
  geom_smooth(method='lm') +
  theme_bw() + 
  facet_grid(genotype ~ strain_type) # OPTIONAL
#! dev.off()

### Now validating the genotype by strain type plot,  by stan_glm function, and plotting the coeeficients results - all strains - regardless of OTU type

weight_by_strain_num_mod <- stan_glm(formula = weight ~ log(strain_abundance)*strain_num + genotype, 
                                          data = dataset_by_strains[dataset_by_strains$treatment%in%c("#1","#2") & dataset_by_strains$exp=="Oct",], seed=12345)
coefficients_stan_summary_weight_by_strain <- c("log(strain_abundance):strain_num2",	"log(strain_abundance):strain_num3",
                                                     "log(strain_abundance):strain_num4",	"log(strain_abundance):strain_num5",
                                                     "log(strain_abundance):strain_num6",	"log(strain_abundance):strain_num7",
                                                "log(strain_abundance):strain_num11","log(strain_abundance):strain_num12","log(strain_abundance):strain_num13",	
                                                        "log(strain_abundance):strain_num14", "log(strain_abundance):strain_num15",
                                                "log(strain_abundance):strain_num16", "log(strain_abundance):strain_num17")
mod_weight_by_strain <- weight_by_strain_num_mod$stan_summary[coefficients_stan_summary_weight_by_strain,c("2.5%","50%","97.5%")]
mod_weight_by_strain <- as.data.frame(mod_weight_by_strain)
mod_weight_by_strain$condition <- rownames(mod_weight_by_strain)

ggplot(data = mod_weight_by_strain ,aes(x=condition,y=`50%` , ymin = `2.5%`, ymax = `97.5%`, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

### Now validating the genotype by strain type plot,  by stan_glm function, and plotting the coeeficients results - within each OTU
# First to OTU5 strains
fit <- lmer(formula = weight ~ log(strain_abundance)*strain_num + exp + (1|genotype) + (1|tray), 
     data = dataset_by_strains[dataset_by_strains$strain_type=="OTU5",])
summary(fit)

weight_by_strain_num_mod_OTU5 <- stan_glm(formula = weight ~ log(strain_abundance)*strain_num + exp + genotype, 
                                          data = dataset_by_strains[dataset_by_strains$strain_type=="OTU5",], seed=12345)
coefficients_stan_summary_weight_by_strain_OTU5 <- c("log(strain_abundance):strain_num2",	"log(strain_abundance):strain_num3",
                                                "log(strain_abundance):strain_num4",	"log(strain_abundance):strain_num5",
                                                "log(strain_abundance):strain_num6",	"log(strain_abundance):strain_num7")
mod_weight_by_strain_OTU5 <- weight_by_strain_num_mod_OTU5$stan_summary[coefficients_stan_summary_weight_by_strain_OTU5,c("2.5%","50%","97.5%")]
mod_weight_by_strain_OTU5 <- as.data.frame(mod_weight_by_strain_OTU5)
mod_weight_by_strain_OTU5$condition <- rownames(mod_weight_by_strain_OTU5)
mod_weight_by_strain_OTU5$strain_OTU <- "OTU5"

# Now to non-OTU5 strains
weight_by_strain_num_mod_nonOTU5 <- stan_glm(formula = weight ~ log(strain_abundance)*strain_num + exp + genotype, 
                                          data = dataset_by_strains[dataset_by_strains$strain_type=="non-OTU5",], seed=12345)
coefficients_stan_summary_weight_by_strain_nonOTU5 <- c("log(strain_abundance):strain_num12","log(strain_abundance):strain_num13",	
                                                        "log(strain_abundance):strain_num14", "log(strain_abundance):strain_num15",
                                                        "log(strain_abundance):strain_num16", "log(strain_abundance):strain_num17")
mod_weight_by_strain_nonOTU5 <- weight_by_strain_num_mod_nonOTU5$stan_summary[coefficients_stan_summary_weight_by_strain_nonOTU5,c("2.5%","50%","97.5%")]
mod_weight_by_strain_nonOTU5 <- as.data.frame(mod_weight_by_strain_nonOTU5)
mod_weight_by_strain_nonOTU5$condition <- rownames(mod_weight_by_strain_nonOTU5)
mod_weight_by_strain_nonOTU5$strain_OTU <- "non-OTU5"

# Now i will present two plots with the results of both strain groups (OTU5 and non-OTU5).

OTU5_plot <- ggplot(data = mod_weight_by_strain_OTU5 ,aes(x=condition,y=`50%` , ymin = `2.5%`, ymax = `97.5%`, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

nonOTU5_plot <- ggplot(data = mod_weight_by_strain_nonOTU5 ,aes(x=condition,y=`50%` , ymin = `2.5%`, ymax = `97.5%`, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)



OTU5_and_nonOTU5_plot <- ggarrange(OTU5_plot, nonOTU5_plot,labels = c("A", "B"),
                    ncol = 1, nrow = 2)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_by_strain_number_correlation_stan_validation.pdf",width = 14,height = 7)
OTU5_and_nonOTU5_plot
#! dev.off()

## now by genotype
dataset_by_strains$genotype <- factor(x =dataset_by_strains$genotype, levels = c("Kus3-1", "Ey15-2", "Schl-7", "Lu3-30", "HE-1", "Tu-Wal-2")) #renaming host genotypes to original 1001 names

stan_results_weight_by_strain_by_gen <- data.frame("strain_num" = NA, "condition" = NA, "X2.5X"=NA,"X50X"=NA,"X97.5X"=NA)
coeffcients_stan_summary_stain_num_by_gen <- c("genotypeHE-1", "genotypeEy15-2",	"genotypeLu3-30",	"genotypeSchl-7",	"genotypeTu-Wal-2",
                                                "log(strain_abundance):genotypeHE-1",	"log(strain_abundance):genotypeEy15-2",
                                                "log(strain_abundance):genotypeLu3-30",	"log(strain_abundance):genotypeSchl-7",	"log(strain_abundance):genotypeTu-Wal-2")

for (strain_num in unique(dataset_by_strains$strain_num)){
  weight_by_strain_num_by_gen_mod <- stan_glm(formula = weight ~ log(strain_abundance)*genotype + exp , data = dataset_by_strains[dataset_by_strains$strain_num==strain_num,], seed=12345)
  mod_results_per_strain <- weight_by_strain_num_by_gen_mod$stan_summary[coeffcients_stan_summary_stain_num_by_gen,c("2.5%","50%","97.5%")]
  mod_results_per_strain <- as.data.frame(mod_results_per_strain)
  colnames(mod_results_per_strain) <- c("X2.5X","X50X","X97.5X")
  mod_results_per_strain$strain_num <- strain_num
  mod_results_per_strain$condition <- rownames(mod_results_per_strain)
  stan_results_weight_by_strain_by_gen <- rbind(stan_results_weight_by_strain_by_gen,mod_results_per_strain)
}

stan_results_weight_by_strain_by_gen <- stan_results_weight_by_strain_by_gen[!is.na(stan_results_weight_by_strain_by_gen$condition),]
stan_results_weight_by_strain_by_gen$condition <- as.factor(stan_results_weight_by_strain_by_gen$condition)
stan_results_weight_by_strain_by_gen$condition_type[grepl(x = stan_results_weight_by_strain_by_gen$condition, pattern = "log")] <- "log"
stan_results_weight_by_strain_by_gen$condition_type[!grepl(x = stan_results_weight_by_strain_by_gen$condition, pattern = "log")] <- "strain_only"

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/weight_by_genotype_by_strain_number_correlation_stan_validation.pdf",width = 7,height = 14)
ggplot(data = stan_results_weight_by_strain_by_gen[stan_results_weight_by_strain_by_gen$condition_type=="log",]
       ,aes(x=condition,y=X50X , ymin = X2.5X, ymax = X97.5X, color = strain_num))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(strain_num ~ .)
#! dev.off()


########################### load by genotype. a.k.a which genotype has higher / lower load, per total load, per group (OTU5/non-OTU5) and per strain
#First modeling, Baesiyan approach


load_by_treat_mod <- stan_glm(formula = log(bacterial_load) ~ treatment + genotype, data = full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug",])

#plot(load_by_treat_mod,plotfun = "mcmc_hist")
#load_by_treat_mod <- lmBF(formula =log_bacterial_load ~ treatment + exp*genotype, data = full_dataset_local_no_out[full_dataset_local_no_out],posterior = T, iter=10000)
#weight_by_load_mod <- stan_glm(formula = weight ~ log(bacterial_load)*treatment*genotype + exp, data = full_dataset_local_no_out)
#load_by_treat_mod<- as.data.frame(summary(load_by_treat_mod)[[2]])

coeffcients_interactions <- c("(Intercept)","treatment#2", "treatment#3")
coeffcients_interactions_summary <- load_by_treat_mod$stan_summary[coeffcients_interactions,c("2.5%","50%","97.5%")] # This is to examplify that the interactions of bacterial load and treamtent are insignificant, Baesiyan style.
coeffcients_interactions_summary <- as.data.frame(coeffcients_interactions_summary)
coeffcients_interactions_summary$condition <- rownames(coeffcients_interactions_summary)

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_treatment_model_stan_validation.pdf")
ggplot(data = coeffcients_interactions_summary, aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() 
#dev.off()


#plotting load of treatments by density plot - ONLY AUG experiment - representative experiment. Combining the results from stan model
dataset_patho <- full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment=="#1",]
dataset_commen <- full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment=="#2",]
dataset_mixed <- full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug" & full_dataset_local_no_out$treatment=="#3",]

#plotting load of treatments by density plot - both experiment.
#dataset_patho <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",]
#dataset_commen <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#2",]
#dataset_mixed <- full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",]




plot_dens <- function(dataset,coeffieint_df, treatment_name_coef, treat_color, x_ax_lables=T){
  xmin <- coeffieint_df$`2.5%`[coeffieint_df$condition==treatment_name_coef]
  xmax <- coeffieint_df$`97.5%`[coeffieint_df$condition==treatment_name_coef]
  median_coef <- coeffieint_df$`50%`[coeffieint_df$condition==treatment_name_coef]
  
  
  plot <- ggplot(aes(x=log(bacterial_load)), data =  dataset)+
    geom_density(colour = treat_color) +
    theme_bw()+ 
    facet_grid(treatment ~ .)
  
  dens <- ggplot_build(plot)$data[[1]]
  dens <- subset(dens, x < xmax)
  dens <- subset(dens, x > xmin)
  
  if (x_ax_lables==F){
    final_plot <- plot + 
      geom_area(data = dens, aes(x=x, y=y), fill= treat_color, alpha=0.5) +
      geom_vline(aes(xintercept= median_coef),
                 color="black", linetype="dashed", size=0.7)+
      xlim( c(-7.5,5))+
      ylim( c(0,0.27))+ 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      ylab("")
  } else{
    final_plot <- plot + 
      geom_area(data = dens, aes(x=x, y=y), fill= treat_color, alpha=0.5) +
      geom_vline(aes(xintercept= median_coef),
                 color="black", linetype="dashed", size=0.7)+
      xlim( c(-7.5,5))+
      ylim( c(0,0.27))+
      ylab("")
      
  }
  return(final_plot)
}

plot_dens_patho <- plot_dens(dataset = dataset_patho,coeffieint_df = coeffcients_interactions_summary,treatment_name_coef = "(Intercept)",treat_color = "#d95f02", x_ax_lables = F)
plot_dens_commen <- plot_dens(dataset = dataset_commen,coeffieint_df = coeffcients_interactions_summary,treatment_name_coef = "treatment#2",treat_color = "#1b9e77",x_ax_lables = F)
plot_dens_mixed <- plot_dens(dataset = dataset_mixed,coeffieint_df = coeffcients_interactions_summary,treatment_name_coef = "treatment#3",treat_color = "#7570b3")


density_treatments_plot <- ggarrange(plot_dens_patho, plot_dens_commen, plot_dens_mixed, ncol = 1)

#! pdf( "~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_treatment_August_density_with_Rstan_model.pdf")
print(density_treatments_plot)
#! dev.off()

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_treatment_August_density_plot.pdf")
ggplot(aes(x=log(bacterial_load), color = treatment), data =  full_dataset_local_no_out[full_dataset_local_no_out$exp=="Aug",])+
  geom_density() +
  theme_bw() +
  geom_ribbon(data=subset(full_dataset_local_no_out,log(bacterial_load)>-1 & log(bacterial_load)<1),aes(ymax=0.2),ymin=0,
              fill="red",colour=NA,alpha=0.5) +
  facet_grid(treatment ~ .) 
#! dev.off()


#facet by genotype
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_treatment_by_gen_August_density_plot.pdf")
ggplot(aes(x=log(load), color = treatment), data =  patho_vs_commensal[!patho_vs_commensal$load_type %in% c("pathogens_mixed","commensals_mixed") 
                                                                       & patho_vs_commensal$exp=="Aug",])+
  geom_density() +
  theme_bw() +
  facet_grid(genotype ~ .) # OPTIONAL
#! dev.off()


#plotting pathogens vs commensals by density plot (by strain type - OTU5 vs non-OTU5) - ONLY AUG experiment - representative experiment. Possible to facet by genotype.

#first - creating two new groups - exclusive communities (PathoCom & CommenCom), vs. the Mixed com
patho_vs_commensal$treatment_type[patho_vs_commensal$treatment %in% c("#1","#2")] <- "exclusive" 
patho_vs_commensal$treatment_type[patho_vs_commensal$treatment %in% c("#3")] <- "mixed" 

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_strain_type_treatment_August_density_plot.pdf")
ggplot(aes(x=log(load), color = strain_type), data =  patho_vs_commensal[!patho_vs_commensal$load_type %in% c("total_load_mixed") 
                                                                         & patho_vs_commensal$exp=="Aug",])+
  geom_density() +
  theme_bw() +
  facet_grid(treatment_type ~ .)
#! dev.off()

#now faceting also by genotype
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_strain_type_treatment_by_gen_August_density_plot.pdf")
ggplot(aes(x=log(load), color = strain_type), data =  patho_vs_commensal[!patho_vs_commensal$load_type %in% c("total_load_mixed") 
                                                                       & patho_vs_commensal$exp=="Aug",])+
  geom_density() +
  theme_bw() +
  facet_grid(genotype ~ treatment_type)
#! dev.off()





#plotting load of each strain by density plot - ONLY AUG experiment MixedCom - representative experiment, first pathogens. 
load_per_strain_by_gen_Aug_patho <- ggplot(aes(x=log(strain_abundance), color = genotype), 
       data =  dataset_by_strains[dataset_by_strains$exp=="Aug" & dataset_by_strains$treatment=="#3" & dataset_by_strains$strain_num%in%c(1:7),])+
  geom_density() +
  theme_bw() +
  facet_grid(strain_num ~ .) +
  coord_cartesian(xlim=c(-10, 2))

#plotting load of each strain by density plot - ONLY AUG experiment MixedCom - representative experiment, now commensals.
load_per_strain_by_gen_Aug_commensal <- ggplot(aes(x=log(strain_abundance), color = genotype), 
       data =  dataset_by_strains[dataset_by_strains$exp=="Aug" & dataset_by_strains$treatment=="#3" & dataset_by_strains$strain_num%in%c(11:17),])+
  geom_density() +
  theme_bw() +
  facet_grid(strain_num ~ .) +
  coord_cartesian(xlim=c(-10, 2))

load_per_strain_by_gen_Aug <- ggarrange(load_per_strain_by_gen_Aug_patho, load_per_strain_by_gen_Aug_commensal,labels = c("A", "B"),
                                   ncol = 2, nrow = 1)

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_strain_MixedCom_by_gen_August_density_plot.pdf", width = 10, height = 14)
load_per_strain_by_gen_Aug
#! dev.off()


#### only for strain 6, in both MixedCom (zoom in from the previous part)
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_strain6_MixedCom_by_gen_August_density_plot.pdf", width = 8, height = 14)
ggplot(aes(x=log(strain_abundance), color = genotype), 
       data =  dataset_by_strains[dataset_by_strains$exp=="Aug" & dataset_by_strains$treatment=="#3" & dataset_by_strains$strain_num%in%c(6),])+
  geom_density() +
  theme_bw() +
  coord_cartesian(xlim=c(-8, 2)) +
  facet_grid(genotype ~ .)
#! dev.off()

#### only for strain 6, in PathoCom (zoom in from the previous part)
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_strain6_MixedCom_by_gen_August_density_plot.pdf", width = 8, height = 14)
ggplot(aes(x=log(strain_abundance), color = genotype), 
       data =  dataset_by_strains[dataset_by_strains$exp=="Aug" & dataset_by_strains$treatment=="#1" & dataset_by_strains$strain_num%in%c(6),])+
  geom_density() +
  theme_bw() +
  coord_cartesian(xlim=c(-8, 4)) +
  facet_grid(genotype ~ .)
#! dev.off()



## now validating strain by strain load, Mixed Com, with BayesFactor library (basically can do an ANOVA test using the Baesiyan technique)
dataset_by_strains$strain_abundance_log <- log(dataset_by_strains$strain_abundance)

strain_load_Eyach_mod <- lmBF(formula = strain_abundance_log~ strain_num*genotype*exp, data = dataset_by_strains[dataset_by_strains$treatment%in% c("#1"),])
chains = posterior(strain_load_Eyach_mod, iterations = 10000)
chains_sum <- as.data.frame(summary(chains)[[2]])

coefficients_strain_load_by_gen_mod <- paste("strain_num-",c(1:7,11:17),sep = "")

mod_strain_load_by_gen <- chains_sum[coefficients_strain_load_by_gen_mod,c("2.5%","50%","97.5%")]
colnames(mod_strain_load_by_gen) <- c("X2.5X","X50X","X97.5X")
mod_strain_load_by_gen <- as.data.frame(mod_strain_load_by_gen)
mod_strain_load_by_gen$condition <- rownames(mod_strain_load_by_gen)
mod_strain_load_by_gen$condition <- factor(mod_strain_load_by_gen$condition, levels =  paste("strain_num-",c(1:7,11:17),sep = ""))
levels(mod_strain_load_by_gen$condition) <- c(1:7,11:17)
mod_strain_load_by_gen$OTU[mod_strain_load_by_gen$condition%in%c(1:7)] <- "OTU5"
mod_strain_load_by_gen$OTU[mod_strain_load_by_gen$condition%in%c(11:17)] <- "non-OTU5"


##now present the Stan model results
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_all_strains_MixedCom_August_density_plot_lmBF_validation.pdf",useDingbats = F)
ggplot(data = mod_strain_load_by_gen ,aes(x=condition,y=X50X , ymin = X2.5X, ymax = X97.5X, color = OTU))+
  geom_pointrange(position=position_dodge(width=1)) +
  scale_color_manual(values=c("#1b9e77","#d95f02")) +
  theme_bw()
#dev.off()




## now validating strain by strain 6 load in Mixed Com, between genotypes with BayesFactor library (basically can do an ANOVA test using the Baesiyan technique)
dataset_by_strains$strain_abundance_log <- log(dataset_by_strains$strain_abundance)
strain_6_load_mod <- lmBF(formula = strain_abundance_log~ genotype*exp, data = dataset_by_strains[dataset_by_strains$treatment%in% c("#3") & dataset_by_strains$strain_num==6,], posterior = T, iter=10000)
strain_6_load_mod_df <- as.data.frame(summary(strain_6_load_mod)[[2]])
coefficients_genotypes_strain_6 <- paste("genotype-",unique(dataset_by_strains$genotype),sep = "")
strain_6_load_mod_df <- strain_6_load_mod_df[coefficients_genotypes_strain_6 ,c("2.5%","50%","97.5%")]
strain_6_load_mod_df$genotype <- gsub("genotype-", "", rownames(strain_6_load_mod_df))

##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_strain6_MixedCom_August_lmBF_validation_dim.pdf",useDingbats = F, width = 8, height = 10)
ggplot(data = strain_6_load_mod_df ,aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1), color = "#a50026") +
  theme_bw() 
#! dev.off()



### now genotype by genotype, per genotype
final_strain_by_strain_by_gen_df <- data.frame()
coefficients_strain_load_by_gen_mod <- paste("strain_num-",c(1:7,11:17),sep = "")

for (gen in unique(dataset_by_strains$genotype)){
  
  strain_load_Eyach_mod <- lmBF(formula = strain_abundance_log~ strain_num, data = dataset_by_strains[dataset_by_strains$treatment%in% c("#3") & dataset_by_strains$genotype==gen & dataset_by_strains$exp=="Aug",])
  chains = posterior(strain_load_Eyach_mod, iterations = 10000)
  chains_sum <- as.data.frame(summary(chains)[[2]])
  temp_df <- as.data.frame(chains_sum[coefficients_strain_load_by_gen_mod,c("2.5%","50%","97.5%")])
  temp_df$genotype <- gen
  temp_df$strain_num <- rownames(temp_df)
  final_strain_by_strain_by_gen_df <- rbind(final_strain_by_strain_by_gen_df,temp_df)
}

final_strain_by_strain_by_gen_df$strain_num <- gsub("strain_num-", "", final_strain_by_strain_by_gen_df$strain_num)
final_strain_by_strain_by_gen_df$strain_num <- factor(final_strain_by_strain_by_gen_df$strain_num,levels = c(1:7,11:17))


final_strain_by_strain_by_gen_df$OTU[final_strain_by_strain_by_gen_df$strain_num%in%c(1:7)] <- "OTU5"
final_strain_by_strain_by_gen_df$OTU[final_strain_by_strain_by_gen_df$strain_num%in%c(11:17)] <- "non-OTU5"

#OPTIONAL: order by strain load
#order_by_strain_load <- final_strain_by_strain_by_gen_df$strain_num[order(final_strain_by_strain_by_gen_df$`50%`[final_strain_by_strain_by_gen_df$genotype=="Ey15-2"])]
library(plyr)
final_strain_by_strain_by_gen_df$strain_num <- factor(final_strain_by_strain_by_gen_df$strain_num,levels = order_by_strain_load)


df2.5_sum<- ddply(final_strain_by_strain_by_gen_df[final_strain_by_strain_by_gen_df$strain_num==6,], "genotype", summarize, X2.5_sum = sum(`2.5%`))
df97.5_sum<- ddply(final_strain_by_strain_by_gen_df[final_strain_by_strain_by_gen_df$strain_num==6,], "genotype", summarize, X97.5_sum = sum(`97.5%`))
df_sum_stan <- df2.5_sum
df_sum_stan$X97.5_sum <- df97.5_sum$X97.5_sum

##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/all_strains_MixedCom_August_by_gen_plot_lmBF_validation.pdf",useDingbats = F)
ggplot(data = final_strain_by_strain_by_gen_df ,aes(x=strain_num,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color = OTU))+
  geom_pointrange(position=position_dodge(width=1)) +
  #geom_path(aes(group=genotype)) +
  scale_color_manual(values=c("#1b9e77","#d95f02")) +
  geom_hline(data = df2.5_sum, aes(yintercept=X2.5_sum)) +
  geom_hline(data = df97.5_sum, aes(yintercept=X97.5_sum)) +
  theme_bw() +
  facet_grid(genotype ~ .)
#! dev.off()
  
dataset_by_strains$strain_abundance_log <- log(dataset_by_strains$strain_abundance)

### now strain by strain, per genotype
final_strain_by_strain_by_gen_df <- data.frame()
coefficients_strain_load_by_gen_mod <- paste("genotype-",unique(dataset_by_strains$genotype),sep = "")

for (strain in c(1:7)){
  
  strain_load_Eyach_mod <- lmBF(formula = strain_abundance_log~ genotype, data = dataset_by_strains[dataset_by_strains$treatment%in% c("#1") & dataset_by_strains$strain_num==strain & dataset_by_strains$exp=="Aug",])
  chains = posterior(strain_load_Eyach_mod, iterations = 10000)
  chains_sum <- as.data.frame(summary(chains)[[2]])
  temp_df <- as.data.frame(chains_sum[coefficients_strain_load_by_gen_mod,c("2.5%","50%","97.5%")])
  temp_df$strain_num <- strain
  temp_df$genotype <- rownames(temp_df)
  final_strain_by_strain_by_gen_df <- rbind(final_strain_by_strain_by_gen_df,temp_df)
}


final_strain_by_strain_by_gen_df <- data.frame()
coefficients_strain_load_by_gen_mod <- paste("genotype-",unique(dataset_by_strains$genotype),sep = "")

for (treat in unique(full_dataset_local_no_out$treatment)){
  total_load_by_gen <- lmBF(formula = log_bacterial_load~ genotype+exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in% treat,],posterior = T,iter=10000)
  total_load_df <- as.data.frame(summary(total_load_by_gen)[[2]])
  total_load_df <- total_load_df[coefficients_strain_load_by_gen_mod,c("2.5%","50%","97.5%")]
  total_load_df$treatment <- treat
  total_load_df$genotype <- rownames(total_load_df)
  final_strain_by_strain_by_gen_df <- rbind(final_strain_by_strain_by_gen_df,total_load_df)
}


final_strain_by_strain_by_gen_df$genotype <- gsub("genotype-", "", final_strain_by_strain_by_gen_df$genotype)
final_strain_by_strain_by_gen_df$genotype <- factor(final_strain_by_strain_by_gen_df$genotype)


##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/total_bacterial_load__per_genotype_all_treatments_posterior_lmBF_validation_dim.pdf", useDingbats = F, width = 7, height = 10)
ggplot(data = final_strain_by_strain_by_gen_df ,aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color=treatment))+
  geom_pointrange(position=position_dodge(width=1)) +
  scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +
  theme_bw() +
  facet_grid(treatment ~ .)
#! dev.off()




### now all strains vs all, abundance, in pathocom
all_strains_abundance_mod <- lmBF(formula = strain_abundance_log ~ strain_num*exp, data = dataset_by_strains[dataset_by_strains$treatment=="#1",], posterior = T, iter=10000)
all_strains_abundance_df <- as.data.frame(summary(all_strains_abundance_mod)[[2]])
all_strains_abundance_df <- all_strains_abundance_df[paste("strain_num-",c(1:7),sep = ""),c("2.5%","50%","97.5%")]
all_strains_abundance_df$strain_num <- as.factor(paste("P",c(1:7),sep = ""))

##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_all_strains_no_genotype_PathoCom_posterior_lmBF_validation.pdf", useDingbats = F, width = 7, height = 10)
ggplot(data = all_strains_abundance_df,aes(x=strain_num,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1), color="#a50026") +
  theme_bw() 
#! dev.off()



### now all strains vs all, abundance, in MixedCom
all_strains_abundance_mod <- lmBF(formula = strain_abundance_log ~ strain_num, data = dataset_by_strains[dataset_by_strains$treatment=="#3" & dataset_by_strains$exp=="Aug",], posterior = T, iter=10000)
all_strains_abundance_df <- as.data.frame(summary(all_strains_abundance_mod)[[2]])
all_strains_abundance_df <- all_strains_abundance_df[paste("strain_num-",c(1:7,11:17),sep = ""),c("2.5%","50%","97.5%")]
all_strains_abundance_df$strain_num <- as.factor(c(paste("P",c(1:7),sep = ""),paste("C",c(1:7),sep = "")))
all_strains_abundance_df$OTU <- c(rep("OTU5",7),rep("non-OTU5",7))
all_strains_abundance_df$strain_num <- factor(all_strains_abundance_df$strain_num, levels = c(paste("P",c(1:7),sep = ""),paste("C",c(1:7),sep = "")))

##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_all_strains_no_genotype_MixedCom_posterior_lmBF_validation_only_Aug.pdf", useDingbats = F)
ggplot(data = all_strains_abundance_df,aes(x=strain_num,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color=OTU))+
  geom_pointrange(position=position_dodge(width=1)) +
  scale_color_manual(values = c("#4393c3","#a50026")) +
  theme_bw() 
#! dev.off()

  
  



### now all strains vs all, abundance, in pathocom
all_strains_abundance_mod <- lmBF(formula = strain_abundance_log ~ strain_num*exp, data = dataset_by_strains[dataset_by_strains$treatment=="#1",], posterior = T, iter=10000)
all_strains_abundance_df <- as.data.frame(summary(all_strains_abundance_mod)[[2]])
all_strains_abundance_df <- all_strains_abundance_df[paste("strain_num-",c(1:7),sep = ""),c("2.5%","50%","97.5%")]
all_strains_abundance_df$strain_num <- as.factor(paste("P",c(1:7),sep = ""))

##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_all_strains_no_genotype_PathoCom_posterior_lmBF_validation.pdf", useDingbats = F, width = 7, height = 10)
ggplot(data = all_strains_abundance_df,aes(x=strain_num,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1), color="#a50026") +
  theme_bw() 
#! dev.off()



### now all strains vs all, abundance, in pathocom
all_strains_abundance_mod <- lmBF(formula = X6 ~ genotype, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3" & full_dataset_local_no_out$exp=="Aug",], posterior = T, iter=10000)
all_strains_abundance_df <- as.data.frame(summary(all_strains_abundance_mod)[[2]])
all_strains_abundance_df <- all_strains_abundance_df[paste("genotype-",unique(full_dataset_local_no_out$genotype),sep = ""),c("2.5%","50%","97.5%")]
all_strains_abundance_df$genotype <- gsub("genotype-", "", rownames(all_strains_abundance_df))
all_strains_abundance_df$genotype <- factor(all_strains_abundance_df$genotype)

##now present the Stan model results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_all_strains_no_genotype_PathoCom_posterior_lmBF_validation.pdf", useDingbats = F, width = 7, height = 10)
ggplot(data = all_strains_abundance_df,aes(x=genotype,y=`50%`, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(position=position_dodge(width=1), color="#a50026") +
  theme_bw() 
#! dev.off()






plot_dens_strain_load <- function(dataset,coeffieint_df, treatment_name_coef, treat_color, x_ax_lables=T, strain=6){
  xmin <- stan_mod_strain_by_gen_MixedCom$X2.5X[stan_mod_strain_by_gen_MixedCom$strain==strain & stan_mod_strain_by_gen_MixedCom$condition==treatment_name_coef]
  xmax <- stan_mod_strain_by_gen_MixedCom$X97.5X[stan_mod_strain_by_gen_MixedCom$strain==strain & stan_mod_strain_by_gen_MixedCom$condition==treatment_name_coef]
  median_coef <- stan_mod_strain_by_gen_MixedCom$X50X[stan_mod_strain_by_gen_MixedCom$strain==strain & stan_mod_strain_by_gen_MixedCom$condition==treatment_name_coef]
  
  
  plot <- ggplot(aes(x=strain_abundance_log), data =  dataset[dataset$strain_num==strain,])+
    geom_density(colour = treat_color) +
    theme_bw()+ 
    xlim(c(-8,3)) +
    facet_grid(treatment ~ .)
  
  dens <- ggplot_build(plot)$data[[1]]
  dens <- subset(dens, x < xmax)
  dens <- subset(dens, x > xmin)
  
  if (x_ax_lables==F){
    final_plot <- plot + 
      geom_area(data = dens, aes(x=x, y=y), fill= treat_color, alpha=0.5) +
      geom_vline(aes(xintercept= median_coef),
                 color="black", linetype="dashed", size=0.7)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      ylab("")
  } else{
    final_plot <- plot + 
      geom_area(data = dens, aes(x=x, y=y), fill= treat_color, alpha=0.5) +
      geom_vline(aes(xintercept= median_coef),
                 color="black", linetype="dashed", size=0.7)+
      ylab("")
    
  }
  return(final_plot)
}


plot_dens_Ey <- plot_dens_strain_load(dataset = dataset_by_strains[dataset_by_strains$treatment%in% c("#3") & dataset_by_strains$genotype=="Ey15-2",],
                                      coeffieint_df = stan_mod_strain_by_gen_MixedCom,treatment_name_coef = "genotype-Ey15-2",treat_color = "#d95f02", x_ax_lables = T,strain = 1)



density_treatments_plot <- ggarrange(plot_dens_TuW, plot_dens_HE, plot_dens_Schl,plot_dens_Ey, plot_dens_Kus,plot_dens_Lu,ncol = 1)











## now validating total abundance by genotype in Mixed Com, with BayesFactor library (basically can do an ANOVA test using the Baesiyan technique)
total_load_by_gen_mod <- lmBF(formula = log_bacterial_load  ~ genotype*exp, data = full_dataset_local_no_out[full_dataset_local_no_out$treatment%in% c("#3"),])
chains = posterior(total_load_by_gen_mod, iterations = 10000)
## 1:13 are the only "interesting" parameters
chains_sum <- as.data.frame(summary(chains)[[2]])

coefficients_total_load_by_gen_mod <- c("genotype-Kus3-1","genotype-Ey15-2", "genotype-Schl-7",	"genotype-Lu3-30",	"genotype-HE-1",	"genotype-Tu-Wal-2")
mod_total_load_by_gen <- chains_sum[coefficients_total_load_by_gen_mod,c("2.5%","50%","97.5%")]
colnames(mod_total_load_by_gen) <- c("X2.5X","X50X","X97.5X")
mod_total_load_by_gen <- as.data.frame(mod_total_load_by_gen)
mod_total_load_by_gen$condition <- rownames(mod_total_load_by_gen)





# Now i will present two plots with the results of both strain groups (OTU5 and non-OTU5).
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_strain6_MixedCom_by_gen_August_density_plot_stan_validation.pdf")
ggplot(data = mod_total_load_by_gen ,aes(x=condition,y=X50X , ymin = X2.5X, ymax = X97.5X, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)
#! dev.off()







#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_all_strains_MixedCom_by_gen_August_density_plot_stan_validation.pdf", width = 7, height = 14)
ggplot(data = stan_mod_strain_by_gen_MixedCom ,aes(x=condition,y=X50X , ymin = X2.5X, ymax = X97.5X, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  facet_grid(strain ~ .)
#! dev.off()


##now validation of which strain is generally most dominant in MixedCom, both experiments.
load_by_strain_Mixed_com_mod <- stan_glm(formula = strain_abundance ~ strain_num + genotype + exp, data = dataset_by_strains[dataset_by_strains$treatment=="#3",])
coefficients_strain_load_by_gen_mod <- c("(Intercept)","strain_num2", "strain_num3", "strain_num4","strain_num5","strain_num6", "strain_num7",
                                         "strain_num11","strain_num12","strain_num13","strain_num14","strain_num15","strain_num16", "strain_num17")
mod_strain_load_by_gen <- load_by_strain_Mixed_com_mod$stan_summary[coefficients_strain_load_by_gen_mod,c("2.5%","50%","97.5%")]
mod_strain_load_by_gen <- as.data.frame(mod_strain_load_by_gen)
mod_strain_load_by_gen$condition <- rownames(mod_strain_load_by_gen)
mod_strain_load_by_gen[mod_strain_load_by_gen$condition=="(Intercept)",1:3] <- mod_strain_load_by_gen[mod_strain_load_by_gen$condition=="(Intercept)",1:3]-mod_strain_load_by_gen$`50%`[mod_strain_load_by_gen$condition=="(Intercept)"]


#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_strain_MixedCom_stan_validation.pdf")
ggplot(data = mod_strain_load_by_gen,aes(x=condition,y=`50%` , ymin = `2.5%`, ymax = `97.5%`, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)
#! dev.off()


##now validation of which strain is generally most dominant in PathoCom, both experiments.
dataset_by_strains$strain_num <- factor(dataset_by_strains$strain_num,levels=c(1:7,11:17))
load_by_strain_Patho_com_mod <- stan_glm(formula = strain_abundance ~ strain_num + genotype + exp, data = dataset_by_strains[dataset_by_strains$treatment=="#1",])
coefficients_strain_load_by_gen_mod <- c("(Intercept)","strain_num2", "strain_num3", "strain_num4","strain_num5","strain_num", "strain_num7")
mod_strain_load_in_pathocom <- load_by_strain_Patho_com_mod$stan_summary[coefficients_strain_load_by_gen_mod,c("2.5%","50%","97.5%")]
mod_strain_load_in_pathocom <- as.data.frame(mod_strain_load_in_pathocom)
mod_strain_load_in_pathocom$condition <- rownames(mod_strain_load_in_pathocom)
mod_strain_load_in_pathocom[mod_strain_load_in_pathocom$condition=="(Intercept)",1:3] <- mod_strain_load_in_pathocom[mod_strain_load_in_pathocom$condition=="(Intercept)",1:3]-mod_strain_load_in_pathocom$`50%`[mod_strain_load_in_pathocom$condition=="(Intercept)"]

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_treatment_and_genotype/load_by_strain_PathoCom_stan_validation.pdf")
ggplot(data = mod_strain_load_in_pathocom,aes(x=condition,y=`50%` , ymin = `2.5%`, ymax = `97.5%`, color = condition))+
  geom_pointrange(position=position_dodge(width=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)
#! dev.off()




########################################################## association between individual strains abundance to treatment*genotype ########################
########################################################## a.k.a the random forest part ##################################################################
###run random forest to predict weight, for all communities, using both exp and all possible predictors
full_dataset_local_no_out$exp <- as.factor(full_dataset_local_no_out$exp)

rf_weight_by_strain_all <- randomForest(y = weight, X1 + X2 + X3 + X4 + X5 + X6 + X7 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + genotype + bacterial_load + tray + pot + exp + treatment, full_dataset_local_no_out,seed = 12345)

                                        
#rf_weight_by_strain <- randomForest(weight ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + bacterial_load + com3_patho + com3_commensal + exp + tray + pot, full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",])
importance(rf_weight_by_strain_all)
#now buildin a dataframe with imporatance results, and order it by importance
rf_importance_all <- data.frame("predictor" = row.names(rf_weight_by_strain_all$importance), "IncNodePurity" = rf_weight_by_strain_all$importance)
order_importance_all <- rf_importance_all$predictor[order(rf_importance_all$IncNodePurity)]
rf_importance_all$predictor <- factor(rf_importance_all$predictor, levels = order_importance_all)

#ploting random forest results
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/random_forest_all_treatments_all_predictors.pdf")
ggplot(aes(x = IncNodePurity, y = predictor),data = rf_importance_all)+
  geom_point() +
  theme_bw()
#! dev.off()


###run random forest to predict weight, only in com#1, using both exp and all possible predictors
full_dataset_local_no_out$exp <- as.factor(full_dataset_local_no_out$exp)
rf_weight_by_strain_com1 <- randomForest(weight ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + bacterial_load + tray + pot + exp, full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#1",],seed = 12345)
#rf_weight_by_strain <- randomForest(weight ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + bacterial_load + com3_patho + com3_commensal + exp + tray + pot, full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",])
importance(rf_weight_by_strain_com1)
#now buildin a dataframe with imporatance results, and order it by importance
rf_importance_com1 <- data.frame("predictor" = row.names(rf_weight_by_strain_com1$importance), "IncNodePurity" = rf_weight_by_strain_com1$importance)
order_importance_com1 <- rf_importance_com1$predictor[order(rf_importance_com1$IncNodePurity)]
rf_importance_com1$predictor <- factor(rf_importance_com1$predictor, levels = order_importance_com1)

#ploting random forest results
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/random_forest_com1_all_predictors_no_genotype.pdf", useDingbats = F)
ggplot(aes(x = IncNodePurity, y = predictor),data = rf_importance_com1)+
  geom_point() +
  theme_bw()
#dev.off()


###run random forest to predict weight, only in com#3, using both exp and all possible predictors
full_dataset_local_no_out$exp <- as.factor(full_dataset_local_no_out$exp)
#rf_weight_by_strain <- randomForest(weight ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + bacterial_load + genotype + com3_patho + com3_commensal + exp + tray + pot, full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",], seed = 12345)
rf_weight_by_strain <- randomForest(weight ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + bacterial_load + com3_patho + com3_commensal + exp + tray + pot, full_dataset_local_no_out[full_dataset_local_no_out$treatment=="#3",], seed = 12345, importance=T)

#now buildin a dataframe with imporatance results, and order it by importance
rf_importance <- data.frame("predictor" = row.names(rf_weight_by_strain$importance), "IncNodePurity" = rf_weight_by_strain$importance)
order_importance <- rf_importance$predictor[order(rf_importance$IncNodePurity..IncMSE)]
rf_importance$predictor <- factor(rf_importance$predictor, levels = order_importance)

#ploting random forest results
#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/random_forest_com3_all_predictors_no_genotype_Inc_MSE.pdf", useDingbats = F)
ggplot(aes(x = IncNodePurity..IncMSE, y = predictor),data = rf_importance)+
  geom_point() +
  theme_bw()
#dev.off()


order_importance <- rf_importance$predictor[order(rf_importance$IncNodePurity.IncNodePurity)]
rf_importance$predictor <- factor(rf_importance$predictor, levels = order_importance)

#pdf("~/ownCloud/My papers/Syncoms_paper/Figures/load_by_weight_correlation/random_forest_com3_all_predictors.pdf", useDingbats = F)
ggplot(aes(x = IncNodePurity.IncNodePurity, y = predictor),data = rf_importance)+
  geom_point() +
  theme_bw()
#dev.off()

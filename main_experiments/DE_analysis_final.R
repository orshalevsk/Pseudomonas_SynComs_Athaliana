library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(qvalue)

#get reads count table
counts_table <- read.csv2('~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/RNA-seq_Exp_July_2019_rep2/results/counts_for_R.csv', sep = '\t')
drops <- c("Chr","Start","End","Strand","Length")
final_table_only_genes_and_counts <- counts_table[,!(names(counts_table) %in% drops)]

counts_table_second <- read.csv2('~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/RNA-seq_Exp_July_2019_rep2/results/counts_second_demultiplex_for_R.csv')
names_updated <- gsub(pattern = "Lane",replacement = "p", names(counts_table_second))
names_updated <- gsub(pattern = ".bam", replacement = ".Aligned.sortedByCoord.out.bam", names_updated)
names(counts_table_second) <- names_updated
drops <- c("Chr","Start","End","Strand","Length")
final_table_only_genes_and_counts_second <- counts_table_second[,!(names(counts_table_second) %in% drops)]

missing_name <- setdiff(names(final_table_only_genes_and_counts),names(final_table_only_genes_and_counts_second))

final_table_only_genes_and_counts_second$p2.F8.Aligned.sortedByCoord.out.bam <- 0

for (name in names(final_table_only_genes_and_counts[,-1])){
  final_table_only_genes_and_counts[,name] <- final_table_only_genes_and_counts[,name] + final_table_only_genes_and_counts_second[,name]
}



#filtering samples with less than X reads
final_table_only_genes<- final_table_only_genes_and_counts[,-1]
only_high_depth <- colSums(final_table_only_genes) > 3000000
only_high_depth_names <- colnames(final_table_only_genes[only_high_depth])

final_table_only_genes_and_counts <- final_table_only_genes_and_counts[,c("Geneid",only_high_depth_names)]


#now create the metadata final table
all_tray_design <- read.csv2('/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/RNA-seq_Exp_July_2019_rep2/plan/all_trays_design.csv')
metadata <- read.csv2('/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/RNA-seq_Exp_July_2019_rep2/results/metadata_sequencing_lane.csv')
final_metadata_table <- merge(x = all_tray_design,y = metadata,by ="tray_position")


extra_wells <- setdiff(colnames(final_table_only_genes_and_counts),(final_metadata_table$sample_file_name))
extra_wells <- extra_wells[extra_wells!="Geneid"]
final_table_only_genes_and_counts <- final_table_only_genes_and_counts[,!(names(final_table_only_genes_and_counts) %in% extra_wells)]


extra_wells <- setdiff((final_metadata_table$sample_file_name),colnames(final_table_only_genes_and_counts))
final_metadata_table <- final_metadata_table[!(final_metadata_table$sample_file_name%in%extra_wells),]


final_metadata_table$Treatment <- factor(final_metadata_table$Treatment, levels = c("Control", "#1","#2", "#2+1", "#2+6", "#3", "#3-s1", "#3-s6"))
final_metadata_table$Perfect.barcode.... <- as.numeric(gsub(pattern = "%", replacement = "", x = as.character(final_metadata_table$Perfect.barcode....)))
final_metadata_table$Trimmed.bases.... <- as.numeric(gsub(pattern = "%", replacement = "", x = as.character(final_metadata_table$Trimmed.bases....)))
final_metadata_table$Mean.quality.score <- as.numeric(gsub(pattern = "%", replacement = "", x = as.character(final_metadata_table$Mean.quality.score)))
final_metadata_table$X....Q30.bases <- as.numeric(gsub(pattern = "%", replacement = "", x = as.character(final_metadata_table$X....Q30.bases)))
final_metadata_table$X..of.lane <- as.numeric(gsub(pattern = "%", replacement = "", x = as.character(final_metadata_table$X..of.lane)))



####subsetting to specific treatments / genotypes
table(final_metadata_table$Treatment)


#final_metadata_table_original <- final_metadata_table
final_metadata_table <- final_metadata_table[final_metadata_table$Treatment%in%c("Control","#1","#2","#3"),]
final_metadata_table$Treatment <- factor(final_metadata_table$Treatment)

extra_col <- setdiff(colnames(final_table_only_genes_and_counts),final_metadata_table$sample_file_name)
extra_col  <- extra_col[extra_col!="Geneid"]

final_table_only_genes_and_counts <- final_table_only_genes_and_counts[,!(names(final_table_only_genes_and_counts) %in% extra_col )]

samples_depth <- colSums(final_table_only_genes_and_counts[-1])
final_metadata_table$sample_file_name
samples_depth_df <- data.frame("sample_file_name"=names(samples_depth), "depth"=as.numeric(samples_depth))

final_metadata_table <- merge(final_metadata_table, samples_depth_df, by="sample_file_name")


final_metadata_table <- final_metadata_table[final_metadata_table$Genotype=="Lu",]

ggplot(data = final_metadata_table, aes(x = Treatment, y=depth, color=Genotype))+
  geom_point()
#now i will subset to the best 8 samples per treatment - 4 from each genotype


final_names <- vector()
for (treat in unique(final_metadata_table$Treatment)){
  for (gen in unique(final_metadata_table$Genotype)){
    final_metadata_table_temp <- final_metadata_table[final_metadata_table$Treatment==treat & final_metadata_table$Genotype==gen,]
    final_names_temp <- final_metadata_table_temp$sample_file_name[order(final_metadata_table_temp$depth,decreasing = T)][1:4]
    final_names <- c(final_names, as.character(final_names_temp))
  }
}


final_names<-final_names[!is.na(final_names)]

final_metadata_table <- final_metadata_table[final_metadata_table$sample_file_name %in% final_names,]
final_table_only_genes_and_counts <- final_table_only_genes_and_counts[,c("Geneid",final_names)]


#remove genes with less than average 5 hits // no hits from any sample 

filter <- apply(final_table_only_genes_and_counts[,-1], 1, mean)
final_table_only_genes_and_counts <- final_table_only_genes_and_counts[filter>5,]

hist(colSums(final_table_only_genes_and_counts[,-1]), breaks = 100)

final_table_only_genes_and_counts[final_table_only_genes_and_counts==0] <- 1


#dds <- DESeqDataSetFromMatrix(countData=final_table_only_genes_and_counts, 
#                              colData=final_metadata_table, 
#                              design=~ Genotype + Treatment + Time_point, tidy = TRUE)


dds <- DESeqDataSetFromMatrix(countData=final_table_only_genes_and_counts, 
                              colData=final_metadata_table, 
                              design= ~ Treatment + Time_point, tidy = TRUE)


dds_res <- DESeq(object = dds)





#------------------------------------ loading data part is done -------------------------------------------#

# plot the expression of specific genes

PR1 <- plotCounts(dds_res, gene="AT2G14610", intgroup="Treatment", returnData=T) #PR1
PR1$gene <- "PR1"
PR1$gene_code <- "AT2G14610"
PR1$process <- "SA/SAR" 

ICS1 <- plotCounts(dds_res, gene="AT1G74710", intgroup="Treatment", returnData=T) #ICS1
ICS1$gene <- "ICS1"
ICS1$gene_code <- "AT1G74710"
ICS1$process <- "SA/SAR" 

EDS1 <- plotCounts(dds_res, gene="AT3G48090", intgroup="Treatment", returnData=T) #EDS1
EDS1$gene <- "EDS1"
EDS1$gene_code <- "AT3G48090"
EDS1$process <- "SA/SAR" 


PR5 <- plotCounts(dds_res, gene="AT1G75040", intgroup="Treatment", returnData=T) #PR5
PR5$gene <- "PR5"
PR5$gene_code <- "AT1G75040"
PR5$process <- "SA/SAR" 

PR2 <- plotCounts(dds_res, gene="AT3G57260", intgroup="Treatment", returnData=T) #PR2
PR2$gene <- "PR2"
PR2$gene_code <- "AT3G57260"
PR2$process <- "SA/SAR" 

WRKY53 <- plotCounts(dds_res, gene="AT4G23810", intgroup="Treatment", returnData=T) #WRKY53
WRKY53$gene <- "WRKY53"
WRKY53$gene_code <- "AT4G23810"
WRKY53$process <- "SA/SAR" 



DEG_of_interest <- rbind(ICS1, EDS1, WRKY53, PR1, PR2, PR5)


samples <- rownames(DEG_of_interest)
samples <- gsub(pattern = "bam.*","bam",samples)
DEG_of_interest$time_point <- as.character(final_metadata_table$Time_point[final_metadata_table$sample_file_name %in% samples])

DEG_of_interest_summary <- data.frame()
for (gene in unique(DEG_of_interest$gene)){
  for (treatment in unique(DEG_of_interest$Treatment)){
    mean_expression <- mean(DEG_of_interest$count[DEG_of_interest$gene==gene & DEG_of_interest$Treatment==treatment])
    sd_expression <- sd(DEG_of_interest$count[DEG_of_interest$gene==gene & DEG_of_interest$Treatment==treatment])
    sem_expression <- sd_expression/sqrt(length(DEG_of_interest$count[DEG_of_interest$gene==gene & DEG_of_interest$Treatment==treatment]))
    process <- unique(DEG_of_interest$process[DEG_of_interest$gene==gene])
    temp <- data.frame(gene, treatment, mean_expression, sem_expression, process)
    DEG_of_interest_summary <- rbind(DEG_of_interest_summary,temp)
  }
}


DEG_of_interest_summary$treatment <- factor(DEG_of_interest_summary$treatment, levels = c("Control","#1","#2","#3"))
levels(DEG_of_interest_summary$treatment) <- c("Control","PathoCom","CommenCom","MixedCom")


#! pdf(file = "~/ownCloud/My papers/Syncoms_paper/Final_edited_figures/figure5/DEG_of_interest_only_SA_transformed.pdf")
ggplot(data=DEG_of_interest_summary[DEG_of_interest_summary$process=="SA/SAR",], aes(x = mean_expression, y = gene, fill=treatment))+
  scale_x_continuous(trans = "log10") +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(xmin=mean_expression-sem_expression, xmax=mean_expression+sem_expression), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#969696","#d95f02", "#1b9e77", "#7570b3")) 
#! dev.off()

#! pdf(file = "~/ownCloud/My papers/Syncoms_paper/Final_edited_figures/figure5/DEG_of_interest_only_SA_transformed.pdf")
ggplot(data=DEG_of_interest_summary[DEG_of_interest_summary$process=="HOUSEKEEPING",], aes(x = mean_expression, y = gene, fill=treatment))+
  scale_x_continuous(trans = "log10") +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(xmin=mean_expression-sem_expression, xmax=mean_expression+sem_expression), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#969696","#d95f02", "#1b9e77", "#7570b3")) 
#! dev.off()

#! pdf(file = "~/ownCloud/My papers/Syncoms_paper/Final_edited_figures/figure5/DEG_of_interest_only_SA_transformed.pdf")
ggplot(data=DEG_of_interest_summary[DEG_of_interest_summary$process=="JA",], aes(x = mean_expression, y = gene, fill=treatment))+
  scale_x_continuous(trans = "log10") +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(xmin=mean_expression-sem_expression, xmax=mean_expression+sem_expression), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#969696","#d95f02", "#1b9e77", "#7570b3")) 
#! dev.off()


# now transformed plot (visually)
#! pdf(file = "~/ownCloud/My papers/Syncoms_paper/Final_edited_figures/figure5/DEG_of_interest_only_SA.pdf")
ggplot(data=DEG_of_interest_summary[DEG_of_interest_summary$process=="SA/SAR",], aes(x = gene, y = mean_expression, fill=treatment))+
  scale_y_continuous(trans = "log10") +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_expression-sem_expression, ymax=mean_expression+sem_expression), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#969696","#d95f02", "#1b9e77", "#7570b3"))
#! dev.off()



# now examind statistical difference in expression levels
data_for_tukey <- DEG_of_interest[DEG_of_interest$process=="SA/SAR",]
for (treat in unique(data_for_tukey$Treatment)){
  for (gene in unique(data_for_tukey$gene)){
    minimal <- min(data_for_tukey$count[data_for_tukey$Treatment==treat & data_for_tukey$gene==gene])
    data_for_tukey <- data_for_tukey[!data_for_tukey$count==minimal,]
  }
}


require(agricolae)
fit <- aov(formula = count ~ Treatment + time_point , data = data_for_tukey_final[data_for_tukey_final$gene=="PR5",])

SNK_res <- SNK.test(fit, "Treatment",alpha = 0.05)
SNK_res$groups

TukeyHSD(x = fit)


##### DEGs analysis starts here
all_pairwise <- function(treatments=c("#1", "#2", "#3"), control = "Control", p_adjust="qval",data=dds_res){
  final_results <- data.frame()
  
  for (treatment in treatments){
    res_pairwise <- results(dds_res, tidy = T, c("Treatment",treatment,control))
    res_pairwise$tested_pair <- paste(treatment,control)
    qobj <- qvalue(p = res_pairwise$pvalue)
    res_pairwise$qval <- qobj$qvalues
    final_results <- rbind(final_results,res_pairwise)
  }
  return(final_results)
}


final_deg_summary <- all_pairwise()
final_deg_summary <- final_deg_summary[!is.na(final_deg_summary$pvalue),]


unique(final_deg_summary$tested_pair)

# now analyze only CommenCom and MixedCom vs PathoCom
DEGs_final_all_vs_control <- function(data=final_deg_summary, adjustmet="qval",val=0.05, file=T, file_name="final_DEG_table_all_vs_Control", 
                                      dir="~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/"){
  
  final_deg_summary_sig <- final_deg_summary[final_deg_summary[,adjustmet] < val,]
  final_deg_summary_sig <- final_deg_summary_sig[!is.na(final_deg_summary_sig[,adjustmet]),]
  
  final_deg_summary_sig_1 <- final_deg_summary_sig[final_deg_summary_sig$tested_pair %in% c("#1 Control"),]
  final_deg_summary_sig_2 <- final_deg_summary_sig[final_deg_summary_sig$tested_pair %in% c("#2 Control"),]
  final_deg_summary_sig_3 <- final_deg_summary_sig[final_deg_summary_sig$tested_pair %in% c("#3 Control"),]
  
  upregulated_genes_1_vs_c <- final_deg_summary_sig_1$row[final_deg_summary_sig_1$log2FoldChange>1]
  downregulated_genes_1_vs_c <- final_deg_summary_sig_1$row[final_deg_summary_sig_1$log2FoldChange<1]
  
  upregulated_genes_2_vs_c <- final_deg_summary_sig_2$row[final_deg_summary_sig_2$log2FoldChange>1]
  downregulated_genes_2_vs_c <- final_deg_summary_sig_2$row[final_deg_summary_sig_2$log2FoldChange<1]
  
  upregulated_genes_3_vs_c <- final_deg_summary_sig_3$row[final_deg_summary_sig_3$log2FoldChange>1]
  downregulated_genes_3_vs_c <- final_deg_summary_sig_3$row[final_deg_summary_sig_3$log2FoldChange<1]

  

  if (length(upregulated_genes_2_vs_c)==0){
    upregulated_genes_2_vs_c <- NA
  }
  if (length(downregulated_genes_2_vs_c)==0){
    downregulated_genes_2_vs_c <- NA
  }
  if (length(upregulated_genes_3_vs_c)==0){
    upregulated_genes_3_vs_c <- NA
  }
  if (length(downregulated_genes_3_vs_c)==0){
    downregulated_genes_3_vs_c <- NA
  }
  if (length(upregulated_genes_1_vs_c)==0){
    upregulated_genes_1_vs_c <- NA
  }
  if (length(downregulated_genes_1_vs_c)==0){
    downregulated_genes_1_vs_c <- NA
  }
  
  
  tmp1 <- cbind("gene"=upregulated_genes_2_vs_c,"pair"="2_vs_c","regulation"="upregulated")
  tmp2 <- cbind("gene"=upregulated_genes_3_vs_c,"pair"="3_vs_c","regulation"="upregulated")
  tmp3 <- cbind("gene"=upregulated_genes_1_vs_c,"pair"="1_vs_c","regulation"="upregulated")
  tmp4 <- cbind("gene"=downregulated_genes_2_vs_c,"pair"="2_vs_c","regulation"="downregulated")
  tmp5 <- cbind("gene"=downregulated_genes_3_vs_c,"pair"="3_vs_c","regulation"="downregulated")
  tmp6 <- cbind("gene"=downregulated_genes_1_vs_c,"pair"="1_vs_c","regulation"="downregulated")
  
  final_DEG_vs_1 <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
  if (file==T){
    final_file_name <- paste(file_name,adjustmet,val,"csv",sep = ".")
    final_file_name_path <- paste(dir,final_file_name, sep = "")
    write.csv2(x = final_DEG_vs_1, file = final_file_name_path)
  }
  return(as.data.frame(final_DEG_vs_1))
}

DEGs<- DEGs_final_all_vs_control(adjustmet = "qval", val = 0.05)

sum(DEGs$pair=="1_vs_c")
sum(DEGs$pair=="2_vs_c")
sum(DEGs$pair=="3_vs_c")


# Expand the DEGs to intersections and unique
DEGs_list <- read.csv2("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/final_DEG_table_all_vs_Control.qval.0.05.csv")


treat2_3_up <- intersect(DEGs$gene[DEGs$pair=="3_vs_c" & DEGs$regulation=="upregulated"],DEGs$gene[DEGs$pair=="2_vs_c" & DEGs$regulation=="upregulated"])
treat2_3_down <- intersect(DEGs$gene[DEGs$pair=="3_vs_c" & DEGs$regulation=="downregulated"],DEGs$gene[DEGs$pair=="2_vs_c" & DEGs$regulation=="downregulated"])
treat2_3_down <- setdiff(treat2_3_down, DEGs$gene[DEGs$pair=="1_vs_c" & DEGs$regulation=="downregulated"])

DEGs_list$intersection[DEGs_list$gene %in% treat2_3_up] <- "treat_2_3_up"
DEGs_list$intersection[DEGs_list$gene %in% treat2_3_down] <- "treat_2_3_down"

#! write.csv(x = DEGs_list, file = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/final_DEG_table_all_vs_Control.qval.0.05.intersect.csv")

# proportional Venn diagram (Euler), requires manual insertion of gene amounts, as recivied by Venn Diagram results
library(eulerr)

up_1<- DEGs$gene[DEGs$pair=="1_vs_c" & DEGs$regulation=="upregulated"]
up_2<- DEGs$gene[DEGs$pair=="2_vs_c" & DEGs$regulation=="upregulated"]
up_3<- DEGs$gene[DEGs$pair=="3_vs_c" & DEGs$regulation=="upregulated"]




t1_2_up_intersect <- intersect(up_1,up_2)
t1_3_up_intersect <- intersect(up_1,up_3)
t2_3_up_intersect <- intersect(up_2,up_3)

A_B <- length(t1_2_up_intersect)
A_C <- length(t1_3_up_intersect)
B_C <- length(t2_3_up_intersect)


A <- (length(up_1)-length(union(t1_2_up_intersect,t1_3_up_intersect)))
B <-  (length(up_2)-length(union(t1_2_up_intersect,t2_3_up_intersect)))
C <-  (length(up_3)-length(union(t1_3_up_intersect,t2_3_up_intersect)))

A_B_C <- length(intersect(intersect(up_1,up_2),up_3))



one_contained <- euler(c("A" = A, "B" = B, "C" = C, 
                         "A&B" = A_B, "A&C" = A_C, "B&C" = B_C,
                         "A&B&C" = A_B_C),
                       shape = "ellipse")

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/Euler_diagrman_DEG_Lu_upregulated.pdf", useDingbats = F)
plot(one_contained, quantities = TRUE, 
     fill = c("#d95f02","#1b9e77","#7570b3"), labels = c("PathoCom", "CommenCom", "MixedCom"), alpha = 0.7)
#! dev.off()

# now the downregulated DEGs

down_1<- DEGs$gene[DEGs$pair=="1_vs_c" & DEGs$regulation=="downregulated"]
down_2<- DEGs$gene[DEGs$pair=="2_vs_c" & DEGs$regulation=="downregulated"]
down_3<- DEGs$gene[DEGs$pair=="3_vs_c" & DEGs$regulation=="downregulated"]

t1_2_down_intersect <- intersect(down_1,down_2)
t1_3_down_intersect <- intersect(down_1,down_3)
t2_3_down_intersect <- intersect(down_2,down_3)

A_B <- length(t1_2_down_intersect)
A_C <- length(t1_3_down_intersect)
B_C <- length(t2_3_down_intersect)

A <- (length(down_1)-length(union(t1_2_down_intersect,t1_3_down_intersect)))
B <-  (length(down_2)-length(union(t1_2_down_intersect,t2_3_down_intersect)))
C <-  (length(down_3)-length(union(t1_3_down_intersect,t2_3_down_intersect)))

A_B_C <- length(intersect(intersect(down_1,down_2),down_3))


one_contained <- euler(c("A" = A, "B" = B, "C" = C, 
                         "A&B" = A_B, "A&C" = A_C, "B&C" = B_C,
                         "A&B&C" = A_B_C),
                       shape = "ellipse")

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/Euler_diagrman_DEG_Lu_downregulated.pdf", useDingbats = F)
plot(one_contained, quantities = TRUE, 
     fill = c("#d95f02","#1b9e77","#7570b3"), labels = c("PathoCom", "CommenCom", "MixedCom"), alpha = 0.7)
#! dev.off()




## now heatmap of z-scores
library(DEFormats)
library("gplots")
library(edgeR)
library(limma)
library(EDASeq)
library(RColorBrewer)

raw_counts<- assay(dds)


final_metadata_table$Treatment_unique <- c("PathoCom.s1", "MixedCom.s1", "MixedCom.s2", 
                                           "Control.s1", "CommenCom.s1", "Control.s2",
                                           "CommenCom.s2", "CommenCom.s3", "MixedCom.s3", 
                                           "CommenCom.s4", "Control.s3", "MixedCom.s4", 
                                           "PathoCom.s2","Control.s4", "PathoCom.s3", "PathoCom.s4")

colnames(raw_counts) <- final_metadata_table$Treatment_unique[match(colnames(raw_counts), final_metadata_table$sample_file_name)]

dds_res_DGEobj <- as.DGEList(dds_res)


# work only with DEGs (in comparison to control)
DEGs_all_vs_control <- read.csv2("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/final_DEG_table_all_vs_Control.qval.0.05.csv")
View(DEGs_all_vs_control)
unique_DEGs_all_vs_control <- unique(DEGs_all_vs_control$gene)[DEGs_all_vs_control$regulation == "upregulated"]


ggplot(data = query_gene, (aes(x = Treatment, y = count)))+
  geom_point()

DEGs_mean <- read.csv("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/DEGs_normalized_by_sample.csv", row.names = 1)

colnames(DEGs_mean) <- c("Control.1","PathoCom.1","CommenCom.1", "MixedCom.1", "Control.2", "PathoCom.2",
                         "CommenCom.2", "MixedCom.2", "Control.3", "PathoCom.3", "CommenCom.3", "MixedCom.3",
                         "Control.4", "PathoCom.4", "CommenCom.4", "MixedCom.4")


# First heatmap of the euclidean distance between samples
#! pdf(file = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/heatmap_euc_distance_DEGs_by_sample.pdf", useDingbats = F)
heatmap.2(as.matrix(dist(t(DEGs_mean))), trace = "none")
#! dev.off()

library(pheatmap)

#heatmap colors
library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

annotation_temp<- gsub(colnames(DEGs_mean), pattern = ".1", replacement = "")
annotation_temp<- gsub(annotation_temp, pattern = ".2", replacement = "")
annotation_temp<- gsub(annotation_temp, pattern = ".3", replacement = "")
annotation_temp<- gsub(annotation_temp, pattern = ".4", replacement = "")





annotation <- data.frame(Var1 = factor(x = annotation_temp))
rownames(annotation) <- colnames(DEGs_mean)


mycolors <- c("#969696",	"#d95f02",	"#1b9e77", "#7570b3")
names(mycolors) <- unique(annotation$Var1)
mycolors <- list(Var1 = mycolors)


#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/heatmap_DEGs_by_sample_scaled.pdf", useDingbats = F)
pheatmap(as.matrix(DEGs_mean), scale = "row", col=rev(morecols(50)), labels_row=NULL, 
         show_rownames=F, annotation = annotation, annotation_colors = mycolors)
#! dev.off()

data <- as.matrix(DEGs_mean)[ , order( colnames( as.matrix(DEGs_mean) ) ) ]
col.order <- c("MixedCom.1", "MixedCom.3", "MixedCom.2", "MixedCom.4", "CommenCom.4", "CommenCom.2", "CommenCom.3", "CommenCom.1", "PathoCom.4", "PathoCom.2", "PathoCom.1", "PathoCom.3", "Control.4", "Control.1", "Control.2", "Control.3")
data <- data[,col.order]

rownames(data[h$tree_row[["order"]],])
temp <- sort(cutree(h$tree_row, k=6))
table(temp)
temp_remove <- c(names(temp[temp==4]), names(temp[temp==2]), names(temp[temp==6]))





data <- read.csv("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/DEGs_normalized_by_sample_corrected.csv", row.names = 1)
data <- as.matrix(data)
          
            
#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/heatmap_DEGs_by_sample_no_cluster_scaled_corrected.pdf", useDingbats = F)
h<- pheatmap(data[!(rownames(data) %in% temp_remove),], scale = "row", col=rev(morecols(50)), labels_row=NULL, cluster_cols=FALSE,
         show_rownames=F, annotation = annotation, annotation_colors = mycolors, Colv=FALSE, dendrogram="row")
#! dev.off()



DEGs_DF_average <- read.csv(file = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/DEGs_DF_average.csv", row.names = 1, sep = ";")


DEGs_DF_average$treatment <- factor(DEGs_DF_average$treatment, levels = c("Control", "PathoCom", "CommenCom", "MixedCom"))
DEGs_DF_average$pair <- factor(DEGs_DF_average$pair, levels = c("PathoCom set", "CommenCom set", "MixedCom set"))


#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/DEGs_by_set_by_treatment_high_bar.pdf", useDingbats = F)
ggplot(DEGs_DF_average, aes(x = treatment, y = scaled_counts, color = treatment, group = treatment))+
  geom_point() +
  facet_grid(regulation ~ pair) +
  scale_color_manual(values = c("#969696",	"#d95f02",	"#1b9e77", "#7570b3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim = c(-0.75,4.4))
#! dev.off()



dds_cpm <- edgeR::cpm(y =  raw_counts)

logcounts <- edgeR::cpm(dds_DGEList,log=TRUE)

groups <- final_metadata_table$Treatment[match(colnames(raw_counts), final_metadata_table$Treatment_unique)]
time_point <- final_metadata_table$Time_point[match(colnames(raw_counts), final_metadata_table$Treatment_unique)]



dds_DGEList <- DGEList(counts = raw_counts, group = groups)

design <- model.matrix(~0 + groups)


v <- voomWithQualityWeights(counts = dds_DGEList, block = time_point, design = design)


par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")





#barplot for REVIGO enrichment summary

library(egg)

top15_bar_plot <- function(input_data, top_hits_number=10, file=F, filename="~/ownCloud/My papers/Syncoms_paper/Final_edited_figures/figure5/3_and_2_vs_1_barplot_top10_qval0.05.pdf", width=8, height=8, full_description = T, shorten = F, title = T){
  
  
  if (shorten == T){
    if (length(input_data$term_ID) < 10){
      temp_full <- data.frame()
      for (i in 1:(10-(length(input_data$term_ID)))){
        temppy <- input_data[1,]
        temppy$description <- i
        temppy$log10.p.value_minus1 <- 0
        temppy$term_ID <- i
        temp_full <- rbind(temp_full, temppy)
      }
      input_data$term_ID <- as.character(input_data$term_ID)
      input_data$description <- as.character(input_data$description)
      input_data <- rbind(input_data, temp_full)
      input_data$term_ID <- as.factor(input_data$term_ID)
      input_data$description <- as.factor(input_data$description)
      
    }
  }

  

  
  if (full_description == T){
    input_data$description <- paste(paste(input_data$description, input_data$term_ID, sep = "\n ("), ")", sep = "")
  }
  
  
  levels_ordered <- as.character(input_data$description[order(input_data$log10.p.value_minus1, decreasing = T)])
  input_data$description <- factor(input_data$description, levels = levels_ordered)
  top <- levels(input_data$description)[1:top_hits_number]
  input_data_top <- input_data[input_data$description %in% top,]
  input_data_top$color <- "#347f99"
  
  if (title == T){
    g <- ggplot(data = input_data_top , aes(x = log10.p.value_minus1, y = description, fill = color))+
      geom_col() + 
      theme_classic() +
      #theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme(text = element_text(size=20)) +
      xlab("Enrichment\n(-log10pvalue)") +
      coord_cartesian(xlim=c(0,10)) + 
      theme(axis.title.x = element_blank()) +
      scale_fill_manual(values = "#347f99") +
      theme(legend.position = "none")
    
  } else {
    
    g <- ggplot(data = input_data_top , aes(x = log10.p.value_minus1, y = description, fill = color))+
      geom_col() + 
      theme_classic() +
      #theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme(text = element_text(size=20)) +
      xlab("Enrichment\n(-log10pvalue)") +
      coord_cartesian(xlim=c(0,10)) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank()) +
      scale_fill_manual(values = "#347f99")+
      theme(legend.position = "none")
    
  }
  

  if(file==F){
    return(g)
  } else {
    pdf(file=filename, width = width, height = height, useDingbats = F)
    print(g)
    dev.off()
  }
}


# MixedCom and CommenCom intersection
summary_REVIGO <- read.csv2(file = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/REVIGO/CommenCom_MixedCom_REVIGO.csv", sep = ",")
summary_REVIGO$log10.p.value
summary_REVIGO$log10.p.value <- as.numeric(as.character(summary_REVIGO$log10.p.value))
summary_REVIGO$log10.p.value_minus1<- summary_REVIGO$log10.p.value*-1
#top15_bar_plot(input_data = summary_REVIGO[summary_REVIGO$eliminated == 0,], top_hits_number = 10, file=T, filename = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/bar_plot_MixedCom_intersect_CommenCom_up.pdf.pdf")
g1 <- top15_bar_plot(input_data = summary_REVIGO[summary_REVIGO$eliminated == 0,], top_hits_number = 10, file=F, full_description = T, title = F)


# MixedCom unique
summary_REVIGO <- read.csv2(file = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/REVIGO/MixedCom_unique_REVIGO.csv", sep = ",")
summary_REVIGO$log10.p.value <- as.numeric(as.character(summary_REVIGO$log10.p.value))
summary_REVIGO$log10.p.value_minus1<- summary_REVIGO$log10.p.value*-1
#top15_bar_plot(input_data = summary_REVIGO[summary_REVIGO$eliminated == 0,], top_hits_number = 10, file=T, filename = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/bar_plot_MixedCom_up.pdf.pdf")
g2 <- top15_bar_plot(input_data = summary_REVIGO[summary_REVIGO$eliminated == 0,], top_hits_number = 10, file=F, full_description = T)

# CommenCom unique
summary_REVIGO <- read.csv2(file = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/REVIGO/CommenCom_unique_REVIGO.csv", sep = ",")
summary_REVIGO$log10.p.value <- as.numeric(as.character(summary_REVIGO$log10.p.value))
summary_REVIGO$log10.p.value_minus1<- summary_REVIGO$log10.p.value*-1
#top15_bar_plot(input_data = summary_REVIGO[summary_REVIGO$eliminated == 0,], top_hits_number = 10, file=T, filename = "~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/bar_plot_CommenCom_up.pdf.pdf")
g3 <- top15_bar_plot(input_data = summary_REVIGO[summary_REVIGO$eliminated == 0,], top_hits_number = 10, file=F, full_description = T, title = F)


single_bar_plot <- ggarrange(g1,g3,g2, heights = c(1,0.2,1))

#! pdf("~/ownCloud/My papers/Syncoms_paper/Figures/RNA-seq/bar_plot_all_together_with_labeles.pdf", useDingbats = F, height = 20, width = 15)
print(single_bar_plot)
#! dev.off()

#now annotations
annotations_TAIR10 <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/RNA-seq_Exp_July_2019_rep2/results/ATH_GO_GOSLIM_updated.txt", sep= "\t", header = F)
colnames(annotations_TAIR10) <- c("locus name",	"TAIR accession",	"object name",	"relationship type",	"GO term",	"GO ID",	"TAIR Keyword", "ID	Aspect",	"GOslim term",	"Evidence code",	"Evidence description",	"Evidence with",	'Reference',	"Annotator",	"Date annotated")


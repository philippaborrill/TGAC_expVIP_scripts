# Plot gene expression data for manuscript
# Philippa Borrill
# 18.7.2016


setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\Wheat_Expression")


##get data as tpm
tpm.data<-read.table("edited_final_output_tpm.txt",sep="\t",header=T)  

head(tpm.data)
#summary(tpm.data)
dim(tpm.data)


# it didn't correctly read in the headers therefore adjust for this:
rownames(tpm.data) <- tpm.data[,1]
tpm <- tpm.data[,-1]
head(tpm)
dim(tpm)
head(rownames(tpm))
head(colnames(tpm))

# Expressed genes:

# Genes expressed in any sample (>1 read mapping)
# genes with no expression get a TRUE
noexpr <- (rowSums(tpm)) == 0
#count number of true (ie genes not expressed)
sum(noexpr)
# number of genes expressed 
273739-sum(noexpr)
(273739-sum(noexpr))/273739

# want to know how many genes are expressed at least tpm > 2 in one sample
# get max number from each row
rowmax <- apply(tpm,1,max)
# make logical vector whether that max value is over 2
rowmax_over2 <- rowmax>2
head(rowmax_over2)
#count up all which were TRUE (ie over 2)
sum(rowmax_over2)
sum(rowmax_over2)/273739

# count number of genes expressed per sample with tpm > 2 (arbitrary cut off)
genes_really_expr_per_sample <- colSums(tpm>2)
mean(genes_really_expr_per_sample)
mean(genes_really_expr_per_sample)/273739

# find sample with min genes expressed
genes_really_expr_per_sample[grep(min(genes_really_expr_per_sample),genes_really_expr_per_sample)]

# find sample with max genes expressed
genes_really_expr_per_sample[grep(max(genes_really_expr_per_sample),genes_really_expr_per_sample)]

####### want to find out the numbers of genes with high, medium and low confidence which are expressed at >2 tpm #####
# first read in info about gene confidence level
transcript_confidence_level <-read.csv(file= "Y:\\expression_browser\\TGAC_assembly\\transcript_confidence_level.csv", header=T)
head(transcript_confidence_level)

# we already have whether each gene is expressed over 2 tpm (rowmax_over2) but need to convert to data.frame
rowmax_over2_df <- as.data.frame(rowmax_over2)
head(rowmax_over2_df)
rownames(rowmax_over2_df)

# merge the confidence level with TRUE/FALSe of gene expressed over 2 tpm
merged_confidence <- merge(transcript_confidence_level, rowmax_over2_df,by.x = "transcript_ID", by.y= 0)
head(merged_confidence)
summary(merged_confidence)

# number of total genes 
num_total_genes <- 273739


# extract just low confidence genes
low_confidence <- merged_confidence[merged_confidence$gene_confidence=="Low",]
head(low_confidence)
# count number of low confidence genes
length(low_confidence$transcript_ID)

# % low confidence transcripts
length(low_confidence$transcript_ID)/num_total_genes

# only select low confidence genes with expr > 2 tpm
low_confidence_over2tpm <- (low_confidence[low_confidence$rowmax_over2=="TRUE",])
head(low_confidence_over2tpm)
# count number of low confidence genes with expr > 2 tpm
length(low_confidence_over2tpm$transcript_ID)
# % of low confidence genes with expr > 2 tpm
length(low_confidence_over2tpm$transcript_ID)/length(low_confidence$transcript_ID)

# extract just high confidence genes
high_confidence <- merged_confidence[merged_confidence$gene_confidence=="High",]
head(high_confidence)
#count number of high confidence genes
length(high_confidence$transcript_ID)

# % high confidence transcripts
length(high_confidence$transcript_ID)/num_total_genes

# only select high confidence genes with expr >2 tpm
high_confidence_over2tpm <- (high_confidence[high_confidence$rowmax_over2=="TRUE",])
head(high_confidence_over2tpm)
# count number of high confidence genes with expr > 2tpm
length(high_confidence_over2tpm$transcript_ID)

# % high conf genes with >2 tpm
length(high_confidence_over2tpm$transcript_ID)/length(high_confidence$transcript_ID)



################# Figure 1 ########################
#### want to make some nices graphs with expression levels shown by metadata group
##### first need to make a nice dataframe to work with

# get genes "really expressed" tpm
genes_really_expr_per_sample_tpm <- colSums(tpm>2)

#convert to dataframes
df.genes_really_expr_per_sample_tpm <- as.data.frame(genes_really_expr_per_sample_tpm)
head(df.genes_really_expr_per_sample_tpm)
rownames(df.genes_really_expr_per_sample_tpm)

# read in metadata
# NB had to make a copy of the edited_metadata_output.txt which was save in notepad++ with UNIX EOL.
meta_data<-read.table("Y:\\expression_browser\\kallisto_analysis\\edited_metadata_output2.txt",sep="\t",header=T)
dim(meta_data)
colnames(meta_data)
rownames(meta_data) <- meta_data[,2]
rownames(meta_data)


# combine expression data with metadata
expr_with_metadata <- merge(df.genes_really_expr_per_sample_tpm, meta_data, by = "row.names")
dim(meta_data)
dim(df.genes_really_expr_per_sample_tpm)
dim(expr_with_metadata)
rownames(expr_with_metadata)<- expr_with_metadata[,1]
expr_with_metadata <- expr_with_metadata[,-1]
head(expr_with_metadata)
dim(expr_with_metadata)
colnames(expr_with_metadata)

# expr_with_metadata is complete dataset

#plot graph
install.packages("ggplot2")
library("ggplot2")
# want to make some nice graphs - use ggplot2
options(scipen=5)

colnames(expr_with_metadata)
data(expr_with_metadata)
theme_set(theme_grey())

jpeg(file='genes_expr_tpm_over2_vs_mapped_reads_per_study.jpg', height=500, width=800)
par(mar=c(4,4,4,2)+0.1,mgp=c(3,1,0))
qplot(Mapped.reads, genes_really_expr_per_sample_tpm, data=expr_with_metadata, colour = secondary_study_accession) +
  scale_color_manual(values=c("light blue","blue","violet", "purple", "yellow", "orange", "red","dark red", "light green", "green", "dark green", "brown", "black", "grey", "pink", "turquoise")) + ylim(c(0,55000))
dev.off()

postscript(file='genes_expr_tpm_over2_vs_mapped_reads_per_study.eps', height=500, width=800)
par(mar=c(4,4,4,2)+0.1,mgp=c(3,1,0))
qplot(Mapped.reads, genes_really_expr_per_sample_tpm, data=expr_with_metadata, colour = secondary_study_accession) + 
  scale_color_manual(values=c("light blue","blue","violet", "purple", "yellow", "orange", "red","dark red", "light green", "green", "dark green", "brown", "black", "grey", "pink", "turquoise")) + ylim(c(0,55000))
dev.off()

# basic graphics tpm vs mapped reads (just to check ggplot looks ok)
jpeg(file='genes_expressed_over_2tpm_vs_mapped_reads.jpg', height=500, width=500)
par(mar=c(4,4,4,2)+0.1,mgp=c(3,1,0))
plot((expr_with_metadata[,15]),expr_with_metadata[,1], xlab= "mapped reads (millions)", ylab="number of expressed genes tpm > 2", 
     pch = 19, col="red", ylim=c(0,50000))
#abline(lm(expr_with_metadata[,2] ~ expr_with_metadata[,4]), col="blue")
dev.off()

# if want to get R-squared value
summary(lm(expr_with_metadata[,15] ~ expr_with_metadata[,1]))



########## Figure 2 #####################
### look at grouped data

setwd("Y:\\expression_browser\\kallisto_analysis\\kallisto_results_edited")

tpm_grouped <- read.table("mean_tpm_per_group_high_level_tissue_no_nullitetra.txt", header=T, sep="\t")
dim(tpm_grouped)
rownames(tpm_grouped) <- tpm_grouped[,1]
tpm_grouped <- tpm_grouped[,-1]
head(tpm_grouped)

#convert to matrix
tpm_grouped_matrix <- as.matrix(tpm_grouped)
dim(tpm_grouped_matrix)


# want to select only genes which are expressed over tpm 2 in at least 1 sample

genes_really_expr_grouped_data <- tpm_grouped_matrix[(rowSums(tpm_grouped_matrix>2)>0),]
head(genes_really_expr_grouped_data)
dim(genes_really_expr_grouped_data)
length(genes_really_expr_grouped_data)
nrow(genes_really_expr_grouped_data)
#select random rows from matrix (R can't plot a heatmap of all genes because there are too many)
# we have 85272 genes which are expressed over 2 tpm in at least one sample, want to select 1000 to plot

set.seed(30)
random_genes_tpm_grouped_matrix <- (genes_really_expr_grouped_data[(floor(runif(1000,min=1, max= (nrow(genes_really_expr_grouped_data)+1)))),])
dim(random_genes_tpm_grouped_matrix)
head(random_genes_tpm_grouped_matrix)
write.csv(rownames(random_genes_tpm_grouped_matrix),"1000_random_genes_Figure2.csv")

# make heatmap
library(RColorBrewer)
library(gplots)
### this is the graph for the paper

postscript(file="heatmap_1000_random_genes_expressed_tpm_over2_in_at_least_1_sample_grouped_data2.eps", height=10, width=10)
par(mar=c(1,1,1,1)+0.1,mgp=c(6,1,0))
heatmap.2(random_genes_tpm_grouped_matrix, col=rev(heat.colors(75)), Rowv=TRUE, Colv=TRUE, dendrogram= "both", key=TRUE,
          keysize=0.5, cexCol=2, trace="none", na.rm =TRUE, density.info="none", scale="row",margins = c(15,15), main= "Similarity between expression of genes between tissues (1000 random expressed genes)")
dev.off()

jpeg(file="heatmap_1000_random_genes_expressed_tpm_over2_in_at_least_1_sample_grouped_data2.jpg", height=1000, width=5000)
par(mar=c(1,1,1,1)+0.1,mgp=c(6,1,0))
heatmap.2(random_genes_tpm_grouped_matrix, col=rev(heat.colors(75)), Rowv=F, Colv=TRUE, dendrogram= "column", key=TRUE,
          keysize=0.5, cexCol=2, trace="none", na.rm =TRUE, density.info="none", scale="row",margins = c(15,15), main= "Similarity between expression of genes between tissues (1000 random expressed genes)")
dev.off()

# want to see if genes most expressed make sense
# take grain samples:

grain_samples_genes_really_expr_grouped_data <- genes_really_expr_grouped_data[,grep("grain",colnames(genes_really_expr_grouped_data))]
dim(grain_samples_genes_really_expr_grouped_data)
average_grain_expr<- rowMeans(grain_samples_genes_really_expr_grouped_data)
head(average_grain_expr)
sorted_grain_expr<- sort(average_grain_expr)
tail(sorted_grain_expr,20)
write.csv(tail(sorted_grain_expr,10),file="top10_most_expr_genes_in_grain.csv")

# take leaf samples
colnames(genes_really_expr_grouped_data)
dim(genes_really_expr_grouped_data)
grep("lea",colnames(genes_really_expr_grouped_data))

leaf_samples_genes_really_expr_grouped_data <- genes_really_expr_grouped_data[,grep("leaves",colnames(genes_really_expr_grouped_data))]
dim(leaf_samples_genes_really_expr_grouped_data)

average_leaf_expr<- rowMeans(leaf_samples_genes_really_expr_grouped_data)
head(average_leaf_expr)
sorted_leaf_expr<- sort(average_leaf_expr)
tail(sorted_leaf_expr,20)
write.csv(tail(sorted_leaf_expr,10),file="top10_most_expr_genes_in_leaf.csv")

# take root samples ## don't use this in paper doesn't look great!
grep("root",colnames(genes_really_expr_grouped_data))

root_samples_genes_really_expr_grouped_data <- genes_really_expr_grouped_data[,grep("root",colnames(genes_really_expr_grouped_data))]
dim(root_samples_genes_really_expr_grouped_data)

average_root_expr<- rowMeans(root_samples_genes_really_expr_grouped_data)
head(average_root_expr)
sorted_root_expr<- sort(average_root_expr)
tail(sorted_root_expr,10)
write.csv(tail(sorted_root_expr,10),file="top10_most_expr_genes_in_root.csv")

########## try to plot sleuth analysis ###FAILED


files <- as.character(list.files(path="Y:\\expression_browser\\kallisto_analysis\\stress_seedlings_analysis\\separate_controls"), full.names=T)
files
files <- files[grep("condition*",files, perl = TRUE)]
files

setwd("Y:\\expression_browser\\kallisto_analysis\\stress_seedlings_analysis\\separate_controls")
sample <- read.csv("condition1 hour of drought stress.csv", header=TRUE)
head(sample)

qval_sample <- sample[sample$qval<0.05,]
head(qval_sample)

qval_sample <- qval_sample[,c("target_id", "qval", "b")]
head(qval_sample)

sample2 <- read.csv("condition1_hour_heat.csv", header=TRUE)
head(sample2)

qval_sample2 <- sample2[sample2$qval<0.05,]
head(qval_sample2)

qval_sample2 <- qval_sample2[,c("target_id", "qval", "b")]
head(qval_sample2)

merged <- merge.data.frame(qval_sample, qval_sample2, by="target_ID", all=TRUE)


new_sample <- sample[order(sample$target_id),]
head(new_sample)


### Look for genes which are expressed at similar levels across all tissues


#set co.var function
co.var <- function(x) ( 100*sd(x)/mean(x) )

setwd("Y:\\expression_browser\\kallisto_analysis\\kallisto_results_edited")

##get data as tpm
tpm.data<-read.table("final_output_tpm.txt",sep="\t",header=T)  

head(tpm.data)
#summary(tpm.data)
dim(tpm.data)


# it didn't correctly read in the headers therefore adjust for this:
rownames(tpm.data) <- tpm.data[,1]
tpm <- tpm.data[,-1]
head(tpm)

colnames(tpm)
tpm_excl_NT<- tpm[,c(1:300,398:418)]
dim(tpm_excl_NT)

#make new data frame expr_across_samples containing the mean, sd and co.var
expr_across_samples <- as.data.frame(apply(tpm_excl_NT,1,mean))
expr_across_samples[,2] <- as.data.frame(apply(tpm_excl_NT,1,sd))
expr_across_samples[,3] <- as.data.frame(apply(tpm_excl_NT,1,co.var))

head(expr_across_samples)
colnames(expr_across_samples) <- c("mean","sd","co.var")
head(expr_across_samples)
dim(expr_across_samples)

# try to select only genes where all samples have expr over 2
rowmin <- apply(tpm_excl_NT,1,min)
rowmin_over2 <- rowmin>2
expr_in_all_samples <- expr_across_samples[rowmin_over2,]
dim(expr_in_all_samples)
head(expr_in_all_samples)

# sort from smallest covar to largest
expr_in_all_samples_sorted <- expr_in_all_samples[order(expr_in_all_samples$co.var),]
head(expr_in_all_samples_sorted)
dim(expr_in_all_samples_sorted)
write.csv(expr_in_all_samples_sorted,"covar_genes_expr_over2tpm_all_samples.csv")

median(expr_in_all_samples_sorted$co.var)

# check if the stable expressed ref genes are still included (they aren't!)
stable.genes <- expr_in_all_samples_sorted[c("Traes_2DL_B9B3770C9.1", "Traes_5BS_4179BB1C5.1", "Traes_6AS_82E93D872.1", "Traes_2BS_7F29300C6.2", "Traes_5AL_8E706641D.1", "TRAES3BF019000090CFD_t1", "Traes_1DS_A7A785837.1", "Traes_1DL_9926BBAE1.1", "Traes_5DL_C3D011187.2", "Traes_1DL_A0E21D8B9.1", "Traes_3DL_2114C4621.2", "Traes_5BL_EBD050E21.2", "Traes_6DL_B3AD4834A.1"),]
stable.genes

# get expression of ref genes
stable.genes_all_genes <- expr_across_samples[c("Traes_2DL_B9B3770C9.1", "Traes_5BS_4179BB1C5.1", "Traes_6AS_82E93D872.1", "Traes_2BS_7F29300C6.2", "Traes_5AL_8E706641D.1", "TRAES3BF019000090CFD_t1", "Traes_1DS_A7A785837.1", "Traes_1DL_9926BBAE1.1", "Traes_5DL_C3D011187.2", "Traes_1DL_A0E21D8B9.1", "Traes_3DL_2114C4621.2", "Traes_5BL_EBD050E21.2", "Traes_6DL_B3AD4834A.1"),]
head(stable.genes_all_genes)
stable.genes_all_genes_sorted <- stable.genes_all_genes[order(stable.genes_all_genes$co.var),]
stable.genes_all_genes_sorted

write.csv(stable.genes_all_genes_sorted, "stable_genes_sorted.csv")

ref_gene_covar <- read.table("reference_genes_covar_for_R.txt", header=TRUE, sep = "\t")
head(ref_gene_covar)
dim(ref_gene_covar)

jpeg(file=paste("Coefficient of variance per gene using all 321 samples -excludes nullitetras with ref gene annotation.jpg", sep=""), height=1500, width=1500)
par(mar=c(6,10,4,3)+0.1,mgp=c(5,1,0))
hist(expr_in_all_samples_sorted[,3], #xlim=c(1,500), 
     xlab="Co-efficient of variance", ylab="Number of genes", breaks = 100, cex.axis=2, cex.lab = 2, cex.main=3, main= "Co-efficient of variance per gene in tpm")

points(ref_gene_covar[,2],c(410,390,370,350,330,310,290,270,250,230,210,190,170),pch=4, cex=2,lwd=2)
text(c(57,66,72,76,79,93,95,95,103,117,117,119,140),c(410,390,370,350,330,310,290,270,250,230,210,190,170),ref_gene_covar[,1], cex = 2,pos=4)
dev.off()

postscript(file=paste("Coefficient of variance per gene using all 321 samples -excludes nullitetras with ref gene annotation.eps", sep=""), height=750, width=750)
par(mar=c(6,10,4,3)+0.1,mgp=c(5,1,0))
hist(expr_in_all_samples_sorted[,3], #xlim=c(1,500), 
     xlab="Co-efficient of variance", ylab="Number of genes", breaks = 100, cex.axis=2, cex.lab = 2, cex.main=3, main= "Co-efficient of variance per gene in tpm")

points(ref_gene_covar[,2],c(410,390,370,350,330,310,290,270,250,230,210,190,170),pch=4, cex=2,lwd=2)
text(c(57,66,72,76,79,93,95,95,103,117,117,119,140),c(410,390,370,350,330,310,290,270,250,230,210,190,170),ref_gene_covar[,1], cex = 1.5,pos=4)
dev.off()


# make line graphs of most consistently expressed and ref genes

ref_gene_expr <- tpm_excl_NT[c("Traes_2DL_B9B3770C9.1", "Traes_5BS_4179BB1C5.1", "Traes_6AS_82E93D872.1", "Traes_2BS_7F29300C6.2", "Traes_5AL_8E706641D.1", "TRAES3BF019000090CFD_t1", "Traes_1DS_A7A785837.1", "Traes_1DL_9926BBAE1.1", "Traes_5DL_C3D011187.2", "Traes_1DL_A0E21D8B9.1", "Traes_3DL_2114C4621.2", "Traes_5BL_EBD050E21.2", "Traes_6DL_B3AD4834A.1"),]
head(ref_gene_expr)
dim(ref_gene_expr)
rownames(ref_gene_expr)
head(rowMeans(ref_gene_expr), n=5)

#make ref gene expr relative to average expr across all tissues
ref_gene_expr<- ref_gene_expr/rowMeans(ref_gene_expr)

install.packages(reshape)
library("reshape2")
head(ref_gene_expr)
dim(ref_gene_expr)
#assign ID to melt by
ref_gene_expr$id <- rownames(ref_gene_expr)
#melt data to re-arrange for ggplot
melted<-melt(ref_gene_expr)
head(melted,20)
dim(melted)

#get most fold difference
tail(melted[(order(melted$value)),])

t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  axis.text.x = element_text(angle = 90, hjust = 1),
  axis.ticks=element_blank()
)
postscript(file="Reference_gene_expression_across_321_samples.eps", height=750, width=750)
par(mar=c(6,10,4,3)+0.1,mgp=c(5,1,0))
ggplot(melted,aes(x=variable,y=value,colour=id, group=id))+ geom_point() + geom_line()+labs(y= "Relate expression", x= NULL) + ylim(0,8) +t1
dev.off()

## for consistently expr

stable_genes_top20 <-tpm_excl_NT[c(head(rownames(expr_in_all_samples_sorted),20)),]
stable_genes_top20<- stable_genes_top20/rowMeans(stable_genes_top20)

library("reshape2")
head(stable_genes_top20)
dim(stable_genes_top20)
#assign ID to melt by
stable_genes_top20$id <- rownames(stable_genes_top20)
#melt data to re-arrange for ggplot
melted_stable<-melt(stable_genes_top20)
head(melted_stable,20)
dim(melted_stable)

#get most fold difference
tail(melted_stable[(order(melted_stable$value)),])

postscript(file="Top20_stable_gene_expression_across_321_samples.eps", height=750, width=750)
par(mar=c(6,10,4,3)+0.1,mgp=c(5,1,0))
ggplot(melted_stable,aes(x=variable,y=value,colour=id, group=id))+ geom_point() + geom_line()+labs(y= "Relative expression", x= NULL) + ylim(0,8) +t1
dev.off()


# plot heatmap of stable genes and ref genes previously known:

# new stable genes:
stable_genes_top20 <-tpm_excl_NT[c(head(rownames(expr_in_all_samples_sorted),20)),]
dim(stable_genes_top20)

# known ref genes:
ref_gene_expr <- tpm_excl_NT[c("Traes_2DL_B9B3770C9.1", "Traes_5BS_4179BB1C5.1", "Traes_6AS_82E93D872.1", "Traes_2BS_7F29300C6.2", "Traes_5AL_8E706641D.1", "TRAES3BF019000090CFD_t1", "Traes_1DS_A7A785837.1", "Traes_1DL_9926BBAE1.1", "Traes_5DL_C3D011187.2", "Traes_1DL_A0E21D8B9.1", "Traes_3DL_2114C4621.2", "Traes_5BL_EBD050E21.2", "Traes_6DL_B3AD4834A.1"),]
dim(ref_gene_expr)

genes_for_heatmap <- rbind(ref_gene_expr,stable_genes_top20)
head(genes_for_heatmap)
dim(genes_for_heatmap)

matrix_genes_for_heatmap <- as.matrix(genes_for_heatmap)

# draws heatmap
library(RColorBrewer)
library(gplots)
jpeg(file="known ref genes plus stable genes.jpg", height=1000, width=1000)
par(mar=c(20,11,4,3)+0.1,mgp=c(6,1,0))

heatmap.2(matrix_genes_for_heatmap, col=rev(heat.colors(75)), Rowv=FALSE, Colv=FALSE, dendrogram= "none", key=TRUE,
          keysize=0.5,trace="none",  density.info="none", scale="row",margins = c(15,10), main= "Ref genes (top 11 rows) and 20 stable genes")
dev.off()

rownames(genes_for_heatmap)
## plot seedling stress

setwd("Y:\\expression_browser\\kallisto_analysis\\stress_seedlings_analysis\\separate_controls")
sample <- read.csv("genes_expr_10_conditions.csv", header=TRUE)

head(sample)

rownames(sample) <- sample[,1]
sample <- sample[,-1]
head(sample)

sample_matrix <- as.matrix(sample)
head(sample_matrix)

library(RColorBrewer)
library(gplots)

postscript(file="genes_expr_10_conditions_heatmap.eps", height=100, width=100)
par(mar=c(15,15,15,15)+0.1,mgp=c(1,1,0))

heatmap.2(sample_matrix, col=rev(heat.colors(75)), Rowv=T, Colv=TRUE, dendrogram= "both", key=TRUE,
          keysize=0.5,trace="none",  density.info="none", scale="row",cexRow =0.8, margins = c(15,15), main= "Genes expressed in 10 conditions", na.color ="grey")
dev.off()


### plot phosphate rice

setwd("Y:\\expression_browser\\kallisto_analysis\\rice_phosphate\\RAP")
data <- read.csv("shoots_fold_change_two_studies.csv")
head(data)

rownames(data) <- data[,1]
data <- data[,-1]
head(data)

postscript(file="scatterplot_shoots.eps", height=1000, width=1000)
par(mar=c(1,1,1,1)+0.1,mgp=c(6,1,0))
plot(data[,1],data[,2], ylim=c(-4,8))
abline(lm(data[,2]~data[,1]), col="black", lty=2)
dev.off()

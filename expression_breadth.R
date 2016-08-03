# Plot graphs of gene expression according to cM position for TGAC genome paper
# Philippa Borrill
# 02.08.2016

# Want to plot expression in all samples across the chromosomes 
# to see if there is positional variation in expression breadth (number of conditions in which genes are expressed)

# set working directory where the tpm file is

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

# want to exclude genes with expression <2 tpm across all samples

# get max number from each row
rowmax <- apply(tpm,1,max)
# make logical vector whether that max value is over 2
rowmax_over2 <- rowmax>2
head(rowmax_over2,20)
#check how many genes are expressed over 2 tpm
sum(rowmax_over2)
length(rowmax_over2)
# keep only genes where the max tpm was over 2
tpm <- tpm[rowmax_over2,]
head(rownames(tpm),20)
dim(tpm)

# need to put a true or false in each sample whether a gene is expressed to over 2 tpm
tpm_over_2 <- tpm >2
head(tpm_over_2)
num_cond_gene_expr <- as.data.frame(apply(tpm_over_2,1,sum))
head((num_cond_gene_expr))
head(rownames(num_cond_gene_expr))
colnames(num_cond_gene_expr) <- c("num_cond_expr")
head(num_cond_gene_expr)

# remove tpm.data which I don't need
rm(tpm.data)

# load in the cM position information
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\ordered_genes")

cm.data<-read.table("combined_transcripts_and_cM_position_with_header_sorted.txt",sep="\t",header=T)  
head(cm.data)
dim(cm.data)

# want to make the rownames of the cm.data the transcript ID
rownames(cm.data) <- cm.data[,2]
cm <- cm.data[,-2]
head(cm)
dim(cm)

# merge together the expression breadth with the cM position info
tpm_cm <- merge(num_cond_gene_expr, cm, by=0) # by=0 means by row.names
head(tpm_cm)
dim(tpm_cm)


# try calculating median expression breadth per genetic bin and plotting 

# first make a new column with the chrom_cM to use as a factor
tpm_cm_uniq_pos <- (transform(tpm_cm,chr_cM=paste0(chromosome,"_",cM)))
head(tpm_cm_uniq_pos)

# now try averaging per cM bin
agg_conditions_cm_median <- aggregate(tpm_cm_uniq_pos[,2], by = list(tpm_cm_uniq_pos$chr_cM), FUN=median)
head(agg_conditions_cm_median)
dim(agg_conditions_cm_median)

agg_conditions_cm_mean <- aggregate(tpm_cm_uniq_pos[,2], by = list(tpm_cm_uniq_pos$chr_cM), FUN=mean)
head(agg_conditions_cm_mean)
dim(agg_conditions_cm_mean)

# count genes per bin
agg_conditions_cm_count <- aggregate(tpm_cm_uniq_pos[,2], by = list(tpm_cm_uniq_pos$chr_cM), FUN=length)
head(agg_conditions_cm_count)
dim(agg_conditions_cm_count)

# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_conditions_cm_median$Group.1),'_') 
head(split_cols)
agg_conditions_cm_median <- data.frame(agg_conditions_cm_median[,2], do.call(rbind, split_cols))
colnames(agg_conditions_cm_median)
head(agg_conditions_cm_median)
colnames(agg_conditions_cm_median) <- c("conditions_per_gene", "chromosome", "cM")
head(agg_conditions_cm_median)
dim(agg_conditions_cm_median)
agg_conditions_cm_median <-as.data.frame (agg_conditions_cm_median)

# need to convert cM to "numeric" to get scaling to work on plot
agg_conditions_cm_median$cM=as.numeric(levels(agg_conditions_cm_median$cM))[agg_conditions_cm_median$cM]


# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_conditions_cm_mean$Group.1),'_') 
head(split_cols)
agg_conditions_cm_mean <- data.frame(agg_conditions_cm_mean[,2], do.call(rbind, split_cols))
colnames(agg_conditions_cm_mean)
head(agg_conditions_cm_mean)
colnames(agg_conditions_cm_mean) <- c("conditions_per_gene", "chromosome", "cM")
head(agg_conditions_cm_mean)
dim(agg_conditions_cm_mean)
agg_conditions_cm_mean <-as.data.frame (agg_conditions_cm_mean)

# need to convert cM to "numeric" to get scaling to work on plot
agg_conditions_cm_mean$cM=as.numeric(levels(agg_conditions_cm_mean$cM))[agg_conditions_cm_mean$cM]


# read in file which has centromere approx position
centromere_pos<-read.table("get_centromere_position.txt",sep="\t")  
head(centromere_pos)
dim(centromere_pos)

colnames(centromere_pos) <- c("chromosome","centromere_pos")
head(centromere_pos)

library(ggplot2)

# try plotting graph
pdf(file="All_samples_cM_excl_under_2_tpm_genes_median_expression_breadth_per_bin_with_centromere.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(agg_conditions_cm_median, aes(cM,conditions_per_gene)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values="black") +
  scale_x_continuous(breaks=c(0,50,100, 150,200))  + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()


# try plotting graph
pdf(file="All_samples_cM_excl_under_2_tpm_genes_mean_expression_breadth_per_bin_with_centromere.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(agg_conditions_cm_mean, aes(cM,conditions_per_gene)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values="black") +
  scale_x_continuous(breaks=c(0,50,100, 150,200))  + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()


### try replotting only for bins with > 3 genes #####

head(agg_conditions_cm_count)
agg_conditions_cm_count <- agg_conditions_cm_count[agg_conditions_cm_count$x >2 , ]
head(agg_conditions_cm_count)
dim(agg_conditions_cm_count)

agg_conditions_cm_median <- aggregate(tpm_cm_uniq_pos[,2], by = list(tpm_cm_uniq_pos$chr_cM), FUN=median)
head(agg_conditions_cm_median)
dim(agg_conditions_cm_median)

head(agg_conditions_cm_median)
agg_conditions_cm_median <- merge(agg_conditions_cm_median, agg_conditions_cm_count, by= "Group.1")
head(agg_conditions_cm_median)
agg_conditions_cm_median <- agg_conditions_cm_median[1:2]
head(agg_conditions_cm_median)

# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_conditions_cm_median$Group.1),'_') 
head(split_cols)
agg_conditions_cm_median <- data.frame(agg_conditions_cm_median[,2], do.call(rbind, split_cols))
colnames(agg_conditions_cm_median)
head(agg_conditions_cm_median)
colnames(agg_conditions_cm_median) <- c("conditions_per_gene", "chromosome", "cM")
head(agg_conditions_cm_median)
dim(agg_conditions_cm_median)
agg_conditions_cm_median <-as.data.frame (agg_conditions_cm_median)

# need to convert cM to "numeric" to get scaling to work on plot
agg_conditions_cm_median$cM=as.numeric(levels(agg_conditions_cm_median$cM))[agg_conditions_cm_median$cM]

library(ggplot2)

# try plotting graph
pdf(file="All_samples_cM_excl_under_2_tpm_genes_median_expression_breadth_per_bin_with_centromere_3genes_per_bin.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(agg_conditions_cm_median, aes(cM,conditions_per_gene)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values="black") +
  scale_x_continuous(breaks=c(0,50,100, 150,200))  + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()



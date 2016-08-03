# GOSEQ analysis of triads for seedling stress experiment for TGAC genome paper
# Philippa Borrill
# 14.07.2016

## README ##
# uses the output of the scripts:
# mergeing_and_analysing_sleuth_results_seedlings.R



#########GOSEQ analysis #################

# requires as input:

# expr data; named vector of gene names and differentially expressed genes (1 for diff expr, 0 for not diff expr)
# lengths of genes: numeric vector, same length as main vector, each entry gives length of corresponding gene in bp, if data unavailable for some genes set to NA
# category mapping: either data.frame of 1st column gene name, 2nd column category name (GO) (will have many rows per gene)
# OR list, where names of list entries are gene ID and entries themselves are a vector of category names

## want to do 2 main comparisons
# TRIADS vs WHOLE GENOME
# CONSERVED OR DIVERGENT TRIAD GENES vs TRIAD GENES

# for each of these comparisons I will need to make separate inputs as described above 
# i.e. for 1st comparison I will need the WHOLE GENOME as background
# for 2nd comparison I will need TRIADS as background

# I have from Luca a tsv with GO term information

#setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05")


# install GOseq
source("https://bioconductor.org/biocLite.R")
biocLite("goseq")

library(goseq)
library(reshape2)

######## first generate category mapping (GO information) ##############

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\GO_analysis")
transcripts_GO <- read.csv(file="transcripts_GO.csv")
head(transcripts_GO)
dim(transcripts_GO)
colnames(transcripts_GO)

# also need list of genes which are found in bins with 3 or more genes
gene_info <- read.csv("Y:\\expression_browser\\TGAC_assembly\\analysis\\ordered_genes\\all_genes_in_bins_with_3_genes.csv",header=T)
head(gene_info)
dim(gene_info)
gene_info <- gene_info[,2:3]
head(gene_info)

#already have 1 gene per line


# just select transcripts_GO which are included in the the gene_info

gene_info_transcripts_GO <- merge(gene_info,transcripts_GO, by.x ="Row.names", by.y = "transcript_id" )
head(gene_info_transcripts_GO)
dim(gene_info_transcripts_GO)

str(gene_info_transcripts_GO)

# need to convert factors to characters to do melt
i <- sapply(gene_info_transcripts_GO, is.factor)
gene_info_transcripts_GO[i] <- lapply(gene_info_transcripts_GO[i], as.character)
head(gene_info_transcripts_GO)
dim(gene_info_transcripts_GO)
str(gene_info_transcripts_GO)

# dont need the "chr_cM" column so remove it
colnames(gene_info_transcripts_GO)
gene_info_transcripts_GO <- cbind(gene_info_transcripts_GO$Row.names,gene_info_transcripts_GO[,3:63])
head(gene_info_transcripts_GO)


#melt to have 1 gene per line with 1 GO term
library(reshape2)
melted_gene_info_transcripts_GO <- melt(gene_info_transcripts_GO, id=c("gene_info_transcripts_GO$Row.names" ))
head(melted_gene_info_transcripts_GO)

colnames(melted_gene_info_transcripts_GO) <- c("transcript","GO_group","GO_term")
head(melted_gene_info_transcripts_GO,20)
dim(melted_gene_info_transcripts_GO)

# select only rows with a GO_term (it kept GO1 to GO61 for all transcripts even though most don't have GO10 +)

melted_gene_info_complete <- melted_gene_info_transcripts_GO[ grep("GO",melted_gene_info_transcripts_GO$GO_term) ,]
head(melted_gene_info_complete)
dim(melted_gene_info_complete)

# make dataframe to use in GOseq

GO_data.frame <- as.data.frame(cbind(as.character(melted_gene_info_complete$transcript),melted_gene_info_complete$GO_term))
colnames(GO_data.frame) <- c("id", "GO")
head(GO_data.frame)
dim(GO_data.frame)
str(GO_data.frame)
is.data.frame(GO_data.frame)

############## second generate list of gene lengths #######################

# already have gene length info in .fai

# read in gene length info

all_transcript_lengths <- read.table (file="Y:\\TGACv1_annotation_CS42_ensembl_release\\Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.cdna.fa.fai")
head(all_transcript_lengths)

# to select only transcripts which have GO terms and are in bins with 3 or more genes
# need to get list of unique transcripts with GO terms in bins with 3 or more genes

transcripts_with_GO_in_triads_list <- as.data.frame(unique(GO_data.frame$id))
head(transcripts_with_GO_in_triads_list)
dim(transcripts_with_GO_in_triads_list)

# merge this unique list with the gene length
 transcripts_with_GO_in_triads_lengths <- merge(transcripts_with_GO_in_triads_list, all_transcript_lengths[,1:2], by.x ="unique(GO_data.frame$id)", by.y ="V1")
 head(transcripts_with_GO_in_triads_lengths)
 dim(transcripts_with_GO_in_triads_lengths)
 
 colnames(transcripts_with_GO_in_triads_lengths) <- c("transcript","length")

 # need numeric vector of lengths
 
# how to from https://biosupport.se/p/126/
length.vector.names <- as.vector(transcripts_with_GO_in_triads_lengths$length)
names(length.vector.names) <- transcripts_with_GO_in_triads_lengths$transcript
head(length.vector.names)


######## third generate named vector of whether each gene is differentially expressed ############

# read in list of genes in "hotspots"
# melt the list so it one column of gene IDs

# merge the transcripts_with_GO_in_triads_lengths with the melted_conserved_genes 
####(make sure to keep all transcripts in the transcripts_with_GO_in_triads_lengths)
# convert it to a named vector http://stackoverflow.com/questions/19265172/converting-two-columns-of-a-data-frame-to-a-named-vector 

#move to folder where list of gene "hotspots" are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\ordered_genes\\")

hotspot_genes <- read.csv("hotspot_genes_over_10tpm.csv", header=T)
head(hotspot_genes,20)

list_of_regions <- c("TGAC_leaf", "TGAC_root", "TGAC_seed", "TGAC_seedling", "TGAC_spike", "TGAC_stem")


# get library for melting
library(reshape2)

# get library for goseq
library(goseq)
i=4
for (i in 1:6) {

# select just information about genes in cM of interest
gene_table <- hotspot_genes[ grep(list_of_regions[i],hotspot_genes$sample_with_median_tpm_over_10_in_bin),]

# or don't select just that region  
#gene_table <- hotspot_genes # use this line if want to use all genes in hotspots
# i = "all_genes_in_hotspots" # use this line if want to use all genes in hotspots

head(gene_table)
dim(gene_table)
str(gene_table)

# make vector of genes in the cM
genes.vector <- as.integer(transcripts_with_GO_in_triads_lengths$transcript%in%gene_table$transcript)
names(genes.vector) <- transcripts_with_GO_in_triads_lengths$transcript
head(genes.vector)
sesstail(head(genes.vector,5000) )
table(genes.vector)

# now run GO seq itself
# You can provide user-specified gene length information to nullp function of goseq Bioconductor package:
# gene_pwf=nullp(gene.vector.up, bias.data=length.vector.names)
pwf <- nullp(genes.vector,bias.data=length.vector.names)

head(pwf)

# calculate GO enrichment using default method
GO.WALL <- goseq(pwf, gene2cat=GO_data.frame)
head(GO.WALL)
dim(GO.WALL)

# need to make a p-value cut-off taking into account mutiple hypothesis testing correction
# e.g. GO categories over enriched using a 0.05 FDR cut-off (Benjamin and Hochberg 1995)

over_rep_GO <- GO.WALL[p.adjust(GO.WALL$over_represented_pvalue,method="BH")<0.05,]
head(over_rep_GO)
dim(over_rep_GO)

under_rep_GO <- GO.WALL[p.adjust(GO.WALL$under_represented_pvalue,method="BH")<0.05,]
head(under_rep_GO)
dim(under_rep_GO)


# write an output file for each conserved and divergent GO for each condition
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\ordered_genes\\")

write.table(over_rep_GO,file=paste("over_rep_GO_",list_of_regions[i],".txt"))
#write.table(over_rep_GO,file=paste("over_rep_GO_",i,".txt")) # use this line if want to use all genes in hotspots
write.table(under_rep_GO,file=paste("under_rep_GO_",list_of_regions[i],".txt"))
#write.table(under_rep_GO,file=paste("under_rep_GO_",i,".txt")) # use this line if want to use all genes in hotspots


}
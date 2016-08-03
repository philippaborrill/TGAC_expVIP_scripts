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

######## first generate category mapping (GO information) ##############

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\GO_analysis")
transcripts_GO <- read.csv(file="transcripts_GO.csv")
head(transcripts_GO)
dim(transcripts_GO)
colnames(transcripts_GO)

# also need triad info
triad_info <- read.table("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05\\triads_over_90_percent.txt",sep="\t", header=T)
head(triad_info)
dim(triad_info)

#melt triad info to have 1 gene per line
library(reshape2)
melted_triad_info <- melt(triad_info, id=c("Group"))
head(melted_triad_info)
dim(melted_triad_info)

# just select transcripts_GO which are included in the triads

triads_transcripts_GO <- merge(melted_triad_info,transcripts_GO, by.x ="value", by.y = "transcript_id" )
head(triads_transcripts_GO)
dim(triads_transcripts_GO)

str(triads_transcripts_GO)

# need to convert factors to characters to do melt
i <- sapply(triads_transcripts_GO, is.factor)
triads_transcripts_GO[i] <- lapply(triads_transcripts_GO[i], as.character)
head(triads_transcripts_GO)
dim(triads_transcripts_GO)
str(triads_transcripts_GO)

#melt to have 1 gene per line with 1 GO term
library(reshape2)
melted_triads_transcripts_GO <- melt(triads_transcripts_GO, id=c("value", "Group", "variable" ))
head(melted_triads_transcripts_GO)

colnames(melted_triads_transcripts_GO) <- c("transcript","Group", "gene", "GO_group","GO_term")
head(melted_triads_transcripts_GO)
dim(melted_triads_transcripts_GO)

# select only rows with a GO_term (it kept GO1 to GO61 for all transcripts even though most don't have GO10 +)

melted_triads_complete <- melted_triads_transcripts_GO[ grep("GO",melted_triads_transcripts_GO$GO_term) ,]
head(melted_triads_complete)
dim(melted_triads_complete)

# make dataframe to use in GOseq

GO_data.frame <- as.data.frame(cbind(melted_triads_complete$transcript,melted_triads_complete$GO_term))
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

# to select only transcripts which have GO terms and are in triads
# need to get list of unique transcripts with GO terms in triads

transcripts_with_GO_in_triads_list <- as.data.frame(unique(GO_data.frame$transcript))
head(transcripts_with_GO_in_triads_list)
dim(transcripts_with_GO_in_triads_list)

# merge this unique list with the gene length
 transcripts_with_GO_in_triads_lengths <- merge(transcripts_with_GO_in_triads_list, all_transcript_lengths[,1:2], by.x ="unique(GO_data.frame$transcript)", by.y ="V1")
 head(transcripts_with_GO_in_triads_lengths)
 dim(transcripts_with_GO_in_triads_lengths)
 
 colnames(transcripts_with_GO_in_triads_lengths) <- c("transcript","length")

 # need numeric vector of lengths
 
# how to from https://biosupport.se/p/126/
length.vector.names <- as.vector(transcripts_with_GO_in_triads_lengths$length)
names(length.vector.names) <- transcripts_with_GO_in_triads_lengths$transcript
head(length.vector.names)


######## third generate named vector of whether each gene is differentially expressed ############

# read in list of e.g. conserved genes in 1 h drought
# melt the list so it one column of gene IDs
# merge the transcripts_with_GO_in_triads_lengths with the melted_conserved_genes 
####(make sure to keep all transcripts in the transcripts_with_GO_in_triads_lengths)
# convert it to a named vector http://stackoverflow.com/questions/19265172/converting-two-columns-of-a-data-frame-to-a-named-vector 

#move to folder where conserved/divergent gene results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05\\")
list_of_files_cons_div <- list.files(pattern="*expr.txt")
list_of_files_cons_div
length(list_of_files_cons_div)

#remove the yellow rust 48 hour diverged because there are no genes in it
list_of_files_cons_div <- c(list_of_files_cons_div[1:21],list_of_files_cons_div[23:24])
list_of_files_cons_div

# get library for melting
library(reshape2)

# get library for goseq
library(goseq)


for (i in 1:23) {

 
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05\\")
  
gene_table <- read.table(file=list_of_files_cons_div[i], header=T)
 # read.table(file=list_of_files_cons_div[i])
head(gene_table)
dim(gene_table)
str(gene_table)
# need to convert factors to characters to do melt
j <- sapply(gene_table, is.factor)
gene_table[j] <- lapply(gene_table[j], as.character)
head(gene_table)
dim(gene_table)
str(gene_table)

melted_gene_table <- melt(gene_table, id=c("Group", "A", "B", "D"))
head(melted_gene_table)
head(transcripts_with_GO_in_triads_lengths)

genes.vector <- as.integer(transcripts_with_GO_in_triads_lengths$transcript%in%melted_gene_table$value)
names(genes.vector) <- transcripts_with_GO_in_triads_lengths$transcript
head(genes.vector)
tail(head(genes.vector,2000) )
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
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\GO_analysis\\")

write.table(over_rep_GO,file=paste("over_rep_GO_",list_of_files_cons_div[i]))

write.table(under_rep_GO,file=paste("under_rep_GO_",list_of_files_cons_div[i]))


}
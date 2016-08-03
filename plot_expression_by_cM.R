# Plot graphs of gene expression according to cM position for TGAC genome paper
# Philippa Borrill
# 12.07.2016

# Want to plot expression in 6 TGAC RNA samples across the chromosomes 
# to see if there is positional variation in expression levels

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

# only need the TGAC samples
tpm_TGAC <- tpm.data[,420:425]
head(tpm_TGAC)
head(rownames(tpm_TGAC))

# remove all tpm data which I don't need from the workspace
rm(tpm.data)
rm(tpm)


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

# merge together the expression tpm with the cM position info
tpm_cm <- merge(tpm_TGAC, cm, by=0) # by=0 means by row.names
head(tpm_cm)

# want to plot a separate graph for each chromosome (21) for each TGAC sample (6)
# will use a for loop to get a separate plot for each TGAC sample
# will use facet_grid in ggplot to get a separate graph for each chromosome

# first need to figure out how the final graph should look
library(ggplot2)

# check factors are correct for tpm_cm
str(tpm_cm)

# this makes a graph of the raw tpms across the different chromosomes
# quite messy

jpeg(file="TGAC_leaf_cM.jpg", height=480, width=480, units="px")
ggplot(tpm_cm, aes(cM,TGAC_leaf)) + geom_line() + facet_wrap(~ chromosome,ncol=3)
dev.off()

# try calculating median gene expression per genetic bin and plotting 

# first make a new column with the chrom_cM to use as a factor
tpm_cm_uniq_pos <- (transform(tpm_cm,chr_cM=paste0(chromosome,"_",cM)))

head(tpm_cm_uniq_pos)

# want to exclude genes with expression <2 tpm across all 6 samples
tpm_cm_uniq_pos <- tpm_cm_uniq_pos[ which(tpm_cm_uniq_pos$TGAC_leaf >2 
                         | tpm_cm_uniq_pos$TGAC_root > 2 | tpm_cm_uniq_pos$TGAC_seed > 2 
                         | tpm_cm_uniq_pos$TGAC_seedling > 2 | tpm_cm_uniq_pos$TGAC_spike > 2
                         | tpm_cm_uniq_pos$TGAC_stem > 2), ]
head(tpm_cm_uniq_pos)

### CAN NOW TRY TO PLOT EXPRESSION BREADTH - SEE SECTION AT END #####

# try using this filtered data to plot original graph (not averaged per bin)

jpeg(file="TGAC_leaf_cM_excl_under_2_tpm_genes.jpg", height=480, width=480, units="px")
ggplot(tpm_cm_uniq_pos, aes(cM,TGAC_leaf)) + geom_line() + facet_wrap(~ chromosome,ncol=3)
dev.off()

# now try averaging per cM bin
agg_tpm_cm_median <- aggregate(tpm_cm_uniq_pos[,2:7], by = list(tpm_cm_uniq_pos$chr_cM), FUN=median)
head(agg_tpm_cm_median)
dim(agg_tpm_cm_median)

agg_tpm_cm_mean <- aggregate(tpm_cm_uniq_pos[,2:7], by = list(tpm_cm_uniq_pos$chr_cM), FUN=mean)
head(agg_tpm_cm_mean)
dim(agg_tpm_cm_mean)

# for interest lets count the number of genes in each bin
agg_tpm_cm_count <- aggregate(tpm_cm_uniq_pos[,2:7], by = list(tpm_cm_uniq_pos$chr_cM), FUN= length)
head(agg_tpm_cm_count)
dim(agg_tpm_cm_count)

##### AT THIS POINT CAN GO TO THE END OF THE DOCUMENT TO SEE HOW TO PLOT ONLY INCLUDING BINS WITH 3 OR MORE GENES ###############

# plot the distribution of number of genes per genetic bin
hist(agg_tpm_cm_count$TGAC_leaf, breaks=100, xlab="number of genes per bin", main="Number of genes per bin TGAC_leaf")
hist(agg_tpm_cm_count$TGAC_leaf, breaks=1000, xlim=c(0,50), xlab="number of genes per bin", main="Number of genes per bin TGAC_leaf")

# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_tpm_cm_median$Group.1),'_') 
head(split_cols)
agg_tpm_cm_median <- data.frame(agg_tpm_cm_median[,2:7], do.call(rbind, split_cols))
colnames(agg_tpm_cm_median) <- c(colnames(agg_tpm_cm_median[1:6]), "chromosome", "cM")
head(agg_tpm_cm_median)
dim(agg_tpm_cm_median)
agg_tpm_cm_median <-as.data.frame (agg_tpm_cm_median)

# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_tpm_cm_mean$Group.1),'_') 
head(split_cols)
agg_tpm_cm_mean <- data.frame(agg_tpm_cm_mean[,2:7], do.call(rbind, split_cols))
colnames(agg_tpm_cm_mean) <- c(colnames(agg_tpm_cm_mean[1:6]), "chromosome", "cM")
head(agg_tpm_cm_mean)
dim(agg_tpm_cm_mean)
agg_tpm_cm_mean <-as.data.frame (agg_tpm_cm_mean)




# check factors are correct in new dataframes
str(agg_tpm_cm_mean)
str(agg_tpm_cm_median)

# now plot median or mean of tpm per cM bin
jpeg(file="TGAC_leaf_cM_excl_under_2_tpm_genes_median_per_bin.jpg", height=480, width=700, units="px")
ggplot(agg_tpm_cm_median, aes(cM,TGAC_leaf,group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3, scales="free_y") # have to add "group=1" to make it plot the lines because the data series needs to belong to a group
dev.off()

jpeg(file="TGAC_leaf_cM_excl_under_2_tpm_genes_mean_per_bin.jpg", height=480, width=700, units="px")
ggplot(agg_tpm_cm_mean, aes(x = cM,y = TGAC_leaf, group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3, scales="free_y")  
dev.off()

# without free scale
jpeg(file="TGAC_leaf_cM_excl_under_2_tpm_genes_median_per_bin_uniform.jpg", height=480, width=700, units="px")
ggplot(agg_tpm_cm_median, aes(cM,TGAC_leaf,group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3) # have to add "group=1" to make it plot the lines because the data series needs to belong to a group
dev.off()

jpeg(file="TGAC_leaf_cM_excl_under_2_tpm_genes_mean_per_bin_uniform.jpg", height=480, width=700, units="px")
ggplot(agg_tpm_cm_mean, aes(x = cM,y = TGAC_leaf, group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3)  
dev.off()

# try plotting with data smoothing
# Smoothed symmetrically:
# average of current sample, 2 future samples, and 2 past samples 
f5 <- rep(1/5,5)
f5

smoothed_y <- filter(agg_tpm_cm_mean$TGAC_leaf, f5, sides=2)

jpeg(file="test_TGAC_leaf_cM_excl_under_2_tpm_genes_mean_per_bin_uniform.jpg", height=480, width=700, units="px")
ggplot(agg_tpm_cm_mean, aes(x = cM,y = TGAC_leaf, group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3)  
ggplot(agg_tpm_cm_mean, aes(x = cM,y = smoothed_y, group=1)) + geom_line(colour ="red") + facet_wrap(~ chromosome,ncol=3)  
dev.off()

# I think the smoothing just makes it confusing and looks like large regions are "high expression" 
# when in actual fact it's just one bin with high expression


# instead of plotting one by one let's plot using a loop
samples <- colnames(agg_tpm_cm_median)[1:6]
samples

for (i in 1:6) {
#  print(i)
#  print (samples[i])
  current_sample <- samples[[i]]
  print(current_sample)
  
  # now plot median or mean of tpm per cM bin
  jpeg(file=paste(current_sample,"_cM_excl_under_2_tpm_genes_median_per_bin.jpg", sep=""), height=480, width=700, units="px")
  print(ggplot(agg_tpm_cm_median, aes_string(x="cM",y=current_sample,group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3, scales="free_y") )
  # have to add "group=1" to make it plot the lines because the data series needs to belong to a group
  # have to use aes_string or else ggplot won't recognise current_sample as a variable
  dev.off()
  
  jpeg(file=paste(current_sample,"_cM_excl_under_2_tpm_genes_mean_per_bin.jpg", sep=""), height=480, width=700, units="px")
  print(ggplot(agg_tpm_cm_mean, aes_string(x="cM", y=current_sample, group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3, scales="free_y") )
  dev.off()
  
  # do with uniform scales
  jpeg(file=paste(current_sample,"_cM_excl_under_2_tpm_genes_median_per_bin_uniform.jpg", sep=""), height=480, width=700, units="px")
  print(ggplot(agg_tpm_cm_median, aes_string(x="cM",y=current_sample,group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3) )
  # have to add "group=1" to make it plot the lines because the data series needs to belong to a group
  # have to use aes_string or else ggplot won't recognise current_sample as a variable
  dev.off()
  
  jpeg(file=paste(current_sample,"_cM_excl_under_2_tpm_genes_mean_per_bin_uniform.jpg", sep=""), height=480, width=700, units="px")
  print(ggplot(agg_tpm_cm_mean, aes_string(x="cM", y=current_sample, group=1)) + geom_line() + facet_wrap(~ chromosome,ncol=3) )
  dev.off()
  
}

# now want to try "melting" the data so I can plot all 6 tissues on the same graphs
# data to use is the binned data

head(agg_tpm_cm_mean)
head(agg_tpm_cm_median)

library(reshape2)
library(RColorBrewer)

melted_agg_tpm_cm_mean <- melt(agg_tpm_cm_mean, id=c("chromosome", "cM"))
head(melted_agg_tpm_cm_mean)
tail(melted_agg_tpm_cm_mean)

# need to convert cM to "numeric" to get scaling to work on plot
melted_agg_tpm_cm_mean$cM=as.numeric(levels(melted_agg_tpm_cm_mean$cM))[melted_agg_tpm_cm_mean$cM]

# store colour-blind friendly pallete
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make as eps
postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_mean_per_bin.eps")
ggplot(melted_agg_tpm_cm_mean, aes(cM,value,group=variable,color=variable)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()


# to make graph with colour brewer colours
jpeg(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_mean_per_bin_rainbow.jpg", height=480, width=900, units="px")
ggplot(melted_agg_tpm_cm_mean, aes(cM,value,group=variable,color=variable)) + geom_line() + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_brewer(palette="Set1")
dev.off()

# do for median as well

melted_agg_tpm_cm_median <- melt(agg_tpm_cm_median, id=c("chromosome", "cM"))
head(melted_agg_tpm_cm_median)
tail(melted_agg_tpm_cm_median)

# need to convert cM to "numeric" to get scaling to work on plot
melted_agg_tpm_cm_median$cM=as.numeric(levels(melted_agg_tpm_cm_median$cM))[melted_agg_tpm_cm_median$cM]

# store colour-blind friendly pallete
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make as eps
postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_per_bin.eps")
ggplot(melted_agg_tpm_cm_median, aes(cM,value,group=variable,color=variable)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()

##### PLOTTING INCLUDING ONLY BINS WITH 3 OR MORE GENES ###############

# we have aggregated mean and median per bin
# now try averaging per cM bin
head(agg_tpm_cm_mean)
head(agg_tpm_cm_median,10)

# and we have the number of genes in each bin
head(agg_tpm_cm_count,10)

# the original data (not averaged per cm was)
head(tpm_cm_uniq_pos)

# we want to only keep cM bins which have 3 or more genes

agg_tpm_cm_median_filt <- agg_tpm_cm_median[agg_tpm_cm_count$TGAC_leaf>2,]
head(agg_tpm_cm_median_filt,10)

agg_tpm_cm_mean_filt <- agg_tpm_cm_mean[agg_tpm_cm_count$TGAC_leaf>2,]
head(agg_tpm_cm_mean_filt,10)

# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_tpm_cm_median_filt$Group.1),'_') 
head(split_cols)
agg_tpm_cm_median_filt <- data.frame(agg_tpm_cm_median_filt[,2:7], do.call(rbind, split_cols))
colnames(agg_tpm_cm_median_filt) <- c(colnames(agg_tpm_cm_median_filt[1:6]), "chromosome", "cM")
head(agg_tpm_cm_median_filt)
dim(agg_tpm_cm_median_filt)
agg_tpm_cm_median_filt <-as.data.frame (agg_tpm_cm_median_filt)

# need to split 1st column back into chrom and cM
split_cols <- strsplit(as.character(agg_tpm_cm_mean_filt$Group.1),'_') 
head(split_cols)
agg_tpm_cm_mean_filt <- data.frame(agg_tpm_cm_mean_filt[,2:7], do.call(rbind, split_cols))
colnames(agg_tpm_cm_mean_filt) <- c(colnames(agg_tpm_cm_mean_filt[1:6]), "chromosome", "cM")
head(agg_tpm_cm_mean_filt)
dim(agg_tpm_cm_mean_filt)
agg_tpm_cm_mean_filt <-as.data.frame (agg_tpm_cm_mean_filt)

# now try plotting the graphs using this "filt" data which doesn't contain bins with under 3 genes in them
install.packages("reshape2")
install.packages("stringi")
library(reshape2)
library(RColorBrewer)

melted_agg_tpm_cm_mean_filt <- melt(agg_tpm_cm_mean_filt, id=c("chromosome", "cM"))
head(melted_agg_tpm_cm_mean_filt)
tail(melted_agg_tpm_cm_mean_filt)

# need to convert cM to "numeric" to get scaling to work on plot
melted_agg_tpm_cm_mean_filt$cM=as.numeric(levels(melted_agg_tpm_cm_mean_filt$cM))[melted_agg_tpm_cm_mean_filt$cM]

# store colour-blind friendly pallete
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make as eps
postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_mean_per_bin_filt_3_genes_per_bin.eps")
ggplot(melted_agg_tpm_cm_mean_filt, aes(cM,value,group=variable,color=variable)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()

postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_mean_per_bin_filt_3_genes_per_bin_points.eps")
ggplot(melted_agg_tpm_cm_mean_filt, aes(cM,value,group=variable,color=variable)) + geom_point(size=1.2) + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()

# do for median as well

melted_agg_tpm_cm_median_filt <- melt(agg_tpm_cm_median_filt, id=c("chromosome", "cM"))
head(melted_agg_tpm_cm_median_filt)
tail(melted_agg_tpm_cm_median_filt)
dim(melted_agg_tpm_cm_median_filt)

#### after I plotted the graph (see below) I realised I wanted to look at bins which are expressed at over 20 tpm ####
# want to know what genes are in them- is there something special about these "hotspots"

bins_over_20_tpm_median <- melted_agg_tpm_cm_median_filt[melted_agg_tpm_cm_median_filt$value>20 , ]
dim(bins_over_20_tpm_median)
(bins_over_20_tpm_median)
sorted_bins_over_20_tpm_median <- bins_over_20_tpm_median[ order(bins_over_20_tpm_median[,1], bins_over_20_tpm_median[,2]), ]
dim(sorted_bins_over_20_tpm_median)
sorted_bins_over_20_tpm_median
# write list of cM bin with median expression > 20 tpm 
write.csv(sorted_bins_over_20_tpm_median, file="bins_median_over_20_tpm_TGAC_tissues.csv")

# how many genes are in these bins
# need to give each bin a unique ID
sorted_bins_over_20_tpm_median_uniq_pos <- (transform(sorted_bins_over_20_tpm_median,chr_cM=paste0(chromosome,"_",cM)))
head(sorted_bins_over_20_tpm_median_uniq_pos)
head(agg_tpm_cm_count,20)

# merge the list of bins with >20tpm with the count list of genes per bin
number_of_genes_per_20tpm_median_bin <- merge(sorted_bins_over_20_tpm_median_uniq_pos,agg_tpm_cm_count, by.x="chr_cM", by.y = "Group.1")
head(number_of_genes_per_20tpm_median_bin)

# get rid of separate column for each tissue - only need for 1 tissue because all tissues have same number of genes per bin
number_of_genes_per_20tpm_median_bin <- number_of_genes_per_20tpm_median_bin[,1:6]
colnames(number_of_genes_per_20tpm_median_bin) <- c(colnames(number_of_genes_per_20tpm_median_bin)[1:3], "tissue", "median_tpm", "genes_in_bin")
number_of_genes_per_20tpm_median_bin
hist(number_of_genes_per_20tpm_median_bin$genes_in_bin)
hist(agg_tpm_cm_count$TGAC_leaf, xlim = c(0,100), breaks=1000)

# are the number of genes per bin statistically significantly different in "high median tpm" bins from all bins
t.test(number_of_genes_per_20tpm_median_bin$genes_in_bin,agg_tpm_cm_count$TGAC_leaf)

# what genes are in these bins
# need to use my list of bins with >20tpm to select from list of genes which are in those bins

# have list of bins with >20tpm median
head(sorted_bins_over_20_tpm_median_uniq_pos)

# have list of genes with bin location
head(tpm_cm_uniq_pos)

# now select only genes within the "hotspot" bins
hotspot_genes <- merge(tpm_cm_uniq_pos,sorted_bins_over_20_tpm_median_uniq_pos,by.x= "chr_cM", by.y ="chr_cM")
dim(tpm_cm_uniq_pos)
dim(hotspot_genes)
head(hotspot_genes)
colnames(hotspot_genes) <- c("chr_cM","transcript", colnames(hotspot_genes)[3:16],"sample_with_median_tpm_over_20_in_bin", "median_tpm_in_bin" )
head(hotspot_genes)

# write csv with the "hotspot" genes
write.csv(hotspot_genes,file="hotspot_genes.csv")

#### after I plotted the graph (see below) I realised I wanted to look at bins which are expressed at over 10 tpm ####
# want to know what genes are in them- is there something special about these "hotspots"

bins_over_10_tpm_median <- melted_agg_tpm_cm_median_filt[melted_agg_tpm_cm_median_filt$value>10 , ]
dim(bins_over_10_tpm_median)
head(bins_over_10_tpm_median)
sorted_bins_over_10_tpm_median <- bins_over_10_tpm_median[ order(bins_over_10_tpm_median[,1], bins_over_10_tpm_median[,2]), ]
dim(sorted_bins_over_10_tpm_median)
head(sorted_bins_over_10_tpm_median)
# write list of cM bin with median expression > 20 tpm 
write.csv(sorted_bins_over_10_tpm_median, file="bins_median_over_10_tpm_TGAC_tissues.csv")

# how many genes are in these bins
# need to give each bin a unique ID
sorted_bins_over_10_tpm_median_uniq_pos <- (transform(sorted_bins_over_10_tpm_median,chr_cM=paste0(chromosome,"_",cM)))
head(sorted_bins_over_10_tpm_median_uniq_pos)
head(agg_tpm_cm_count,20)

# merge the list of bins with >10tpm with the count list of genes per bin
number_of_genes_per_10tpm_median_bin <- merge(sorted_bins_over_10_tpm_median_uniq_pos,agg_tpm_cm_count, by.x="chr_cM", by.y = "Group.1")
head(number_of_genes_per_10tpm_median_bin)

# get rid of separate column for each tissue - only need for 1 tissue because all tissues have same number of genes per bin
number_of_genes_per_10tpm_median_bin <- number_of_genes_per_10tpm_median_bin[,1:6]
colnames(number_of_genes_per_10tpm_median_bin) <- c(colnames(number_of_genes_per_10tpm_median_bin)[1:3], "tissue", "median_tpm", "genes_in_bin")
number_of_genes_per_10tpm_median_bin
hist(number_of_genes_per_10tpm_median_bin$genes_in_bin)
hist(agg_tpm_cm_count$TGAC_leaf, xlim = c(0,100), breaks=1000)

# calc mean transcripts per bin in 
mean(number_of_genes_per_10tpm_median_bin$genes_in_bin)
median(number_of_genes_per_10tpm_median_bin$genes_in_bin)

# calc mean transcripts per bin all bins
mean(agg_tpm_cm_count$TGAC_leaf)
median(agg_tpm_cm_count$TGAC_leaf)

# are the number of genes per bin statistically significantly different in "high median tpm" bins from all bins
t.test(number_of_genes_per_10tpm_median_bin$genes_in_bin,agg_tpm_cm_count$TGAC_leaf)

# what genes are in these bins
# need to use my list of bins with >10tpm to select from list of genes which are in those bins

# have list of bins with >10tpm median
head(sorted_bins_over_10_tpm_median_uniq_pos)

# have list of genes with bin location
head(tpm_cm_uniq_pos)

# now select only genes within the "hotspot" bins
hotspot_genes <- merge(tpm_cm_uniq_pos,sorted_bins_over_10_tpm_median_uniq_pos,by.x= "chr_cM", by.y ="chr_cM")
dim(tpm_cm_uniq_pos)
dim(hotspot_genes)
head(hotspot_genes)
colnames(hotspot_genes) <- c("chr_cM","transcript", colnames(hotspot_genes)[3:16],"sample_with_median_tpm_over_10_in_bin", "median_tpm_in_bin" )
head(hotspot_genes)

# write csv with the "hotspot" genes
write.csv(hotspot_genes,file="hotspot_genes_over_10tpm.csv")


# need to also get the list of genes which are found within bins with 3 or more genes

# already have list of bins iwth 3 or more genes:
head(agg_tpm_cm_median_filt,10)
all_genes_in_bins_with_3_genes <- merge(tpm_cm_uniq_pos,agg_tpm_cm_median_filt,by.x= "chr_cM", by.y ="Group.1")
head(all_genes_in_bins_with_3_genes)
dim(all_genes_in_bins_with_3_genes)

# write csv with all genes in bins with 3 genes:
write.csv(all_genes_in_bins_with_3_genes, file="all_genes_in_bins_with_3_genes.csv")

# need to convert cM to "numeric" to get scaling to work on plot
melted_agg_tpm_cm_median_filt$cM=as.numeric(levels(melted_agg_tpm_cm_median_filt$cM))[melted_agg_tpm_cm_median_filt$cM]

# store colour-blind friendly pallete
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make as eps
postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_per_bin_filt_3_genes_per_bin.eps")
ggplot(melted_agg_tpm_cm_median_filt, aes(cM,value,group=variable,color=variable)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()

postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_per_bin_filt_3_genes_per_bin_points.eps")
ggplot(melted_agg_tpm_cm_median_filt, aes(cM,value,group=variable,color=variable)) + geom_point(size=1.2) + facet_wrap(~ chromosome,ncol=3, scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()

postscript(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_per_bin_filt_3_genes_per_bin_points_uniform.eps")
ggplot(melted_agg_tpm_cm_median_filt, aes(cM,value,group=variable,color=variable)) + geom_point(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + scale_y_continuous(breaks=c(0,10,20,30), limits=c(0,35)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + theme(axis.text.y=element_text(size=10))
dev.off()

# want to add in centromere position

# read in file which has centromere approx position
centromere_pos<-read.table("get_centromere_position.txt",sep="\t")  
head(centromere_pos)
dim(centromere_pos)

colnames(centromere_pos) <- c("chromosome","centromere_pos")
head(centromere_pos)

# try plotting graph
pdf(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_per_bin_filt_3_genes_per_bin_points_with_centromere.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(melted_agg_tpm_cm_median_filt, aes(cM,value,group=variable,color=variable)) + geom_point(size=1.2) + facet_wrap(~ chromosome,ncol=3,scales="free_y") + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + scale_y_continuous(breaks=c(0,10,20,30), limits=c(0,35)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()

# if wanted to plot geom_vline underneath the ggplot should put that line of code first

# with uniform scale
pdf(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_per_bin_filt_3_genes_per_bin_points_with_centromere_uniform.pdf", height=200, width=200)
ggplot(melted_agg_tpm_cm_median_filt, aes(cM,value,group=variable,color=variable)) + geom_point(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks=c(0,50,100, 150,200)) + scale_y_continuous(breaks=c(0,10,20,30), limits=c(0,35)) + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()





### CAN NOW TRY TO PLOT EXPRESSION BREADTH - SEE SECTION AT END #####
head(tpm_cm_uniq_pos)

# need to put a true or false in each tissue whether a gene is expressed to over 2 tpm
head(tpm_cm_uniq_pos)


#get a logical vector for each tissue and join them together in a new dataframe

tpm_over_2_logical <- cbind.data.frame(tpm_cm_uniq_pos$Row.names, tpm_cm_uniq_pos$TGAC_leaf >2, tpm_cm_uniq_pos$TGAC_root >2, tpm_cm_uniq_pos$TGAC_seed >2, 
                                       tpm_cm_uniq_pos$TGAC_seedling >2, tpm_cm_uniq_pos$TGAC_spike >2, tpm_cm_uniq_pos$TGAC_stem >2, tpm_cm_uniq_pos$scaffold,
                                       tpm_cm_uniq_pos$gene, tpm_cm_uniq_pos$biotype, tpm_cm_uniq_pos$confidence, tpm_cm_uniq_pos$chromosome, 
                                       tpm_cm_uniq_pos$cM, tpm_cm_uniq_pos$chr_cM )
head(tpm_over_2_logical)
 dim(tpm_over_2_logical)               
# add up the number of trues in for each gene (i.e. number of tissues expressed in)
   
true_per_gene <- rowSums(tpm_over_2_logical[,2:7])
length(true_per_gene)

tpm_over_2_logical_summary <- cbind.data.frame(tpm_cm_uniq_pos$Row.names, true_per_gene, tpm_cm_uniq_pos$scaffold,
                                       tpm_cm_uniq_pos$gene, tpm_cm_uniq_pos$biotype, tpm_cm_uniq_pos$confidence, tpm_cm_uniq_pos$chromosome, 
                                       tpm_cm_uniq_pos$cM, tpm_cm_uniq_pos$chr_cM )
head(tpm_over_2_logical_summary)
colnames(tpm_over_2_logical_summary)
#assign correct colnames
colnames(tpm_over_2_logical_summary) <- c("Row.names", "conditions_per_gene", "scaffold", "gene" , "biotype","confidence","chromosome", "cM", "chr_cM")
colnames(tpm_over_2_logical_summary)
head(tpm_over_2_logical_summary)

# now try averaging per cM bin
agg_conditions_cm_median <- aggregate(tpm_over_2_logical_summary[,2], by = list(tpm_over_2_logical_summary$chr_cM), FUN=median)
head(agg_conditions_cm_median)
dim(agg_conditions_cm_median)

agg_conditions_cm_mean <- aggregate(tpm_over_2_logical_summary[,2], by = list(tpm_over_2_logical_summary$chr_cM), FUN=mean)
head(agg_conditions_cm_mean)
dim(agg_conditions_cm_mean)


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



# try plotting graph
pdf(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_median_expression_breadth_per_bin_with_centromere.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(agg_conditions_cm_median, aes(cM,conditions_per_gene)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values="black") +
  scale_x_continuous(breaks=c(0,50,100, 150,200))  + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()


# try plotting graph
pdf(file="TGAC_all_tissues_cM_excl_under_2_tpm_genes_mean_expression_breadth_per_bin_with_centromere.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(agg_conditions_cm_mean, aes(cM,conditions_per_gene)) + geom_line(size=1.2) + facet_wrap(~ chromosome,ncol=3) + scale_colour_manual(values="black") +
  scale_x_continuous(breaks=c(0,50,100, 150,200))  + theme(axis.text.x=element_text(size=15)) +theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.y=element_text(size=10)) + geom_vline(data=centromere_pos,aes(xintercept=as.numeric(centromere_pos)),size=1,colour="grey",alpha=0.5)
dev.off()
# Comparison of promoters from triads for seedling stress experiment for TGAC genome paper
# Philippa Borrill
# 14.07.2016

## README ##
# uses the output of the scripts:
# compare_promoter_similarity_triads.pl

########### for mRNA 500 bp conserved vs diverged ##################

## get promoter perc ID for all triads

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity")

triads_with_PID <- read.csv(file = "list_of_triads_with_perc_id.csv", header=F)
head(triads_with_PID)
dim(triads_with_PID)

# substitute away the "mRNA.500.bp" from the ID column, add ID header, also calculate Means column with mean per row

triad_mean_PID <- data.frame(ID=gsub(".mRNA.500.bp","",triads_with_PID[,1]),Means=rowMeans(triads_with_PID[,-1]))
head(triad_mean_PID)

hist(triad_mean_PID[,2], xlab="mean percentage ID within triad", main="Promoter (500 bp upstream of mRNA)
     percentage ID within triads", xlim=c(0,100), ylim=c(0,1200))


# calculate promoter perc ID for conserved vs diverged triads for each stress


#move to folder where conserved/divergent gene results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\")
list_of_files_cons_div <- list.files(pattern="*expr.txt")
list_of_files_cons_div
length(list_of_files_cons_div)

#remove the yellow rust 48 hour diverged because there are no genes in it
list_of_files_cons_div <- c(list_of_files_cons_div[1:21],list_of_files_cons_div[23:24])
list_of_files_cons_div

results <- c("pattern","mean_perc_id")

i=1

for (i in 1:23) {
  
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\")
  
  gene_table <- read.table(file=list_of_files_cons_div[i], header=T)
  # read.table(file=list_of_files_cons_div[i])
  head(gene_table)
  dim(gene_table)
  str(gene_table)
  
  merged_gene_table <- merge(gene_table,triad_mean_PID,by.x="Group", by.y="ID")
  head(merged_gene_table)
  dim(merged_gene_table)
  mean_perc_id <- mean(merged_gene_table$Means)
  
  results <- rbind(results,c(sub("_expr.txt","",list_of_files_cons_div[i]),mean_perc_id))
  results
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity\\")
  
  jpeg(file=paste((sub("_expr.txt","",list_of_files_cons_div[i])),".jpg"), height=480, width=480, units="px")
hist(merged_gene_table$Means, xlab="percentage ID within triad", main=paste("Promoter (500 bp upstream of mRNA)
     percentage ID within", (sub("_expr.txt","",list_of_files_cons_div[i])), "triads"), xlim=c(0,100))
  dev.off()
  
}  
results

write.table(results,file="mean_perc_id_sim_500bp_upstream_of_mRNA.txt")
  
########### for mRNA 500 bp opposite vs 3 UP vs 3 DOWN ##################

## get promoter perc ID for all triads

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity")

triads_with_PID <- read.csv(file = "list_of_triads_with_perc_id.csv", header=F)
head(triads_with_PID)
dim(triads_with_PID)

# substitute away the "mRNA.500.bp" from the ID column, add ID header, also calculate Means column with mean per row

triad_mean_PID <- data.frame(ID=gsub(".mRNA.500.bp","",triads_with_PID[,1]),Means=rowMeans(triads_with_PID[,-1]))
head(triad_mean_PID)
mean(triad_mean_PID$Means)
median(triad_mean_PID$Means)

hist(triad_mean_PID[,2], xlab="mean percentage ID within triad", main="Promoter (500 bp upstream of mRNA)
     percentage ID within triads", xlim=c(0,100), ylim=c(0,1200))


# calculate promoter perc ID for conserved vs diverged triads for each stress


#move to folder where 3UP, 3DOWN, opposite group results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\groups_up_down_opposite")
list_of_files_cons_div <- list.files(pattern="*homoeologues.txt")
list_of_files_cons_div
length(list_of_files_cons_div)

#remove the samples with no groups in them
list_of_files_cons_div <- c(list_of_files_cons_div[3:12],list_of_files_cons_div[14:18], 
                            list_of_files_cons_div[21],list_of_files_cons_div[24],list_of_files_cons_div[27],list_of_files_cons_div[30],list_of_files_cons_div[36])
list_of_files_cons_div

results <- c("pattern","mean_perc_id")

i=1

for (i in 1:20) {
  
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\groups_up_down_opposite")
  
  gene_table <- read.table(file=list_of_files_cons_div[i], header=T)
  # read.table(file=list_of_files_cons_div[i])
  head(gene_table)
  dim(gene_table)
  str(gene_table)
  
  merged_gene_table <- merge(gene_table,triad_mean_PID,by.x=0, by.y="ID")
  head(merged_gene_table)
  dim(merged_gene_table)
  mean_perc_id <- mean(merged_gene_table$Means)
  
  results <- rbind(results,c(sub("_homoeologues.txt","",list_of_files_cons_div[i]),mean_perc_id))
  results
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity\\")
  
  jpeg(file=paste((sub("_homoeologues.txt","",list_of_files_cons_div[i])),".jpg"), height=480, width=480, units="px")
  hist(merged_gene_table$Means, xlab="percentage ID within triad", main=paste("Promoter (500 bp upstream of mRNA)
     percentage ID within", (sub("_homoeologues.txt","",list_of_files_cons_div[i])), "triads"), xlim=c(0,100))
  dev.off()
  
}  
results

write.table(results,file="3up_3down_opp_mean_perc_id_sim_500bp_upstream_of_mRNA.txt")


########### for mRNA 500 bp opposite vs 3 UP vs 3 DOWN merged across stresses ##################

## get promoter perc ID for all triads

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity")

triads_with_PID <- read.csv(file = "list_of_triads_with_perc_id.csv", header=F)
head(triads_with_PID)
dim(triads_with_PID)

# substitute away the "mRNA.500.bp" from the ID column, add ID header, also calculate Means column with mean per row

triad_mean_PID <- data.frame(ID=gsub(".mRNA.500.bp","",triads_with_PID[,1]),Means=rowMeans(triads_with_PID[,-1]))
head(triad_mean_PID)
mean(triad_mean_PID$Means)
median(triad_mean_PID$Means)

hist(triad_mean_PID[,2], xlab="mean percentage ID within triad", main="Promoter (500 bp upstream of mRNA)
     percentage ID within triads", xlim=c(0,100), ylim=c(0,1200))


# calculate promoter perc ID for conserved vs diverged triads for each stress


#move to folder where 3UP, 3DOWN, opposite group results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\groups_up_down_opposite")

up_abiotic <- read.table("up_abiotic.txt")
head(up_abiotic)

down_abiotic <- read.table("down_abiotic.txt")
head(down_abiotic)

opposite_abiotic <- read.table("opposite_abiotic.txt")
head(opposite_abiotic)

# check if overlap between lists
num_triads_in_3UP_list_and_3DOWN_list <- merge(up_abiotic,down_abiotic,by="V1")
num_triads_in_3UP_list_and_3DOWN_list
num_triads_in_3UP_list_and_opposite_list <- merge(up_abiotic,opposite_abiotic,by="V1")
num_triads_in_3UP_list_and_opposite_list
num_triads_in_3DOWN_list_and_opposite_list <- merge(down_abiotic,opposite_abiotic,by="V1")
num_triads_in_3DOWN_list_and_opposite_list

# find out percentage ID of promoters within each list

merged_up_abiotic <- merge(up_abiotic,triad_mean_PID,by.x="V1", by.y="ID")
head(merged_up_abiotic)
dim(merged_up_abiotic)
mean_perc_id_up_abiotic <- mean(merged_up_abiotic$Means)
mean_perc_id_up_abiotic

median_perc_id_up_abiotic <- median(merged_up_abiotic$Means)
median_perc_id_up_abiotic

merged_down_abiotic <- merge(down_abiotic,triad_mean_PID,by.x="V1", by.y="ID")
head(merged_down_abiotic)
dim(merged_down_abiotic)
mean_perc_id_down_abiotic <- mean(merged_down_abiotic$Means)
mean_perc_id_down_abiotic

median_perc_id_down_abiotic <- median(merged_down_abiotic$Means)
median_perc_id_down_abiotic

merged_opposite_abiotic <- merge(opposite_abiotic,triad_mean_PID,by.x="V1", by.y="ID")
head(merged_opposite_abiotic)
dim(merged_opposite_abiotic)
mean_perc_id_opposite_abiotic <- mean(merged_opposite_abiotic$Means)
mean_perc_id_opposite_abiotic

median_perc_id_opposite_abiotic <- median(merged_opposite_abiotic$Means)
median_perc_id_opposite_abiotic

t.test(merged_up_abiotic$Means, merged_down_abiotic$Means, var.equal = T)

t.test(merged_up_abiotic$Means, merged_opposite_abiotic$Means, var.equal = T)

t.test(merged_opposite_abiotic$Means, merged_down_abiotic$Means, var.equal = T)



########### for mRNA 1000 bp conserved vs diverged ##################

## get promoter perc ID for all triads

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity")

triads_with_PID <- read.csv(file = "list_of_triads_with_perc_id_1000_mRNA.csv", header=F)
head(triads_with_PID)
dim(triads_with_PID)

# substitute away the "mRNA.1000.bp" from the ID column, add ID header, also calculate Means column with mean per row

triad_mean_PID <- data.frame(ID=gsub(".mRNA.1000.bp","",triads_with_PID[,1]),Means=rowMeans(triads_with_PID[,-1]))
head(triad_mean_PID)

jpeg(file="histogram_of_mean_perc_ID_within_triad_mRNA_1000bp_upstream.jpg", height=500, width =500)
hist(triad_mean_PID[,2], xlab="mean percentage ID within triad", main="Promoter (1000 bp upstream of mRNA)
     percentage ID within triads", xlim=c(0,100), ylim=c(0,1500))
dev.off()

mean(triad_mean_PID$Means)
median(triad_mean_PID$Means)

# calculate promoter perc ID for conserved vs diverged triads for each stress


#move to folder where conserved/divergent gene results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\")
list_of_files_cons_div <- list.files(pattern="*expr.txt")
list_of_files_cons_div
length(list_of_files_cons_div)

#remove the yellow rust 48 hour diverged because there are no genes in it
list_of_files_cons_div <- c(list_of_files_cons_div[1:21],list_of_files_cons_div[23:24])
list_of_files_cons_div

results <- c("pattern","mean_perc_id")

i=1

for (i in 1:23) {
  
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\")
  
  gene_table <- read.table(file=list_of_files_cons_div[i], header=T)
  # read.table(file=list_of_files_cons_div[i])
  head(gene_table)
  dim(gene_table)
  str(gene_table)
  
  merged_gene_table <- merge(gene_table,triad_mean_PID,by.x="Group", by.y="ID")
  head(merged_gene_table)
  dim(merged_gene_table)
  mean_perc_id <- mean(merged_gene_table$Means)
  
  results <- rbind(results,c(sub("_expr.txt","",list_of_files_cons_div[i]),mean_perc_id))
  results
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity\\")
  
  jpeg(file=paste((sub("_expr.txt","",list_of_files_cons_div[i])),".1000bp.jpg"), height=480, width=480, units="px")
  hist(merged_gene_table$Means, xlab="percentage ID within triad", main=paste("Promoter (1000 bp upstream of mRNA)
                                                                              percentage ID within", (sub("_expr.txt","",list_of_files_cons_div[i])), "triads"), xlim=c(0,100))
  dev.off()
  
}  
results

write.table(results,file="mean_perc_id_sim_1000bp_upstream_of_mRNA.txt")

########### for mRNA 1000 bp opposite vs 3 UP vs 3 DOWN ##################

## get promoter perc ID for all triads

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity")

triads_with_PID <- read.csv(file = "list_of_triads_with_perc_id_1000_mRNA.csv", header=F)
head(triads_with_PID)
dim(triads_with_PID)

# substitute away the "mRNA.1000.bp" from the ID column, add ID header, also calculate Means column with mean per row

triad_mean_PID <- data.frame(ID=gsub(".mRNA.1000.bp","",triads_with_PID[,1]),Means=rowMeans(triads_with_PID[,-1]))
head(triad_mean_PID)
mean(triad_mean_PID$Means)
median(triad_mean_PID$Means)

hist(triad_mean_PID[,2], xlab="mean percentage ID within triad", main="Promoter (1000 bp upstream of mRNA)
     percentage ID within triads", xlim=c(0,100), ylim=c(0,1200))


# calculate promoter perc ID for conserved vs diverged triads for each stress


#move to folder where 3UP, 3DOWN, opposite group results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\groups_up_down_opposite")
list_of_files_cons_div <- list.files(pattern="*homoeologues.txt")
list_of_files_cons_div
length(list_of_files_cons_div)

#remove the samples with no groups in them
list_of_files_cons_div <- c(list_of_files_cons_div[3:12],list_of_files_cons_div[14:18], 
                            list_of_files_cons_div[21],list_of_files_cons_div[24],list_of_files_cons_div[27],list_of_files_cons_div[30],list_of_files_cons_div[36])
list_of_files_cons_div

results <- c("pattern","mean_perc_id")

i=1

for (i in 1:20) {
  
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\groups_up_down_opposite")
  
  gene_table <- read.table(file=list_of_files_cons_div[i], header=T)
  # read.table(file=list_of_files_cons_div[i])
  head(gene_table)
  dim(gene_table)
  str(gene_table)
  
  merged_gene_table <- merge(gene_table,triad_mean_PID,by.x=0, by.y="ID")
  head(merged_gene_table)
  dim(merged_gene_table)
  mean_perc_id <- mean(merged_gene_table$Means)
  
  results <- rbind(results,c(sub("_homoeologues.txt","",list_of_files_cons_div[i]),mean_perc_id))
  results
  
  setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity\\")
  
  jpeg(file=paste((sub("_homoeologues.txt","",list_of_files_cons_div[i])),"1000bp.jpg"), height=480, width=480, units="px")
  hist(merged_gene_table$Means, xlab="percentage ID within triad", main=paste("Promoter (1000 bp upstream of mRNA)
                                                                              percentage ID within", (sub("_homoeologues.txt","",list_of_files_cons_div[i])), "triads"), xlim=c(0,100))
  dev.off()
  
}  
results

write.table(results,file="3up_3down_opp_mean_perc_id_sim_1000bp_upstream_of_mRNA.txt")


########### for mRNA 1000 bp opposite vs 3 UP vs 3 DOWN merged across stresses ##################

## get promoter perc ID for all triads

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\promoter_similarity")

triads_with_PID <- read.csv(file = "list_of_triads_with_perc_id_1000_mRNA.csv", header=F)
head(triads_with_PID)
dim(triads_with_PID)

# substitute away the "mRNA.1000.bp" from the ID column, add ID header, also calculate Means column with mean per row

triad_mean_PID <- data.frame(ID=gsub(".mRNA.1000.bp","",triads_with_PID[,1]),Means=rowMeans(triads_with_PID[,-1]))
head(triad_mean_PID)
mean(triad_mean_PID$Means)
median(triad_mean_PID$Means)

hist(triad_mean_PID[,2], xlab="mean percentage ID within triad", main="Promoter (1000 bp upstream of mRNA)
     percentage ID within triads", xlim=c(0,100), ylim=c(0,1200))


# calculate promoter perc ID for conserved vs diverged triads for each stress


#move to folder where 3UP, 3DOWN, opposite group results are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001\\groups_up_down_opposite")

up_abiotic <- read.table("up_abiotic.txt")
head(up_abiotic)

down_abiotic <- read.table("down_abiotic.txt")
head(down_abiotic)

opposite_abiotic <- read.table("opposite_abiotic.txt")
head(opposite_abiotic)

# check if overlap between lists
num_triads_in_3UP_list_and_3DOWN_list <- merge(up_abiotic,down_abiotic,by="V1")
num_triads_in_3UP_list_and_3DOWN_list
num_triads_in_3UP_list_and_opposite_list <- merge(up_abiotic,opposite_abiotic,by="V1")
num_triads_in_3UP_list_and_opposite_list
num_triads_in_3DOWN_list_and_opposite_list <- merge(down_abiotic,opposite_abiotic,by="V1")
num_triads_in_3DOWN_list_and_opposite_list

# find out percentage ID of promoters within each list

merged_up_abiotic <- merge(up_abiotic,triad_mean_PID,by.x="V1", by.y="ID")
head(merged_up_abiotic)
dim(merged_up_abiotic)
mean_perc_id_up_abiotic <- mean(merged_up_abiotic$Means)
mean_perc_id_up_abiotic

median_perc_id_up_abiotic <- median(merged_up_abiotic$Means)
median_perc_id_up_abiotic

merged_down_abiotic <- merge(down_abiotic,triad_mean_PID,by.x="V1", by.y="ID")
head(merged_down_abiotic)
dim(merged_down_abiotic)
mean_perc_id_down_abiotic <- mean(merged_down_abiotic$Means)
mean_perc_id_down_abiotic

median_perc_id_down_abiotic <- median(merged_down_abiotic$Means)
median_perc_id_down_abiotic

merged_opposite_abiotic <- merge(opposite_abiotic,triad_mean_PID,by.x="V1", by.y="ID")
head(merged_opposite_abiotic)
dim(merged_opposite_abiotic)
mean_perc_id_opposite_abiotic <- mean(merged_opposite_abiotic$Means)
mean_perc_id_opposite_abiotic

median_perc_id_opposite_abiotic <- median(merged_opposite_abiotic$Means)
median_perc_id_opposite_abiotic

t.test(merged_up_abiotic$Means, merged_down_abiotic$Means, var.equal = T)

t.test(merged_up_abiotic$Means, merged_opposite_abiotic$Means, var.equal = T)

t.test(merged_opposite_abiotic$Means, merged_down_abiotic$Means, var.equal = T)






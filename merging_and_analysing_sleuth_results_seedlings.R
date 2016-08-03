# Merging and analysing results for seedling stress experiment for TGAC genome paper
# Philippa Borrill
# 07.07.2016

## README ##
# uses the output of the scripts:
# sleuth_seedling_stress_indiv_controls_drought_and_heat.R
# sleuth_seedling_stress_indiv_controls_powdery_mildew.R
# sleuth_seedling_stress_indiv_controls_stripe_rust.R

### First need to filter according to desired q-value
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05")

list_of_files <- list.files(pattern="*csv")
list_of_files

filter_stringency <- 0.001

for (k in 1:12) {
  
  data_unfilt <- read.csv(file=list_of_files[k],header=T)
  head(data_unfilt)
  dim(data_unfilt)
  tail(data_unfilt)
  data_filt <- data_unfilt[data_unfilt$qval<filter_stringency,]
  head(data_filt)
  dim(data_filt)
  tail(data_filt)
  data_filt <- data_filt[,2:12]

  write.csv(data_filt,file=paste("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_",filter_stringency,"\\",filter_stringency,list_of_files[k],sep=""))
  
}


### First aim is to make a merged table of the b-fold changes for all the stresses together

# set working directory where the csv files are
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05")
#OR
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.01")
#OR
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.001")



# list the csv files in the directory
list_of_files <- list.files(pattern="*csv")
list_of_files

# colClasses means you don't read in the NULL columns. The NA means that R guesses itself what type of data is in the column

drought_1h <- read.csv(list_of_files[1], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
drought_heat_1h <- read.csv(list_of_files[2], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
heat_1h <- read.csv(list_of_files[3], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
drought_6h <- read.csv(list_of_files[4], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
drought_heat_6h <-read.csv(list_of_files[5], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
heat_6h <- read.csv(list_of_files[6], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
mildew_24h <- read.csv(list_of_files[7], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
mildew_48h <- read.csv(list_of_files[8], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
mildew_72h <- read.csv(list_of_files[9], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
yellow_rust_24h <- read.csv(list_of_files[10], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
yellow_rust_48h <- read.csv(list_of_files[11], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))
yellow_rust_72h <- read.csv(list_of_files[12], colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL"))

head(yellow_rust_72h)

# merge together the datasets
# need to do this step wise adding 1 dataset at a time:

all_datasets <- merge(drought_1h, heat_1h, by="target_id", all=T)
all_datasets <- merge(all_datasets, drought_heat_1h, by="target_id", all=T)
all_datasets <- merge(all_datasets, drought_6h, by="target_id", all=T)
all_datasets <- merge(all_datasets, heat_6h, by="target_id", all=T)
all_datasets <- merge(all_datasets, drought_heat_6h, by="target_id", all=T)
all_datasets <- merge(all_datasets, mildew_24h,  by="target_id", all=T)
all_datasets <- merge(all_datasets, mildew_48h, by="target_id", all=T)
all_datasets <- merge(all_datasets, mildew_72h, by="target_id", all=T)
all_datasets <- merge(all_datasets, yellow_rust_24h, by="target_id", all=T)
all_datasets <- merge(all_datasets, yellow_rust_48h, by="target_id", all=T)
all_datasets <- merge(all_datasets, yellow_rust_72h, by="target_id", all=T)

# write out the new column names for the final all_datasets
column_names_for_all_datasets <- c("target_id","drought_1h", "heat_1h", "drought_heat_1h", "drought_6h", "heat_6h", "drought_heat_6h", "mildew_24h", 
                      "mildew_48h", "mildew_72h", "yellow_rust_24h", "yellow_rust_48h", "yellow_rust_72h")
column_names_for_all_datasets

# replace the meaningless column names with meaningful column names
colnames(all_datasets) <- column_names_for_all_datasets

head(all_datasets)
dim(all_datasets)


# checked that this all worked as expected by grepping the b values in the original files for each sample

# NOTE the b values (fold changes) are in natural log base so any positive value is a fold increase and any negative value is a fold decrease


####### General analysis of differential expresion ########
# want to count the number of numerical values (i.e. differentially expressed genes) in each column
# do this using length(which(!is.na(a)))

# get count of numbers (NOT NA) in column 2
length(which(!is.na(all_datasets[,2])))

# get count of NA in column 2
length(which(is.na(all_datasets[,2])))

# get the count of numbers (differentially expressed genes) in each column i.e. per condition and put them into a table

colnames(all_datasets[,2:13])

num_diff_expr_genes <- c(length(which(!is.na(all_datasets[,2]))), length(which(!is.na(all_datasets[,3]))), length(which(!is.na(all_datasets[,4]))),  
                         length(which(!is.na(all_datasets[,5]))), length(which(!is.na(all_datasets[,6]))), length(which(!is.na(all_datasets[,7]))),
                         length(which(!is.na(all_datasets[,8]))), length(which(!is.na(all_datasets[,9]))), length(which(!is.na(all_datasets[,10]))),
                         length(which(!is.na(all_datasets[,11]))), length(which(!is.na(all_datasets[,12]))), length(which(!is.na(all_datasets[,13]))))
num_diff_expr_genes
num_diff_expr_genes <- rbind(colnames(all_datasets[,2:13]),num_diff_expr_genes)
num_diff_expr_genes

colnames(num_diff_expr_genes) <- num_diff_expr_genes[1,]
num_diff_expr_genes <- num_diff_expr_genes[-1,]
num_diff_expr_genes<- as.data.frame(num_diff_expr_genes)
num_diff_expr_genes
colnames(num_diff_expr_genes)
rownames(num_diff_expr_genes)

# get count of diff expr genes in each column over 2 fold change
log(2)
log(0.5)

num_diff_expr_genes_2fold <- c(length(which((all_datasets[,2]>log(2)))) + length(which((all_datasets[,2]<log(0.5)))), 
                               length(which((all_datasets[,3]>log(2)))) + length(which((all_datasets[,3]<log(0.5)))),
                               length(which((all_datasets[,4]>log(2)))) + length(which((all_datasets[,4]<log(0.5)))),
                               length(which((all_datasets[,5]>log(2)))) + length(which((all_datasets[,5]<log(0.5)))),
                               length(which((all_datasets[,6]>log(2)))) + length(which((all_datasets[,6]<log(0.5)))),
                               length(which((all_datasets[,7]>log(2)))) + length(which((all_datasets[,7]<log(0.5)))),
                               length(which((all_datasets[,8]>log(2)))) + length(which((all_datasets[,8]<log(0.5)))),
                               length(which((all_datasets[,9]>log(2)))) + length(which((all_datasets[,9]<log(0.5)))),
                               length(which((all_datasets[,10]>log(2)))) + length(which((all_datasets[,10]<log(0.5)))),
                               length(which((all_datasets[,11]>log(2)))) + length(which((all_datasets[,11]<log(0.5)))),
                               length(which((all_datasets[,12]>log(2)))) + length(which((all_datasets[,12]<log(0.5)))),
                               length(which((all_datasets[,13]>log(2)))) + length(which((all_datasets[,13]<log(0.5)))))
                               
num_diff_expr_genes_2fold <- rbind(colnames(all_datasets[,2:13]),num_diff_expr_genes_2fold)
num_diff_expr_genes_2fold

# plot the number of genes differentially expressed in each sample
barplot(t(num_diff_expr_genes))

##### want to plot a graph with up and downreg genes - so need to count number of up and down reg genes in each condition#####
upreg_genes <- c(length(which((all_datasets[,2])>0)), length(which((all_datasets[,3])>0)), length(which((all_datasets[,4])>0)),  
                         length(which((all_datasets[,5])>0)), length(which((all_datasets[,6])>0)), length(which((all_datasets[,7])>0)),
                         length(which((all_datasets[,8])>0)), length(which((all_datasets[,9])>0)), length(which((all_datasets[,10])>0)),
                         length(which((all_datasets[,11])>0)), length(which((all_datasets[,12])>0)), length(which((all_datasets[,13])>0)))
upreg_genes

downreg_genes <- c(length(which((all_datasets[,2])<0)), length(which((all_datasets[,3])<0)), length(which((all_datasets[,4])<0)),  
                   length(which((all_datasets[,5])<0)), length(which((all_datasets[,6])<0)), length(which((all_datasets[,7])<0)),
                   length(which((all_datasets[,8])<0)), length(which((all_datasets[,9])<0)), length(which((all_datasets[,10])<0)),
                   length(which((all_datasets[,11])<0)), length(which((all_datasets[,12])<0)), length(which((all_datasets[,13])<0)))

downreg_genes

upreg_genes+downreg_genes

# make dataframe including names and up and down reg genes
diff_expr_genes <- rbind(colnames(all_datasets[,2:13]),upreg_genes,downreg_genes)
diff_expr_genes
colnames(diff_expr_genes) <- diff_expr_genes[1,]
diff_expr_genes <- diff_expr_genes[-1,]
diff_expr_genes<- as.matrix(diff_expr_genes)
diff_expr_genes
diff_expr_genes_transposed <- as.data.frame(t(diff_expr_genes))
diff_expr_genes_transposed
rownames(diff_expr_genes_transposed)

write.table(diff_expr_genes_transposed,"diff_expressed_genes_q_0.001_all_conditions.txt")


###### want to plot a graph with 2-fold up and downreg genes - so need to count number of up and down reg genes in each condition#####
upreg_genes_2fold <- c(length(which((all_datasets[,2])>log(2))), length(which((all_datasets[,3])>log(2))), length(which((all_datasets[,4])>log(2))),  
                 length(which((all_datasets[,5])>log(2))), length(which((all_datasets[,6])>log(2))), length(which((all_datasets[,7])>log(2))),
                 length(which((all_datasets[,8])>log(2))), length(which((all_datasets[,9])>log(2))), length(which((all_datasets[,10])>log(2))),
                 length(which((all_datasets[,11])>log(2))), length(which((all_datasets[,12])>log(2))), length(which((all_datasets[,13])>log(2))))
upreg_genes_2fold

downreg_genes_2fold <- c(length(which((all_datasets[,2])<log(0.5))), length(which((all_datasets[,3])<log(0.5))), length(which((all_datasets[,4])<log(0.5))),  
                   length(which((all_datasets[,5])<log(0.5))), length(which((all_datasets[,6])<log(0.5))), length(which((all_datasets[,7])<log(0.5))),
                   length(which((all_datasets[,8])<log(0.5))), length(which((all_datasets[,9])<log(0.5))), length(which((all_datasets[,10])<log(0.5))),
                   length(which((all_datasets[,11])<log(0.5))), length(which((all_datasets[,12])<log(0.5))), length(which((all_datasets[,13])<log(0.5))))

downreg_genes_2fold

upreg_genes_2fold+downreg_genes_2fold

# make dataframe including names and up and down reg genes
diff_expr_genes_2fold <- rbind(colnames(all_datasets[,2:13]),upreg_genes_2fold,downreg_genes_2fold)
diff_expr_genes_2fold
colnames(diff_expr_genes_2fold) <- diff_expr_genes_2fold[1,]
diff_expr_genes_2fold <- diff_expr_genes_2fold[-1,]
diff_expr_genes_2fold<- as.matrix(diff_expr_genes_2fold)
diff_expr_genes_2fold
diff_expr_genes_2fold_transposed <- as.data.frame(t(diff_expr_genes_2fold))
diff_expr_genes_2fold_transposed
rownames(diff_expr_genes_2fold_transposed)

write.table(diff_expr_genes_2fold_transposed,"diff_expressed_genes_q_0.001_2fold_all_conditions.txt")



##### identify which gene triads have similar or different expression between homoeologues ###########

# read in list of triads

triad_info <- read.table("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\sleuth_results_q_value_0.05\\triads_over_90_percent.txt",sep="\t", header=T)
head(triad_info)
dim(triad_info)

# get separate lists of A, B and D genome genes within triads
A_genome_genes <- as.data.frame(triad_info[,2])
head(A_genome_genes)
B_genome_genes <- as.data.frame(triad_info[,3])
D_genome_genes <- as.data.frame(triad_info[,4])

colnames(A_genome_genes) <- "target_id"
colnames(B_genome_genes) <- "target_id"
colnames(D_genome_genes) <- "target_id"

head(A_genome_genes)
head(B_genome_genes)
head(D_genome_genes)

# merge the A_genome_gene list with all_datasets, 
# use all.x=T to only keep expression for genes in the A_genome_gene list
A_genome_expression <- merge(A_genome_genes, all_datasets, by="target_id", all.x=T)
head(A_genome_expression)
tail(A_genome_expression)
dim(A_genome_genes)
dim(A_genome_expression)

# count 2 fold upreg and 2fold down reg A genome genes

upreg_genes_2fold_A <- c(length(which((A_genome_expression[,2])>log(2))), length(which((A_genome_expression[,3])>log(2))), length(which((A_genome_expression[,4])>log(2))),  
                       length(which((A_genome_expression[,5])>log(2))), length(which((A_genome_expression[,6])>log(2))), length(which((A_genome_expression[,7])>log(2))),
                       length(which((A_genome_expression[,8])>log(2))), length(which((A_genome_expression[,9])>log(2))), length(which((A_genome_expression[,10])>log(2))),
                       length(which((A_genome_expression[,11])>log(2))), length(which((A_genome_expression[,12])>log(2))), length(which((A_genome_expression[,13])>log(2))))
upreg_genes_2fold_A

downreg_genes_2fold_A <- c(length(which((A_genome_expression[,2])<log(0.5))), length(which((A_genome_expression[,3])<log(0.5))), length(which((A_genome_expression[,4])<log(0.5))),  
                         length(which((A_genome_expression[,5])<log(0.5))), length(which((A_genome_expression[,6])<log(0.5))), length(which((A_genome_expression[,7])<log(0.5))),
                         length(which((A_genome_expression[,8])<log(0.5))), length(which((A_genome_expression[,9])<log(0.5))), length(which((A_genome_expression[,10])<log(0.5))),
                         length(which((A_genome_expression[,11])<log(0.5))), length(which((A_genome_expression[,12])<log(0.5))), length(which((A_genome_expression[,13])<log(0.5))))

downreg_genes_2fold_A
colnames(all_datasets[2:13])
by_homoeologue_2_fold <- rbind(colnames(all_datasets[2:13]), upreg_genes_2fold_A, downreg_genes_2fold_A)
by_homoeologue_2_fold
colnames(by_homoeologue_2_fold) <- by_homoeologue_2_fold[1,]
by_homoeologue_2_fold <- by_homoeologue_2_fold[-1,]
by_homoeologue_2_fold
barplot(by_homoeologue_2_fold, beside=TRUE)


# merge the B_genome_gene list with all_datasets, 
# use all.x=T to only keep expression for genes in the B_genome_gene list
B_genome_expression <- merge(B_genome_genes, all_datasets, by="target_id", all.x=T)
head(B_genome_expression)
tail(B_genome_expression)
dim(B_genome_genes)
dim(B_genome_expression)

# count 2 fold upreg and 2fold down reg B genome genes

upreg_genes_2fold_B <- c(length(which((B_genome_expression[,2])>log(2))), length(which((B_genome_expression[,3])>log(2))), length(which((B_genome_expression[,4])>log(2))),  
                         length(which((B_genome_expression[,5])>log(2))), length(which((B_genome_expression[,6])>log(2))), length(which((B_genome_expression[,7])>log(2))),
                         length(which((B_genome_expression[,8])>log(2))), length(which((B_genome_expression[,9])>log(2))), length(which((B_genome_expression[,10])>log(2))),
                         length(which((B_genome_expression[,11])>log(2))), length(which((B_genome_expression[,12])>log(2))), length(which((B_genome_expression[,13])>log(2))))
upreg_genes_2fold_B

downreg_genes_2fold_B <- c(length(which((B_genome_expression[,2])<log(0.5))), length(which((B_genome_expression[,3])<log(0.5))), length(which((B_genome_expression[,4])<log(0.5))),  
                           length(which((B_genome_expression[,5])<log(0.5))), length(which((B_genome_expression[,6])<log(0.5))), length(which((B_genome_expression[,7])<log(0.5))),
                           length(which((B_genome_expression[,8])<log(0.5))), length(which((B_genome_expression[,9])<log(0.5))), length(which((B_genome_expression[,10])<log(0.5))),
                           length(which((B_genome_expression[,11])<log(0.5))), length(which((B_genome_expression[,12])<log(0.5))), length(which((B_genome_expression[,13])<log(0.5))))

downreg_genes_2fold_B

by_homoeologue_2_fold <- rbind(by_homoeologue_2_fold, upreg_genes_2fold_B, downreg_genes_2fold_B)
by_homoeologue_2_fold

# merge the D_genome_gene list with all_datasets, 
# use all.x=T to only keep expression for genes in the D_genome_gene list
D_genome_expression <- merge(D_genome_genes, all_datasets, by="target_id", all.x=T)
head(D_genome_expression)
tail(D_genome_expression)
dim(D_genome_genes)
dim(D_genome_expression)

# count 2 fold upreg and 2fold down reg D genome genes

upreg_genes_2fold_D <- c(length(which((D_genome_expression[,2])>log(2))), length(which((D_genome_expression[,3])>log(2))), length(which((D_genome_expression[,4])>log(2))),  
                         length(which((D_genome_expression[,5])>log(2))), length(which((D_genome_expression[,6])>log(2))), length(which((D_genome_expression[,7])>log(2))),
                         length(which((D_genome_expression[,8])>log(2))), length(which((D_genome_expression[,9])>log(2))), length(which((D_genome_expression[,10])>log(2))),
                         length(which((D_genome_expression[,11])>log(2))), length(which((D_genome_expression[,12])>log(2))), length(which((D_genome_expression[,13])>log(2))))
upreg_genes_2fold_D

downreg_genes_2fold_D <- c(length(which((D_genome_expression[,2])<log(0.5))), length(which((D_genome_expression[,3])<log(0.5))), length(which((D_genome_expression[,4])<log(0.5))),  
                           length(which((D_genome_expression[,5])<log(0.5))), length(which((D_genome_expression[,6])<log(0.5))), length(which((D_genome_expression[,7])<log(0.5))),
                           length(which((D_genome_expression[,8])<log(0.5))), length(which((D_genome_expression[,9])<log(0.5))), length(which((D_genome_expression[,10])<log(0.5))),
                           length(which((D_genome_expression[,11])<log(0.5))), length(which((D_genome_expression[,12])<log(0.5))), length(which((D_genome_expression[,13])<log(0.5))))

downreg_genes_2fold_D

by_homoeologue_2_fold <- rbind(by_homoeologue_2_fold, upreg_genes_2fold_D, downreg_genes_2fold_D)
by_homoeologue_2_fold

write.csv(by_homoeologue_2_fold,file="triad_genes_up_down_reg_2_fold_q0.001_split_by_homoeologue.csv")

# add group info (i.e. which triad the gene belongs to) 
# as the rownames in each of A, B and D _genome_expression

rownames(A_genome_expression) <- triad_info[,1] 
head(rownames(A_genome_expression))

rownames(B_genome_expression) <- triad_info[,1] 
head(rownames(B_genome_expression))

rownames(D_genome_expression) <- triad_info[,1] 
head(rownames(D_genome_expression))

# now make separate dataframes which contain the A, B and D genome expression for each condition

head(B_genome_expression)

D1h_homoeologues <- cbind(A_genome_expression[,2],B_genome_expression[,2],D_genome_expression[,2]) 
rownames(D1h_homoeologues) <- triad_info[,1]
colnames(D1h_homoeologues) <- c("A","B","D")

#check it looks ok:
head(D1h_homoeologues)

head(A_genome_expression)
head(B_genome_expression)
head(D_genome_expression)

# do the same for all the conditions:

H1h_homoeologues <- cbind(A_genome_expression[,3],B_genome_expression[,3],D_genome_expression[,3]) 
rownames(H1h_homoeologues) <- triad_info[,1]
colnames(H1h_homoeologues) <- c("A","B","D")

DH1h_homoeologues <- cbind(A_genome_expression[,4],B_genome_expression[,4],D_genome_expression[,4]) 
rownames(DH1h_homoeologues) <- triad_info[,1]
colnames(DH1h_homoeologues) <- c("A","B","D")

D6h_homoeologues <- cbind(A_genome_expression[,5],B_genome_expression[,5],D_genome_expression[,5]) 
rownames(D6h_homoeologues) <- triad_info[,1]
colnames(D6h_homoeologues) <- c("A","B","D")

H6h_homoeologues <- cbind(A_genome_expression[,6],B_genome_expression[,6],D_genome_expression[,6]) 
rownames(H6h_homoeologues) <- triad_info[,1]
colnames(H6h_homoeologues) <- c("A","B","D")

DH6h_homoeologues <- cbind(A_genome_expression[,7],B_genome_expression[,7],D_genome_expression[,7]) 
rownames(DH6h_homoeologues) <- triad_info[,1]
colnames(DH6h_homoeologues) <- c("A","B","D")

PM24h_homoeologues <- cbind(A_genome_expression[,8],B_genome_expression[,8],D_genome_expression[,8]) 
rownames(PM24h_homoeologues) <- triad_info[,1]
colnames(PM24h_homoeologues) <- c("A","B","D")

PM48h_homoeologues <- cbind(A_genome_expression[,9],B_genome_expression[,9],D_genome_expression[,9]) 
rownames(PM48h_homoeologues) <- triad_info[,1]
colnames(PM48h_homoeologues) <- c("A","B","D")

PM72h_homoeologues <- cbind(A_genome_expression[,10],B_genome_expression[,10],D_genome_expression[,10]) 
rownames(PM72h_homoeologues) <- triad_info[,1]
colnames(PM72h_homoeologues) <- c("A","B","D")

SR24h_homoeologues <- cbind(A_genome_expression[,11],B_genome_expression[,11],D_genome_expression[,11]) 
rownames(SR24h_homoeologues) <- triad_info[,1]
colnames(SR24h_homoeologues) <- c("A","B","D")

SR48h_homoeologues <- cbind(A_genome_expression[,12],B_genome_expression[,12],D_genome_expression[,12]) 
rownames(SR48h_homoeologues) <- triad_info[,1]
colnames(SR48h_homoeologues) <- c("A","B","D")

SR72h_homoeologues <- cbind(A_genome_expression[,13],B_genome_expression[,13],D_genome_expression[,13]) 
rownames(SR72h_homoeologues) <- triad_info[,1]
colnames(SR72h_homoeologues) <- c("A","B","D")


head(SR72h_homoeologues)
tail(SR72h_homoeologues)
dim(SR72h_homoeologues)

# Want to categorise the triads in each condition in three categories:
## same expression change 
#e.g up up up, down down down and flat flat flat for A, B and D respectively
## similar expression change 
#e.g. up flat flat, down down flat, flat down flat etc. (e.g. flat and down or flat and up exclusively)
## opposite expression change 
#e.g. up down up, up up down, down up flat, flat down up etc. (i.e. always contain at least 1 up and at least 1 down)

# to do this need to categorise expression change in each vector as up, down or flat.
# cut-offs will be up is over 2 fold upreg, down is 0.5 fold change and flat is between 0.5 to 2 fold.
# The values in the table are natural logged so need to be careful
# 2 fold change in natural log is:
log(2)
 #= 0.6931472
log(0.5)
#=-0.6931472

# Will do this for each of the 12 conditions. Will start as a trial with D1h_homoeologues
# need to convert to data frame
D1h_homoeologues <- as.data.frame(D1h_homoeologues)

# this statement says that if the value is over log(2) write UP, if not see if it is below log(0.5)- if it is write DOWN
# Otherwise write FLAT
D1h_homoeologues_recoded <- ifelse(D1h_homoeologues > log(2), c("UP"), ifelse(D1h_homoeologues < log(0.5),c("DOWN"), c("FLAT"))) 
head(D1h_homoeologues)
head(D1h_homoeologues_recoded)

D1h_homoeologues_recoded <- as.data.frame(D1h_homoeologues_recoded)
# get the rows which are the same in A, B and D
D1h_homoeologues_recoded_same <- D1h_homoeologues_recoded[ which(D1h_homoeologues_recoded$A == D1h_homoeologues_recoded$B & D1h_homoeologues_recoded$B == D1h_homoeologues_recoded$D),]
head(D1h_homoeologues_recoded)
head(D1h_homoeologues_recoded_same)

# get rows which are opposite in A, B and D
D1h_homoeologues_recoded_opposite <- 
  D1h_homoeologues_recoded[ which ( (D1h_homoeologues_recoded$A == 'UP' | D1h_homoeologues_recoded$B == 'UP'| D1h_homoeologues_recoded$D == 'UP') & (D1h_homoeologues_recoded$A == 'DOWN' | D1h_homoeologues_recoded$B == 'DOWN'| D1h_homoeologues_recoded$D == 'DOWN') ),]
head(D1h_homoeologues_recoded_opposite)
dim(D1h_homoeologues_recoded_opposite)

# get rows which are UP reg in all 3 homoeologues
D1h_homoeologues_recoded_up <- 
  D1h_homoeologues_recoded[ which ( (D1h_homoeologues_recoded$A == 'UP' & D1h_homoeologues_recoded$B == 'UP'& D1h_homoeologues_recoded$D == 'UP')  ),]
head(D1h_homoeologues_recoded_up)
dim(D1h_homoeologues_recoded_up)

# get rows which are DOWN reg in all 3 homoeologues
D1h_homoeologues_recoded_down <- 
  D1h_homoeologues_recoded[ which ( (D1h_homoeologues_recoded$A == 'DOWN' & D1h_homoeologues_recoded$B == 'DOWN'& D1h_homoeologues_recoded$D == 'DOWN')  ),]
head(D1h_homoeologues_recoded_down)
dim(D1h_homoeologues_recoded_down)


D1h_results <- c(length(D1h_homoeologues_recoded_opposite$A),length(D1h_homoeologues_recoded_up$A), length(D1h_homoeologues_recoded_down$A))
D1h_results



## NOW want to make a loop to go through each of the conditions

## list of conditions with homoeologue expression

list_cond_homoeologues <- list(a=D1h_homoeologues,b=H1h_homoeologues,c=DH1h_homoeologues,
                            d=D6h_homoeologues,e=H6h_homoeologues,f=DH6h_homoeologues,
                            g=PM24h_homoeologues,h=PM48h_homoeologues,i=PM72h_homoeologues,
                            j=SR24h_homoeologues,k=SR48h_homoeologues,l=SR72h_homoeologues)


#list_cond_homoeologues <- list(D1h_homoeologues=D1h_homoeologues)
#list_cond_homoeologues
# Will do this for each of the 12 conditions. Will start as a trial with D1h_homoeologues
# need to convert to data frame

# set up empty results vector
results <- c("opposite", "All_3_UP", "All_3_DOWN", "2_UP", "2_DOWN", "UP or DOWN 1 homoeologue")
for (i in 1:12) {
  print(i)

homoeologues <- list_cond_homoeologues[[i]]  
test <- head(homoeologues)
print(test)


# this statement says that if the value is over log(2) write UP, if not see if it is below log(0.5)- if it is write DOWN
# Otherwise write FLAT
homoeologues_recoded <- ifelse(homoeologues > log(2), c("UP"), ifelse(homoeologues < log(0.5),c("DOWN"), c("FLAT"))) 
#homoeologues_recoded <- ifelse(homoeologues > 0, c("UP"), ifelse(homoeologues < 0,c("DOWN"), c("FLAT"))) 
head(homoeologues)
head(homoeologues_recoded)

homoeologues_recoded <- as.data.frame(homoeologues_recoded)
dim(homoeologues_recoded)
# get rows which are opposite in A, B and D
homoeologues_recoded_opposite <- 
  as.data.frame(homoeologues_recoded[ which ( (homoeologues_recoded$A == 'UP' | homoeologues_recoded$B == 'UP'| homoeologues_recoded$D == 'UP') & (homoeologues_recoded$A == 'DOWN' | homoeologues_recoded$B == 'DOWN'| homoeologues_recoded$D == 'DOWN') ),])


# get rows which are UP reg in all 3 homoeologues
homoeologues_recoded_up <- 
  as.data.frame(homoeologues_recoded[ which ( (homoeologues_recoded$A == 'UP' & homoeologues_recoded$B == 'UP'& homoeologues_recoded$D == 'UP')  ),])

# get rows which are DOWN reg in all 3 homoeologues
homoeologues_recoded_down <- 
  as.data.frame(homoeologues_recoded[ which ( (homoeologues_recoded$A == 'DOWN' & homoeologues_recoded$B == 'DOWN'& homoeologues_recoded$D == 'DOWN')  ),])

# get rows which are DOWN reg in 2 out of 3 homoeologues (and the third one isn't UP)
homoeologues_recoded_down_in_2 <- 
  as.data.frame(homoeologues_recoded[ which ( (homoeologues_recoded$A == 'DOWN' & homoeologues_recoded$B == 'DOWN'& (homoeologues_recoded$D != 'UP'| is.na(homoeologues_recoded$D)) | 
                                                 homoeologues_recoded$A == 'DOWN' & homoeologues_recoded$D == 'DOWN' & (homoeologues_recoded$B != 'UP'| is.na(homoeologues_recoded$B)) |
                                                 homoeologues_recoded$B == 'DOWN' & homoeologues_recoded$D == 'DOWN' & (homoeologues_recoded$A != 'UP'| is.na(homoeologues_recoded$A)))  ),])

# get rows which are UP reg in 2 out of 3 homoeologues (and the third one isn't DOWN)
homoeologues_recoded_up_in_2 <- 
  as.data.frame(homoeologues_recoded[ which ( (homoeologues_recoded$A == 'UP' & homoeologues_recoded$B == 'UP'& (homoeologues_recoded$D != 'DOWN'| is.na(homoeologues_recoded$D)) | 
                                                 homoeologues_recoded$A == 'UP' & homoeologues_recoded$D == 'UP' & (homoeologues_recoded$B != 'DOWN'| is.na(homoeologues_recoded$B)) |
                                                 homoeologues_recoded$B == 'UP' & homoeologues_recoded$D == 'UP' & (homoeologues_recoded$A != 'DOWN'| is.na(homoeologues_recoded$A)))  ),])


# get rows which are UP or DOWN reg in at least 1 homoeologue
homoeologues_recoded_UP_or_DOWN_1_homoeologue <- 
  as.data.frame(homoeologues_recoded[ which ( (homoeologues_recoded$A == 'DOWN' | 
             homoeologues_recoded$B == 'DOWN'| homoeologues_recoded$D == 'DOWN' | 
             homoeologues_recoded$A == 'UP' | homoeologues_recoded$B == 'UP'| homoeologues_recoded$D == 'UP')) ,])

results <- rbind (results,(c(length(homoeologues_recoded_opposite$A),length(homoeologues_recoded_up$A), length(homoeologues_recoded_down$A),
                             length(homoeologues_recoded_up_in_2$A), length(homoeologues_recoded_down_in_2$A) , length(homoeologues_recoded_UP_or_DOWN_1_homoeologue$A) )))
results

triad_names_conserved_expr <- rbind(homoeologues_recoded_up,homoeologues_recoded_down,homoeologues_recoded_up_in_2,homoeologues_recoded_down_in_2)
head(triad_names_conserved_expr)
dim(triad_names_conserved_expr)

# add in gene info for each triad for conserved genes
merged_triads_conserved_expr <- merge(triad_names_conserved_expr,triad_info, by.x = 'row.names', by.y = "Group")
head(merged_triads_conserved_expr)

colnames(merged_triads_conserved_expr) <- c("Group",colnames(merged_triads_conserved_expr)[2:7])
colnames(merged_triads_conserved_expr)

# add in gene info for each triad for divergent genes
triad_names_divergent_expr <- homoeologues_recoded_opposite

merged_triads_divergent_expr <- merge(triad_names_divergent_expr,triad_info, by.x = 'row.names', by.y = "Group")
head(merged_triads_divergent_expr)
colnames(merged_triads_divergent_expr) <- c("Group",colnames(merged_triads_divergent_expr)[2:7])
colnames(merged_triads_divergent_expr)


# need to specify the names of the 12 conditions to write the output file names correctly
row_names_for_results <- c("drought_1h", "heat_1h", "drought_heat_1h", "drought_6h", "heat_6h", "drought_heat_6h", "mildew_24h", 
                           "mildew_48h", "mildew_72h", "yellow_rust_24h", "yellow_rust_48h", "yellow_rust_72h")

# write an output file for each conserved and divergent list for each condition
write.table(merged_triads_conserved_expr,file=paste(row_names_for_results[i],"_triad_conserved_expr.txt"))
write.table(merged_triads_divergent_expr, file=paste(row_names_for_results[i],"_triad_divergent_expr.txt"))

# write an output file for 3 UP, 3 DOwn and opposite for each condition
write.table(homoeologues_recoded_up, file= paste(row_names_for_results[i],"_3_up_homoeologues.txt"))
write.table(homoeologues_recoded_down, file= paste(row_names_for_results[i],"_3_down_homoeologues.txt"))
write.table(homoeologues_recoded_opposite, file= paste(row_names_for_results[i],"_opposite_homoeologues.txt"))

}
colnames(results) <- results[1,]
results <- results[-1,]
row_names_for_results <- c("drought_1h", "heat_1h", "drought_heat_1h", "drought_6h", "heat_6h", "drought_heat_6h", "mildew_24h", 
                                   "mildew_48h", "mildew_72h", "yellow_rust_24h", "yellow_rust_48h", "yellow_rust_72h")
rownames(results) <- row_names_for_results
results

results_log2_log0.5 <- results
results_log2_log0.5

#results_log1.5_log2_3 <- results
results_log1.5_log2_3

barplot(results_log2_log0.5)
head(results_log2_log0.5)
rownames(results_log2_log0.5)
results_log2_log0.5 <- as.data.frame(results_log2_log0.5)
head(results_log2_log0.5)

write.table(results_log2_log0.5, file="diff_expr_pattern_results_log2_log0.5.txt")

# melt data for ggplot
library(reshape2)
library(RColorBrewer)

results_log2_log0.5_with_rownames <- cbind(rownames(results_log2_log0.5), results_log2_log0.5)
head(results_log2_log0.5_with_rownames)
tail(results_log2_log0.5_with_rownames)

results_log2_log0.5_with_rownames <- as.data.frame(results_log2_log0.5_with_rownames)
str(results_log2_log0.5_with_rownames)

# need to convert "factors" to numeric
colnames(results_log2_log0.5_with_rownames)
results_log2_log0.5_with_rownames$opposite <- as.numeric(as.character(results_log2_log0.5_with_rownames$opposite))
results_log2_log0.5_with_rownames$UP <- as.numeric(as.character(results_log2_log0.5_with_rownames$UP))
results_log2_log0.5_with_rownames$DOWN <- as.numeric(as.character(results_log2_log0.5_with_rownames$DOWN))
results_log2_log0.5_with_rownames$`UP or DOWN 1 homoeologue` <- as.numeric(as.character(results_log2_log0.5_with_rownames$`UP or DOWN 1 homoeologue`))
str(results_log2_log0.5_with_rownames)


melted_results_log2_log0.5 <- melt(results_log2_log0.5_with_rownames, id=c("rownames(results_log2_log0.5)"))
head(melted_results_log2_log0.5)
tail(melted_results_log2_log0.5)

str(melted_results_log2_log0.5)

colnames(melted_results_log2_log0.5) <- c("condition", "expr_pattern", "number_of_triads")
  colnames(melted_results_log2_log0.5)
  head(melted_results_log2_log0.5)

library(scales)

pdf(file="barplot_opp_dir.pdf") # have to do as pdf because eps doesn't support transparent
ggplot(melted_results_log2_log0.5, aes(expr_pattern,number_of_triads)) + geom_bar(stat="identity")  + facet_wrap(~condition) 
dev.off()




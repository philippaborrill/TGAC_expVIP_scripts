# Get tpms for specific samples for David and Luca
# Philippa Borrill
# 18.9.2015


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

# select just columns from Choulet paper
(colnames(tpm[112:135]))
colnames(tpm[413:418])


# select just columns from stress/disease
(colnames(tpm[159:193]))


stress_plus_choulet_data <- cbind(tpm[,112:135], tpm[,413:418], tpm[,159:193])
head(stress_plus_choulet_data)
dim(stress_plus_choulet_data)

# replace names of columns to be meaningful
colnames(stress_plus_choulet_data[1:30])

names <- c("stripe_rust_24h_rep1","stripe_rust_24h_rep2","stripe_rust_24h_rep3",
  "stripe_rust_48h_rep1","stripe_rust_48h_rep2","stripe_rust_48h_rep3",
  "stripe_rust_72h_rep1","stripe_rust_72h_rep2","stripe_rust_72h_rep3",
  "control_for_stripe_rust_and_powdery_mildew_rep1", "control_for_stripe_rust_and_powdery_mildew_rep2", "control_for_stripe_rust_and_powdery_mildew_rep3",
  "powdery_mildew_24h_rep1","powdery_mildew_24h_rep2","powdery_mildew_24h_rep3",
  "powdery_mildew_48h_rep1","powdery_mildew_48h_rep2","powdery_mildew_48h_rep3",
  "powdery_mildew_72h_rep1","powdery_mildew_72h_rep2","powdery_mildew_72h_rep3",
  "control_for_drought_and_or_heat_stress_rep1", "control_for_drought_and_or_heat_stress_rep2", 
  "1h_drought_rep1","1h_drought_rep2",
  "6h_drought_rep1","6h_drought_rep2",
  "1h_heat_rep1","1h_heat_rep2",
  "6h_heat_rep1","6h_heat_rep2",
  "1h_drought_and_heat_rep1","1h_drought_and_heat_rep2",
  "6h_drought_and_heat_rep1","6h_drought_and_heat_rep2")

names
colnames(stress_plus_choulet_data) <- c(colnames(stress_plus_choulet_data[1:30]),names)
colnames(stress_plus_choulet_data)

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\for_David_and_Luca")
write.table(stress_plus_choulet_data, file="stress_plus_choulet_data.txt")



# select only genes which Fu Hao wants

# edited EOL to be unix
# use grep on the command line (i.e. go to the cluster then run these commands)

#grep -f 247.gene.ids.for.synteny.txt choulet_data.txt > 247_genes_in_choulet_data.txt

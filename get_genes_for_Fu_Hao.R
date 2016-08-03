# Plot gene expression data for manuscript
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

choulet_data <- cbind(tpm[,112:135], tpm[,413:418])
head(choulet_data)
dim(choulet_data)

setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\for_Fu_Hao")
write.table(choulet_data,file="choulet_data.txt")

# select just columns which are Chinese Spring with 3 reps
colnames(tpm[1:12])
colnames(tpm[73:102])

replicated_data <- cbind(tpm[,1:12], tpm[,73:102])
head(replicated_data)
dim(replicated_data)
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\for_Fu_Hao")
write.table(replicated_data,file="replicated_data.txt")


# select only genes which Fu Hao wants

# edited EOL to be unix
# use grep on the command line (i.e. go to the cluster then run these commands)

#grep -f 247.gene.ids.for.synteny.txt choulet_data.txt > 247_genes_in_choulet_data.txt

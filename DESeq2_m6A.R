#If running DESeq2 on Aquila, type:
#Rscript-3.5.0 DESeq2_m6A.R
#Remember to change sample name twice!

#For parallelization

library("BiocParallel")
register(MulticoreParam(4))

#### Loading counts table (tab-delimited, header row, row name (gene) in column 1

#Input sample name
genotype = "Mettl3-KO-S2"
#Input version here. Either "all" or "Aonly"
version = "Aonly"


inputname = paste("DESeq2/",genotype,"_for_DESeq2_",version,".txt",sep="")
counts = read.table(inputname,header=T,row.names=1,sep="\t")

#### Set up parameters of the experiment (assuming each column is a sample)

## List the genotype for each sample

condition = c("no","yes","no","yes")#"no",#"yes")

## List the batch for each sample. You can add additional information such as
#   tissue etc with additional vectors
##

batch = c("b1","b1","b2","b2")#,"b3","b3")


## Combine the multiple conditions/expt parameters into a single table

colInfo = data.frame(condition,batch,row.names=colnames(counts))

#### Load DESeq2
#BiocManager::install("DESeq2")
library(DESeq2)

#### Create a DESeq2 dataset using the count table and the expt parameter table

## In this example, we are just comparing samples based on their genotype 

######dds = DESeqDataSetFromMatrix(countData = counts, colData = colInfo, design = ~ condition)

## Explicitly state the "control" for DESeq2

######dds$condition = relevel(dds$condition,"no")

#### Run DESeq2
######dds = DESeq(dds)

#### Obtaining results
## For all results, regardless of significance

######res = results(dds, alpha=1.0)

## For results with an FDR < 0.05. Default "alpha" is 0.1

# res = results(dds,alpha=0.1)

## For filtering results based on log FC

# res = results(dds,lfcThreshold = 1)

#### Output results

######write.table(res,"store/DESeq2/NoVsYes_DESeq2_results.txt",sep="\t",quote=F,row.names=T,col.names=T)


#### Multi-factor analysis (DO THIS!!!)

## This is if you wish to "account" for additional variables that might have
#   a simple analysis. E.g. batch effect. It uses a generalized linear model.
#   It is the simplest approach, but there are other algorithms that might be
#   more robust (but also more complicated)

#### Assuming that you are using the previous dataset, redefining the comparison

## Put the variable that you want to test for at the end, and "confounding" variables before it


## If starting from scratch, use this command

ddss = DESeqDataSetFromMatrix(countData = counts, colData = colInfo, design = ~ batch + condition)

ddss$condition = relevel(ddss$condition,"no")


#### Run DESeq2 as before
ddss = DESeq(ddss,parallel=T)

res = results(ddss, alpha=0.1)

outputname = paste("DESeq2/",genotype,"_NoVsYes_DESeq2_results_condition_and_batch_",version,".txt",sep="")
write.table(res,outputname,sep="\t",quote=F,row.names=T)

# This script is for running DESeq as a part of the differential expression analysis comparing 
# the arils from the Durian cultivars Musang King and Monthing. The script is also for presenting 
# the results in useful ways. 

library( "DESeq2" )

# Musang King Aril
SRR6040094 = read.table(file="/Users/Elin/Documents/Uni/Genomanalys/Labs/SRR6040094_count.txt")
SRR6040097 = read.table(file="/Users/Elin/Documents/Uni/Genomanalys/Labs/SRR6040097_count.txt")
# Monthong
SRR6156066 = read.table(file="/Users/Elin/Documents/Uni/Genomanalys/Labs/SRR6156066_count.txt")
SRR6156067 = read.table(file="/Users/Elin/Documents/Uni/Genomanalys/Labs/SRR6156067_count.txt")

#Create the table of the data for all samples
# Musang King Aril
full_table <- SRR6040094[2]
names(full_table)[1] <- "SRR6040094"
full_table$SRR6040097 <- SRR6040097[,2]
# Monthong Aril
full_table$SRR6156066 <- SRR6156066[,2]
full_table$SRR6156067 <- SRR6156067[,2]
rownames(full_table) <-SRR6040094[,1]

# Remove last rows 
full_table<-full_table[1:(nrow(full_table)-5),]

# Create metadata table
run_ID =c('SRR6040094', 'SRR6040097', 'SRR6156066', 'SRR6156067')
# 1 = Musang King, 2 = Monthong
Sample_ID = c('1', '1', '2', '2')
metaData = data.frame(run_ID, Sample_ID)

#-----Deseq2-------------------------------------

# Create the count table with the Sample_ID as a differentiator
dds <- DESeqDataSetFromMatrix(
  countData = full_table,
  colData = metaData,
  design = ~Sample_ID)
dds

# Run DESeq function
dds <- DESeq(dds)

# Extract the results
res <- results(dds)
mcols(res, use.names=TRUE)

#Summary of results
summary(res)

# How many adjusted p-values were under 0.1? 
sum(res$padj < 0.1, na.rm=TRUE)
# I got 4
# How many adjusted p-values were under 0.05? 
sum(res$padj < 0.05, na.rm=TRUE)
#  I got 4

# Shrinkage of effect size for visualization and ranking 
resLFC <- lfcShrink(dds, coef="Sample_ID_2_vs_1", type="apeglm")

# MA plot 
plotMA( res, ylim = c(-10, 10) )
# MA plot for shrunken log2 fold change
plotMA(resLFC, ylim = c(-10, 10) )

# Volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))

# Add colored point: green if log2FC>1 and padj<0.1)
with(subset(res, padj<.1 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

#Rlog transforms the data
rld <- rlog( dds)
head( assay(rld) )

#PCA comparing the two samples
ramp <- 1:3/3
cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print( plotPCA( rld, intgroup = c( "Sample_ID")) )

#Export all with adjusted p-vaue under 0.05 to csv file
resSig <- subset(res, padj < 0.05)
write.csv(as.data.frame(resSig), 
          file="aril_aril.csv")



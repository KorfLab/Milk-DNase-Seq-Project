library(DESeq2)


# Following steps listed at:
# http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/


# 1. Make list of individual count files
directory <- "HTSeq_count_files/"
sampleFiles <- grep("txt",list.files(directory),value=TRUE)
#sampleFiles


# 2. make array of sample conditions using repetition operator
# these will pair up with corresponding file names in sampleFiles
sampleCondition <- c(rep("Lactation", 16), rep("Virgin", 14))
#sampleCondition


# 3. Make a table of sample names, file names, and conditions
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
#sampleTable


# 4. Now use DESeqDataSetFromHTSeqCount function to load all of the individual files
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
#ddsHTSeq

# 5. I think this step is setting the control condition (Virgin) in relation to the 
# test condition (Lactation). I can see that in the resulting results output, the 
# direction of log fold changes are reversed, so positive changes reflect increases 
# in lactation state
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels=c("Virgin","Lactation"))


# 6. Main DESeq2 functions
dds<-DESeq(ddsHTSeq)
res<-results(dds)
head(res)

# 7. sort results by padj column (adjusted P value)
res<-res[order(res$padj),]

head(res)
summary(res)

# 8. Main plot
plotMA(dds,ylim=c(-5,5), main="DESeq2")

# how to make a subset based on target P value 
# resSig <- subset(res, padj < 0.01)

mcols(res)$description
mcols(res,use.names=TRUE)
write.csv(as.data.frame(res),file="deseq2_output_virgin_vs_lactation.csv")


library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:15]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(1,1))
dev.copy(png,"DESeq2_heatmap1.png")
dev.off()
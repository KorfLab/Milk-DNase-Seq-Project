# run3, r = 50 (default)
data = read.table("run3_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,800), breaks=bins, col="red", main = "r=50", xlab="Insert size between mapped read pairs")
abline(v = 201.8, col="red", lwd=2)


# run13, r = 100
data = read.table("run13_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,800), breaks=bins, col="red", main = "r=100", xlab="Insert size between mapped read pairs")
abline(v = 204.2, col="red", lwd=2)

# run12, r = 200
data = read.table("run12_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,800), breaks=bins, col="red", main = "r=200", xlab="Insert size between mapped read pairs")
abline(v = 251, col="red", lwd=2)

# run14, r = 300
data = read.table("run14_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,800), breaks=bins, col="red", main = "r=300", xlab="Insert size between mapped read pairs")
abline(v = 225.1, col="red", lwd=2)

# run15, r = 400
data = read.table("run15_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,800),breaks=bins, col="red", main = "r=400", xlab="Insert size between mapped read pairs")
abline(v = 214.9, col="red", lwd=2)

# run15b, r = 400
#data = read.table("run15b_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
#summary(data$V2)
#bins=seq(-202,97350,by=5)
#hist(data$V2, add=T, xlim=c(-101,500), ylim=c(0,1500), breaks=bins, col=rgb(0, 1, 0, 0.5), main = "r=400", xlab="Insert size between mapped read pairs")
#abline(v = 222.2)

# run16, r = 400
data = read.table("run16_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
range(data$V2)
bins=seq(-101,116515,5)
hist(data$V2, add=T, col=rgb(0, 1, 0, 0.5), breaks=bins, xlim=c(-101,500), main = "r=400 Transcriptome", xlab="Insert size between mapped read pairs")
abline(v = 410, col="green", lwd=2)


# run15, r = 400
data = read.table("run15_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=10)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "r=400, stdev=20", xlab="Insert size between mapped read pairs")
abline(v = 214.9, col="red", lwd=2)


# run17, r = 400, stdev = 1382
data = read.table("run17_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=10)
hist(data$V2, add=T, xlim=c(-101,2000),breaks=bins, col=rgb(0, 1, 0, 0.5), main = "r=400, stdev=1382", xlab="Insert size between mapped read pairs")
abline(v = 529, col="green", lwd=2)



# run18, Transcriptome, r = 50, stdev = 20
data = read.table("run18_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "Insert sizes from cow RNA-seq data aligned with TopHat", xlab="Insert size between mapped read pairs")
#hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T, r=50, stdev=20", xlab="Insert size between mapped read pairs")
abline(v = 553.9, col="red", lwd=2)

# run19, Transcriptome, pre-map option, r = 50, stdev = 20
data = read.table("run19_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T -M, r=50, stdev=20", xlab="Insert size between mapped read pairs")
abline(v = 554.1, col="red", lwd=2)


# run20, Transcriptome, pre-map option, r = 400, stdev = 20
data = read.table("run20_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T -M, r=400, stdev=20", xlab="Insert size between mapped read pairs")
abline(v = 599.8, col="red", lwd=2)

# run21, Transcriptome, pre-map option, r = 800, stdev = 20
data = read.table("run21_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T -M, r=800, stdev=20", xlab="Insert size between mapped read pairs")
abline(v = 620.1, col="red", lwd=2)


# run22, Transcriptome, pre-map option, r = 800, stdev = 1000
data = read.table("run22_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T -M, r=800, stdev=1000", xlab="Insert size between mapped read pairs")
abline(v = 860, col="red", lwd=2)


# run23, Transcriptome, pre-map option, r = 800, stdev = 2000
data = read.table("run23_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T -M, r=800, stdev=2000", xlab="Insert size between mapped read pairs")
abline(v = 900.9, col="red", lwd=2)


###


# run22, Transcriptome, pre-map option, r = 800, stdev = 1000
data = read.table("run22_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,248500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "-T -M, r=800, stdev=1000", xlab="Insert size between mapped read pairs")
abline(v = 860, col="red", lwd=2)



# Sample B, run24, Transcriptome, pre-map option, r = 800, stdev = 1000
data = read.table("run24_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,160500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "Sample B, -T -M, r=800, stdev=1000", xlab="Insert size between mapped read pairs")
abline(v = 1016, col="red", lwd=2)

# Sample C, run24, Transcriptome, pre-map option, r = 800, stdev = 1000
data = read.table("run25_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,190500,by=5)
hist(data$V2, xlim=c(-101,2000),breaks=bins, col="red", main = "Sample C, -T -M, r=800, stdev=1000", xlab="Insert size between mapped read pairs")
abline(v = 699.4, col="red", lwd=2)



####


# run27, Bowtie2 against Transcriptome-index, X=500
data = read.table("run27_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,300,by=3)
hist(data$V2, xlim=c(-101, 800),breaks=bins, col="red", main = "Reads mapped to Transcriptome with Bowtie 2", xlab="Inner size between mapped read pairs")
abline(v = 33.96, col="red", lwd=2)
sd(data$V2)


# run28, Bowtie2 against Transcriptome-index, X=1000
data = read.table("run28_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
data.subset = subset(data, V2>298)
bins=seq(-101,800,by=3)
hist(data.subset$V2, add=T, xlim=c(101, 800),breaks=bins, col=rgb(0, 1, 0, 0.5), main = "Reads mapped to Transcriptome with Bowtie 2", xlab="Inner size between mapped read pairs")
abline(v = 47, col="green", lwd=2)
sd(data$V2)


Y <- rnorm(100,50,10)
range(Y)
bins=seq(22,70,4)
bins
hist(Y,breaks=bins)


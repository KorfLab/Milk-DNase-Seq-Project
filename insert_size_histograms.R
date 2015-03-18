# run3, r = 50 (default)
data = read.table("run3_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,500), breaks=bins, col="red", main = "r=50", xlab="Insert size between mapped read pairs")
abline(v = 201.8)


# run13, r = 100
data = read.table("run13_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,500), breaks=bins, col="red", main = "r=100", xlab="Insert size between mapped read pairs")
abline(v = 204.2)

# run12, r = 200
data = read.table("run12_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,500), breaks=bins, col="red", main = "r=200", xlab="Insert size between mapped read pairs")
abline(v = 251)

# run14, r = 300
data = read.table("run14_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,500), breaks=bins, col="red", main = "r=300", xlab="Insert size between mapped read pairs")
abline(v = 225.1)

# run15, r = 400
data = read.table("run15_mapping_distances.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
summary(data$V2)
bins=seq(-101,97350,by=5)
hist(data$V2, xlim=c(-101,500), breaks=bins, col="red", main = "r=500", xlab="Insert size between mapped read pairs")
abline(v = 214.9)

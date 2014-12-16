# Want to plot read counts from DNase-Seq data, binned in various intervals, across mouse chromosomes
# Have different samples, each at different DNase concentrations.

# make large bottom outer margin
#par(oma = c(7, 3, 1, 1))

# set layout of a 2x1 plot (and reduce margins of indvididual plots)
par(mfrow=c(4,1), mai=c(0.7, 0.8, 0.7, 0.2))

# read data file into table
data = read.table("Results/chr5_region_of_interest.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)


#################################
# Plot 1: 16P_5
#################################

data_subset = subset(data, data$V1 == "16P_5_mm10_MAPQ24")
plot(x=data_subset$V2,  y=data_subset$V4, type="l",  xlim=c(87000000, 89000000), ylim=c(0,100), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="16P_5")

#################################
# Plot 2: 16P_10
#################################

data_subset = subset(data, data$V1 == "16P_10_mm10_MAPQ24")
plot(x=data_subset$V2,  y=data_subset$V4, type="l",  xlim=c(87000000, 89000000), ylim=c(0,100), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="16P_10")

#################################
# Plot 3: 16P_10
#################################

data_subset = subset(data, data$V1 == "16P_20_mm10_MAPQ24")
plot(x=data_subset$V2,  y=data_subset$V4, type="l",  xlim=c(87000000, 89000000), ylim=c(0,100), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="16P_20")

#################################
# Plot 3: 16P_50
#################################

data_subset = subset(data, data$V1 == "16P_50_mm10_MAPQ24")
plot(x=data_subset$V2,  y=data_subset$V4, type="l",  xlim=c(87000000, 89000000), ylim=c(0,100), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="16P_50")


#################################
# Plot 1: 2_5_mm10_MAPQ24_chr16
#################################
plot(1, type="n", axes=F, xlab="", ylab="")
par(new=TRUE)
title(main="2_5_mm10_MAPQ24_chr16")




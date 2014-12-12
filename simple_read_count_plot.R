# Want to plot read counts from DNase-Seq data, binned in various intervals, across mouse chromosomes
# Have different samples, each at different DNase concentrations.

# make large bottom outer margin
par(oma = c(7, 3, 1, 1))

# set layout of a 2x1 plot (and reduce margins of indvididual plots)
par(mfrow=c(4,2), mai=c(0.7, 0.8, 0.7, 0.2))


#################################
# Plot 1: 2_5_mm10_MAPQ24_chr16
#################################
plot(1, type="n", axes=F, xlab="", ylab="")
par(new=TRUE)
title(main="2_5_mm10_MAPQ24_chr16")


#################################
# Plot 2: 8Mg_5_mm10_MAPQ24_chr16
#################################
data = read.table("Results/8Mg_5_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_5_mm10_MAPQ24_chr16")




#################################
# Plot 3: 2_10_mm10_MAPQ24_chr16
#################################
data = read.table("Results/2_10_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="2_10_mm10_MAPQ24_chr16")


#################################
# Plot 4: 8Mg_10_mm10_MAPQ24_chr16
#################################
data = read.table("Results/8Mg_10_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_10_mm10_MAPQ24_chr16")



#################################
# Plot 5: 2_20_mm10_MAPQ24_chr16
#################################
data = read.table("Results/2_20_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="2_20_mm10_MAPQ24_chr16")

#################################
# Plot 6: 8Mg_20_mm10_MAPQ24_chr16
#################################
data = read.table("Results/8Mg_20_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_20_mm10_MAPQ24_chr16")



#################################
# Plot 7: 2_50_mm10_MAPQ24_chr16
#################################
data = read.table("Results/2_50_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="2_50_mm10_MAPQ24_chr16")

#################################
# Plot 8: 8Mg_50_mm10_MAPQ24_chr16
#################################
data = read.table("Results/8Mg_50_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(0, 100000000), ylim=c(0,1500), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_50_mm10_MAPQ24_chr16")




# make large bottom outer margin
par(oma = c(7, 3, 1, 1))

# set layout of a 2x1 plot (and reduce margins of indvididual plots)
par(mfrow=c(2,1))

#################################
# Plot 9: 8Mg_5_mm10_MAPQ24_chr16
#################################
data = read.table("Results/8Mg_5_mm10_MAPQ24_chr16_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(10000000, 15000000), ylim=c(0,600), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_50_mm10_MAPQ24_chr16_10-15_Mbp Bin size = 1,000 bp")


data = read.table("Results/8Mg_5_mm10_MAPQ24_chr16_100_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(10000000, 15000000), ylim=c(0,300), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_50_mm10_MAPQ24_chr16_10-15_Mbp Bin size = 100 bp")

#################################
# Plot 10: 8Mg_5_mm10_MAPQ24_chr16
#################################
par(mfrow=c(2,1), mai=c(0.7, 0.8, 0.7, 0.2))
par(oma = c(7, 3, 1, 1))
data = read.table("Results/8Mg_5_mm10_MAPQ24_chr16_100_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(10000000, 10100000), ylim=c(0,15), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_50_mm10_MAPQ24_chr16 10.0-10.1_Mbp Bin size = 100 bp")
mtext("At bin sizes 1,000, 100, and 10 bp, % of bins with zero read counts = 13%, 69%, 95%", side=1, outer=T, at=0.5, cex=1.1)

data = read.table("Results/8Mg_5_mm10_MAPQ24_chr16_10_read_counts.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data$V2,  y=data$V4, type="l",  xlim=c(10000000, 10100000), ylim=c(0,15), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
title(main="8Mg_50_mm10_MAPQ24_chr16 10.0-10.1_Mbp Bin size = 10 bp")







# Add series labels
#mtext("Distance of windows from edge of block in bp", side=1, outer=T, at=0.5, cex=1.2)
#mtext("% of window occupied by feature", side=2, outer=T, at=0.5, cex=1.2)

#par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("bottom", c("Duplicated block", "Triplicated block"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n",  col = c("red", "blue"), cex=1.2, lty = c(1, 1), lwd=3)
# xpd = TRUE tells R that it is OK to plot outside the region 
# horiz = TRUE tells R that I want a horizontal legend 
# inset = c(x,y) tells R how to move the legend relative to the 'bottom' location 
# bty = 'n' means that 'no' box will be drawn around it 
# pch and col are the types and colors of points 

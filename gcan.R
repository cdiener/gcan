#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)
if( length(args)<1 ) stop(" Usage: ./gcan.R data_file") 

require(MASS)
require(ggplot2)
require(reshape2)

OD.min = 0.3
OD.max = 1.1

args = args[1]

# Functions
read.gc = function(name)
{
    ftype = tolower( strsplit(args, '\\.')[[1]][2] )
    print(ftype)
    if( ftype=="csv" ) data = read.csv(name, header=T)
    else if( ftype=="txt" ) data = read.table(name, header=1)
    else stop("Not a recognized file type. Must be txt or csv!")

	return(data)
}

normalize = function(data)
{
	ODs = as.matrix(data[,-1])
	ODs = apply(ODs, 1, function(row) row-(ODs[1,]-OD.ref)) 
	data[,-1] = t(ODs)
	return(data)
}

lm_stats = function(OD, Time)
{
	# Calculate the area under the curve
	fun = approxfun(Time, OD, method="linear")
	ABC = integrate(fun, min(Time), max(Time))$value
	
	good = which(OD>OD.min & OD<OD.max)
	if(length(good)>4)
	{
		logOD = log(OD[good])
		Time = Time[good]
		model = lm(logOD ~ Time)
		stats = c( coef(model), anova(model)[1,5], length(good), ABC)
	}
	else stats = c(0, 0, 1, length(good), ABC)
	
	return(stats)
}

lm.stats.multi = function(data)
{
	Time = data$Time
	probes = names(data[,-1])
	stats = t( apply(as.matrix(data[,-1]), 2, function(od) lm_stats(od, Time)) )
	stats = as.data.frame(stats)
	names(stats) = c("Intercept", "Growth.rate", "pval", "N_reg", "area")
	stats$Probe = sub("\\..+","", probes)

	return(stats)
}


cat("Reading growth curves...\n")
gcs = read.gc(args)
OD.ref = gcs[1,2]

#Normalize the data
cat("Normalizing...\n")
gcs = normalize(gcs)

# Get all model statistics
cat("Calculate Linear Regressions...\n")
models = lm.stats.multi(gcs)
rownames(models) <- NULL

# Plots
pdf("linear_parts.pdf", width=8, height=6)
plot(NULL, xlim=c(0,24), ylim=c(0.1,1.8), log="y", xlab="Time [h]", ylab="OD600")
 
gc = gcs[,-1]
dummy = apply(gc, 2, function(col) lines(gcs$Time, col, lwd=1, col="blue"))

abline(h=OD.min, lwd=2)
abline(h=OD.max, lwd=2)
dev.off()	

cat("Plotting...\n")

growth.plot = ggplot(models, aes(x=Probe, y=Growth.rate, group=Probe)) +
		geom_boxplot( aes(fill=Probe) ) + xlab("") + theme_bw()
        
area.plot = ggplot(models, aes(x=Probe, y=area, group=Probe)) +
		geom_boxplot( aes(fill=Probe) ) + xlab("") + theme_bw()

# to get the original growth curves in nice formatting
new = melt(gcs, id.vars="Time", variable.name="Probe", value.name="OD600")
new$Probe = sub("\\..+","", new$Probe)

curves.plot = ggplot(new, aes(x=Time, y=OD600, col=Probe, shape=Probe)) + 
		geom_pointrange(stat="summary", fun.data=mean_sdl, mult=1, size=1 ) +
		geom_line(stat="summary", fun.y=mean, size=1) + theme_bw()

cat(" - curves\n")		
ggsave(curves.plot, file="curves.pdf", width=8, height=6)

cat(" - growth rates\n")
ggsave(growth.plot, file="growth.pdf", width=6, height=4)

cat(" - area under curves\n")
ggsave(area.plot, file="area.pdf", width=6, height=4)

capture.output(print(models), file="results.txt")

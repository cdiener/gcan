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

source('libgc.R')

cat("Reading growth curves...\n")
gcs = read.gc(args)

#Normalize the data
cat("Normalizing...\n")
gcs = normalize(gcs)

# Get all model statistics
cat("Calculate Linear Regressions...\n")
models = lm.stats.multi(gcs)

# Plots
pdf("linear_parts.pdf", width=8, height=6)
plot(NULL, xlim=c(0,24), ylim=c(0.1,1.8), log="y", xlab="Time [h]", ylab="OD600")
 
gc = gcs[,-1]
dummy = apply(gc, 2, function(col) lines(gcs$Time, col, lwd=1, col="blue"))

abline(h=OD.min, lwd=2)
abline(h=OD.max, lwd=2)
dev.off()	

cat("Plotting...\n")

growth.plot = ggplot(models, aes(x=treatment, y=growth.rate)) +
		geom_boxplot( aes(fill=treatment) ) + xlab("") + theme_bw()


area.plot = ggplot(models, aes(x=treatment, y=area)) +
		geom_boxplot( aes(fill=treatment) ) + xlab("") + theme_bw()

# to get the original growth curves in nice formatting
new = melt(gcs, id.vars="Time", variable.name="Probe", value.name="OD600")
fs = names_conv(new$Probe)
new  <- cbind(new, fs)

curves.plot = ggplot(new, aes(x=Time, y=OD600, col=treatment)) + 
		geom_linerange(stat="summary", fun.data=mean_sdl, mult=1, size=1, alpha=0.5 ) +
		geom_line(stat="summary", fun.y=mean, size=1) + facet_grid(~strain) + 
		xlab("Time [h]") + theme_bw()

cat(" - curves\n")		
ggsave(curves.plot, file="curves.svg", width=6, height=4)

cat(" - growth rates\n")
ggsave(growth.plot, file="growth.svg", width=6, height=4)

cat(" - area under curves\n")
ggsave(area.plot, file="area.svg", width=6, height=4)

capture.output(print(models), file="results.txt")

#  libgc.R
#  
#  Copyright 2014 Christian Diener <christian@giiku.com>
#  
#  MIT license. See LICENSE for more information.

read.gc = function(name)
{
    ftype = tolower( strsplit(args, '\\.')[[1]][2] )
    print(ftype)
    if( ftype=="csv" ) data = read.csv(name, header=T)
    else if( ftype=="txt" ) data = read.table(name, header=1)
    else stop("Not a recognized file type. Must be txt or csv!")
	
	data$Time = data$Time*24
	return(data)
}

subsample = function(gcs, n_points)
{
    idx = seq(1:nrow(gcs), length.out=n_points)
    
    return( gcs[idx,] )
}

get_factors = function(name)
{
	name = gsub("\\..+", "", name)
	splits = strsplit(name, "_")[[1]]
	
	if (length(splits)<2) splits = c("WT", "C", NA)
	else if (length(splits)<3) splits = c( splits, NA )
	else splits = c(splits[1], splits[2], 
					paste(splits[3:length(splits)], collapse=' '))
	
	return( splits )
}

names_conv = function(name_list)
{
	splits = sapply(name_list, get_factors)
	df_info = data.frame(strain=splits[1,], treatment=splits[2,], info=splits[3,])
	
	# Resort by appearance
	for( i in 1:ncol(df_info) ) 
		df_info[,i] = factor( df_info[,i], levels=unique(df_info[,i]) )
	
	return( df_info )
}

normalize = function(data)
{
	ODs = as.matrix(data[,-1])
	OD.ref = min(ODs[1,])
	ODs = apply(ODs, 1, function(row) row-(ODs[1,]-OD.ref)) 
	data[,-1] = t(ODs)
	return(data)
}


lm_stats = function(OD, Time)
{
	# Calculate the area under the curve
	fun = approxfun(Time, OD, method="linear")
	ABC = integrate(fun, min(Time), max(Time), rel.tol=1e-3)$value
	
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
	if(max(Time)<3) Time = Time*24
	probes = names(data[,-1])
	stats = t( apply(as.matrix(data[,-1]), 2, function(od) lm_stats(od, Time)) )
	stats = as.data.frame(stats)
	names(stats) = c("intercept", "growth.rate", "pval", "n_reg", "area")
	
	info = names_conv(probes)
	stats = cbind(stats, info)
	rownames(stats) = NULL
	
	return(stats)
}

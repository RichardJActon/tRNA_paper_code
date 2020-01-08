#!/usr/bin/env Rscript

# uses sampledata as a reference for which samples to model and summarise

###################################################################
# Libs
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(optparse))

options(stringsAsFactors = FALSE)

###################################################################
# Args
option_list <- list(
	make_option(
		c("-s", "--sampledata"),
		action = "store",
		type = "character",
		help = "Path to Sample data file"
	),
	make_option(
		c("-c","--ncores"),
		action = "store",
		type = "integer",
		help = "number of threads to use",
		default = as.integer(2)
	),
	make_option(
		c("-t", "--header"),
		action = "store",
		type = "character",
		help = "header file for the data"
	), # chr_Combined_RED_rpm_bedheader.txt
	make_option(
		c("-f", "--modelFactors"),
		action = "store",
		type = "character",
		help = "sting of the form:\n 'list(c(\"sampledata$`covar`\",\"sampledata$lymphocytes\"))'\n multiple models are run depending on the number of covar vectors in the list",
		default = 'list(c(NULL),c("sampledata$neutrophils","sampledata$monocytes","sampledata$eosinophils","sampledata$lymphocytes"))'
	),
	make_option(
		c("-n", "--chunksize"),
		action = "store",
		type = "integer",
		default = as.integer(20),
		help = "number of lines to read form the file at once"
	)
)

opt <- parse_args(OptionParser(option_list = option_list))

## Arg processing

### Get list of model covars #! named lists?
modelFactors <- eval(parse(text = opt$modelFactors))

### Get data header
header <- unlist(as.character(read.table(opt$header)))

### get sample data 
# NB BGIid present in sampledata dictate which sample will be evaluated
sampledata <- read.table(opt$sampledata,header=TRUE)

####################################################################
# Functions
getDescriptives <- function(row,sampledata,factors=NULL){
	values <- as.double(as.character(row))
	normValues <- qqnorm(values,plot.it=FALSE)$x
	descriptives <- NULL
	
	descriptives$mean <- mean(values) %>% unlist()
	descriptives$median <- median(values) %>% unlist()
	descriptives$min <- min(values) %>% unlist()
	descriptives$max <- max(values) %>% unlist()
	descriptives$sd <- sd(values) %>% unlist()
	descriptives$var <- var(values) %>% unlist()
	
	descriptives$Nmean <- mean(normValues) %>% unlist()
	descriptives$Nmedian <- median(normValues) %>% unlist()
	descriptives$Nmin <- min(normValues) %>% unlist()
	descriptives$Nmax <- max(normValues) %>% unlist()
	descriptives$Nsd <- sd(normValues) %>% unlist()
	descriptives$Nvar <- var(normValues) %>% unlist()
	
	descriptives$percentZero <- (
		length(values[values == as.double(0)]) / length(values)
	) * 100 %>% unlist()
	
	# age modelling (multiple models permitted)
	for (i in 1:length(factors)) {
		factor <- c("normValues",factors[[i]])
		model <- as.formula(
			paste(
				"sampledata$ageDNAextraction ~ ", paste(factor, collapse="+")
			)
		)
		result <- lm(model)
		
		descriptives[paste0("slope",i)] <- summary(result)$coeff[2,1] %>% 
			unlist()
		descriptives[paste0("p",i)] <- summary(result)$coeff[2,4] %>% 
			unlist()
		descriptives[paste0("neglog10P",i)] <- 
			-log10(descriptives[[paste0("p",i)]]) %>% unlist()
	}
	return(descriptives %>% as.data.frame())
}

#####################################################################
# File reading loop
# NB reads file N (`opt$chunksize`) lines at a time computes descriptives
# and models for that chunk and moves onto the next
# execution of the `getDescriptives` function is done in parallel for each
# line in the chunk number of threads can be specified with `opt$ncores`

infile <- file("stdin","r")

while(length(lines <- readLines(con = infile, n = opt$chunksize)) > 0) {
	
	# parse lines to data frame
	mat <- lines %>%
		strsplit(.,"\n") %>%
		unlist(.,recursive = FALSE) %>%
		lapply(.,function(x){do.call(rbind,strsplit(x, "\t"))}) %>%
		do.call(rbind.data.frame,.)
	# name columns with header file - permits samplesheet column look-up
	colnames(mat) <- header
	# select only those columns present in the samplesheet
	mat <- mat[,c(1:3,which(colnames(mat) %in% sampledata$BGIid))]
	# apply get descriptive linewise in parallel
	desc <- mat[,4:ncol(mat)] %>% # remove chr/start/stop
		t() %>%
		as.data.frame() %>%
		mclapply(
			.,getDescriptives,
			sampledata = sampledata,
			mc.cores = opt$ncores,
			factors = modelFactors
		) %>% 
		do.call(rbind, .) %>%
		as.data.frame() %>% 
		cbind(mat[,1:3], .) # add back chr/start/stop
	write.table(
		desc, file = stdout(), quote = FALSE, sep = "\t",
		col.names = FALSE, row.names = FALSE
	)
}

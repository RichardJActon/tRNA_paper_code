poolColours <-c(
	"4.07"="#ffffcc","4.09"="#CCCC56",
	"28.07"="#a1dab4","28.23"="#3CA661",
	
	"63.26"="#225ea8","63.4"="#14488C",
	"77.22"="#41b6c4","77.96"="#2899A8"#,
	#"DNAm Age" = "red"
)
#library(ggsignif)
tRNAMethCovPlotter <- function(
	plotdata,
	tRNAs = NULL,
	uniformLims = FALSE,
	RnBeadsPairwiseComparisons = NULL
) {
	compLabDat <- NULL
	if (! is.null(RnBeadsPairwiseComparisons)) {
		compLabDat <- 
		left_join(
			RnBeadsPairwiseComparisons %>%
				extract(col = comparison,into = c("start","end"),regex = "(\\w+) vs. (\\w+) ") %>%
				dplyr::select(tRNAname=name,start,end,comb.p.val,comb.p.adj.fdr,num.sites,mean.mean.diff) %>%
				mutate(
					#comb.p.val = sprintf("p-value: %4.3g",comb.p.val),#mean.mean.diff
					comb.p.val = sprintf("p-value: %4.2e, diff: %3.2f",comb.p.val,(mean.mean.diff*100)),#mean.mean.diff
					comb.p.adj.fdr = sprintf("p-value (BH): %4.2e",comb.p.adj.fdr)
				),
			
			plotdata %>% 
				filter(type == "meth") %>%
				group_by(tRNAname) %>%
				dplyr::summarise(y = max(values) * 1.05),
			by = "tRNAname"
		) %>% 
			na.omit() %>%
			mutate(type = "meth") %>% 
			nest(-tRNAname) %>%
			mutate(
				y = purrr::map(
					data,
					~(function(q){
						#q + (0.4 + seq_along(q))
						#q + ((sin(.$y) + 0.15) + seq_along(q))
						q + ( (0.1*.$y) * seq_along(q) )
					})(.$y)
				)
			) %>%
			unnest() %>%
			dplyr::select(-y1)
	}
	
	# plotData columns
	cols <- c(
		"type", "tRNAname", "values", "tRNAge", "ageBand", "age",
		"tStrand", "chr", "tStart", "tEnd", "pseudo"
	)
	if (all(cols %in% colnames(plotdata))) {
		if (! all(c("meth","coverage") %in% unique(plotdata$type))) {
			stop("type must contain 'meth' & 'coverage'")
		}
	} else {
		missingCols <- paste0(
			cols[!(cols %in% colnames(plotdata))], collapse = ", "
		)
		stop(paste0("Some required Columns are missing: \n", missingCols))
	}
	
	if (is.null(tRNAs)) {
		tRNAs <- unique(plotdata$tRNAname)
	}
	
	if (uniformLims == TRUE) {
		methLimsData <- data.frame( ##!! NB conflict with ggsignif bracket max values - need to increase to accomodate
			type="meth",
			values=range(
				plotdata %>%
					filter(tRNAname %in% tRNAs) %>%
					filter(type=="meth") %>%
					pull(values)
			)
		)
		covLimsData <- data.frame(
			type="coverage",
			values=range(
				plotdata %>%
					filter(tRNAname %in% tRNAs) %>%
					filter(type=="coverage") %>%
					pull(values)
			)
		)
	}

	plots <- plotdata %>%
		filter(tRNAname %in% tRNAs) %>%
		mutate(age = as.double(as.character(age))) %>%
		arrange(age) %>%
		mutate(age=as.character(round(age,2))) %>%
		mutate(age=factor(age,levels=unique(age))) %>%
		mutate(
			ageBand=factor(ageBand,levels=unique(ageBand),ordered = TRUE)
		) %>%
		group_by(tRNAname) %>%
		do(plot = {
			tRNA <- .$tRNAname
			p <- ggplot(., aes(ageBand,values)) + 
				# geom_boxplot(
				# 	data = as.data.frame(.) %>%
				# 		filter(type =="coverage"),
				# 	position = "dodge",
				# 	show.legend = FALSE,
				# 	#position = position_dodge(width = 0.1),
				# 	aes(ageBand, values,
				# 		fill = age
				# 	)
				# ) +
				geom_col(
					data = as.data.frame(.) %>%
						filter(type =="coverage"),
					position = "dodge",
					show.legend = FALSE,
					#position = position_dodge(width = 0.1),
					aes(ageBand, values,
						fill = age
					)
				) +
				geom_errorbar(
					data = as.data.frame(.) %>%
						filter(type =="coverage"),
					position = position_dodge(0.9),#"dodge",
					show.legend = FALSE,
					#position = position_dodge(width = 0.1),
					aes(
						ageBand,
						ymax = values + sd,
						ymin = values - sd,
						fill = age
					),
					width = 0.2
				) +
				geom_violin(
					data = as.data.frame(.) %>% 
						filter(type == "meth"),
					position = "dodge",
					aes(ageBand, values, fill = age),
					draw_quantiles = c(0.25, 0.5, 0.75),
					alpha = 0.8
				) +
				geom_jitter(
					data = as.data.frame(.) %>% 
						filter(type == "meth"),
					position = position_jitterdodge(),
					show.legend = FALSE,
					aes(ageBand, values,fill = age),
					size = 0.6
				) +
				scale_fill_manual(
					values = rep(
						each = 2,
						poolColours
						#c("#ffffcc","#a1dab4","#41b6c4","#225ea8") # Names?
					)
				) +
				scale_colour_manual(
					values = c("TRUE" = "red","FALSE" = "darkgrey"),
					limits = c("TRUE", "FALSE"),
					drop = FALSE
				) +
				#scale_shape_manual(values = 21:24) + 
				stat_summary(
					data = as.data.frame(.) %>% 
						filter(type == "meth"),
					aes(colour = tRNAge, group = 1),#y=value,
					fun.y = mean,
					geom = "line",
					size = 1
				) +
				ylim(0,NA) + 
				labs(
					title = paste0(
						na.omit(.$tRNAname)," (",
						na.omit(.$tStrand),") ",
						tools::toTitleCase(na.omit(.$chr))," : ",
						format(
							na.omit(.$tStart),big.mark=",",trim = TRUE,
							scientific = FALSE
						)
						," - ",
						format(
							na.omit(.$tEnd),big.mark=",",trim = TRUE,
							scientific = FALSE
						)
						," [",
						na.omit(.$pseudo),"]"
					),
					#subtitle = "Targeted BS-seq, 4 replicates from 2 pools of 20-25 people at 4 timepoints",
					#y = "					meth / % methylated		total coverage / N reads",
					y = "					meth / % methylated			coverage / N reads",
					x = "Age Group"
				) +
				guides(fill = guide_legend(title = "Pool Mean\nAge / yrs")) + 
				facet_wrap(~type, scales = "free_y", nrow = 2)
			
			if (uniformLims == TRUE) {
				p <- p + 
					geom_blank(
						data = methLimsData, aes(y = values),
						inherit.aes = FALSE
					) + 
					geom_blank(
						data = covLimsData, aes(y = values),
						inherit.aes = FALSE
					)
			}
			
			if (! is.null(compLabDat)) {
				p <- p + 
					geom_signif(
						inherit.aes = FALSE,
						data = compLabDat %>% filter(tRNAname %in% tRNA),
						aes(
							xmin = start,
							xmax = end,
							annotations = comb.p.val,
							y_position = y
						),
						#step_increase = c(5,10,15,20,20,30),
						textsize = 3,
						vjust = -0.2,
						manual=TRUE
					) + 
					labs(caption = "p-values from RnBeads - limma")
			}
			
			p
		})
	return(plots)
}


plotPrinter <- function(
	plotsDF,
	coverageBoxProp = 0.3,
	printNoSave = TRUE,
	outDir = "../graphics/",
	filePrefix,
	units = "in",
	width = 12,
	height = 6.75,
	res = 384
) {
	plotCol <- ncol(plotsDF)
	Nvar <- plotCol - 1
	
	lapply(seq_along(plotsDF$plot), function(i){
		gt <- plotsDF[[i,plotCol]] %>% ggplot_build() %>% ggplot_gtable()
		gt$heights[8] = coverageBoxProp * gt$heights[8]
		#gt$heights[18] = coverageBoxProp * gt$heights[18]
		if (printNoSave){
			grid.newpage()
			grid.draw(gt)
		} else {
			dir.create(
				paste0(outDir,filePrefix),
				showWarnings = FALSE,
				recursive = TRUE
			)
			png(
				filename = paste0(
					outDir,
					filePrefix,"/",
					filePrefix,"_",
					paste0(plotsDF[i,1:Nvar],collapse = "_"),
					".png"
				),
				units = units,
				width = width,
				height = height,
				res = res
			)
			grid.draw(gt)
			dev.off()
		}
	})
}



######################################################################################
######################################################################################
######################################################################################
# tRNAMethCovPlotter <- function(
# 	plotdata,
# 	tRNAs = NULL,
# 	uniformLims = FALSE,
# 	#comparisons = FALSE,
# ) {
# 	#data 
# 	# cols:
# 	# tRNAname
# 	# type (== "coverage" "meth")
# 	# values (double)
# 	# tRNAge (logical)
# 	# ageBand
# 	# age
# 	cols <- c(
# 		"type", "tRNAname", "values", "tRNAge", "ageBand", "age",
# 		"tStrand", "chr", "tStart", "tEnd", "pseudo"
# 	)
# 	if (all(cols %in% colnames(plotdata))) {
# 		if (! all(c("meth","coverage") %in% unique(plotdata$type))) {
# 			stop("type must contain 'meth' & 'coverage'")
# 		}
# 	} else {
# 		missingCols <- paste0(
# 			cols[!(cols %in% colnames(plotdata))], collapse = ", "
# 		)
# 		stop(paste0("Some required Columns are missing: \n", missingCols))
# 	}
# 	
# 	if (is.null(tRNAs)) {
# 		tRNAs <- unique(plotdata$tRNAname)
# 	}
# 	#tRNAs <- c("tRNA-Ile-AAT-1-1","tRNA-Ile-AAT-5-4","tRNA-Ile-AAT-10-1","tRNA-Ile-AAT-5-5")
# 	#tRNA-Ile-AAT-5-5 %>%#"tRNA-iMet-CAT-1-4",#tRNA-Ile-AAT-10-1
# 	#filter(!tRNAname %in% c("tRNA-Leu-AAG-2-3","tRNA-Ser-AGA-2-5")) %>% # problem plots?
# 	#filter(tRNAname %in% c("tRNA-Leu-AAG-2-3","tRNA-Ser-AGA-2-5")) %>% # problem plots?
# 	
# 	if (uniformLims == TRUE) {
# 		methLimsData <- data.frame(
# 			type="meth",
# 			values=range(
# 				plotdata %>%
# 					filter(tRNAname %in% tRNAs) %>%
# 					filter(type=="meth") %>%
# 					pull(values)
# 			)
# 		)
# 		covLimsData <- data.frame(
# 			type="coverage",
# 			values=range(
# 				plotdata %>%
# 					filter(tRNAname %in% tRNAs) %>%
# 					filter(type=="coverage") %>%
# 					pull(values)
# 			)
# 		)
# 	}
# 	
# 	plots <- plotdata %>%
# 		filter(tRNAname %in% tRNAs) %>%
# 		mutate(age = as.double(as.character(age))) %>%
# 		arrange(age) %>%
# 		mutate(age=as.character(round(age,2))) %>%
# 		mutate(age=factor(age,levels=unique(age))) %>%
# 		mutate(
# 			ageBand=factor(ageBand,levels=unique(ageBand),ordered = TRUE)
# 		) %>%
# 		group_by(tRNAname) %>%
# 		do(plot = {
# 			p <- ggplot(., aes(ageBand,values)) + 
# 				#ggplot(data=.) + 
# 				#geom_point(alpha=0.8) + #aes(fill=pool)
# 				geom_boxplot(
# 					data = as.data.frame(.) %>% 
# 						filter(type =="coverage"),
# 					position = "dodge",
# 					show.legend = FALSE,
# 					#position = position_dodge(width = 0.1),
# 					aes(ageBand, values,
# 						fill = age#,
# 						#group=ageBand
# 					)
# 				) + 
# 				geom_violin(
# 					data = as.data.frame(.) %>% 
# 						filter(type == "meth"),
# 					position = "dodge",
# 					#aes(x=reorder(age,sort(as.numeric(as.character(age)))),values,fill=ageBand),
# 					aes(ageBand, values, fill = age),
# 					draw_quantiles = c(0.25, 0.5, 0.75),
# 					alpha = 0.8
# 				) +
# 				geom_jitter(
# 					data = as.data.frame(.) %>% 
# 						filter(type == "meth"),
# 					position = position_jitterdodge(),
# 					show.legend = FALSE,
# 					#aes(x=reorder(age,sort(as.numeric(as.character(age)))),values),
# 					aes(ageBand, values,fill = age),
# 					#colour="black",
# 					size = 0.6
# 				) +
# 				#scale_fill_brewer(palette = 3,type="qual") +
# 				#scale_fill_brewer(palette = 1) +
# 				# scale_fill_manual(values = c("4"="#EDC947","28"="#49A644",
# 				# 							 "63"="#368FB5","77"="#4E2FB5")) +
# 				scale_fill_manual(
# 					values = rep(
# 						each = 2,
# 						#c("gold2","green3","deepskyblue2","navy")
# 						c("#ffffcc","#a1dab4","#41b6c4","#225ea8")
# 						#c("#EDC947","#49A644","#368FB5","#4E2FB5")
# 					)
# 				) +
# 				scale_colour_manual(
# 					values = c("TRUE" = "red","FALSE" = "darkgrey"),
# 					limits = c("TRUE", "FALSE"),
# 					drop = FALSE
# 				) +
# 				#scale_shape_manual(values = 21:24) + 
# 				stat_summary(
# 					data = as.data.frame(.) %>% 
# 						filter(type == "meth"),
# 					aes(colour = tRNAge, group = 1),#y=value,
# 					fun.y = mean,
# 					geom = "line",
# 					size = 1
# 				) +
# 				labs(
# 					title = paste0(
# 						na.omit(.$tRNAname)," (",
# 						na.omit(.$tStrand),") ",
# 						tools::toTitleCase(na.omit(.$chr))," : ",
# 						format(
# 							na.omit(.$tStart),big.mark=",",trim = TRUE,
# 							scientific = FALSE
# 						)
# 						," - ",
# 						format(
# 							na.omit(.$tEnd),big.mark=",",trim = TRUE,
# 							scientific = FALSE
# 						)
# 						," [",
# 						na.omit(.$pseudo),"]"
# 					),#.$tRNAname,
# 					#title="Two Isodecoders One Previously Identified as Hypermethylating With Age",
# 					subtitle = "Targeted BS-seq, 4 replicates from 2 pools of 25 people at 4 timepoints",
# 					y = "					meth / % methylated			coverage / N reads",
# 					x = "Age Group"
# 				) +
# 				guides(fill = guide_legend(title = "Pool Mean\nAge / yrs")) + 
# 				facet_wrap(~type, scales = "free_y", nrow = 2)
# 			#facet_grid(data~tRNAname,scales = "free")
# 			
# 			if (uniformLims == TRUE) {
# 				p <- p + 
# 					geom_blank(
# 						data = methLimsData, aes(y = values),
# 						inherit.aes = FALSE
# 					) + 
# 					geom_blank(
# 						data = covLimsData, aes(y = values),
# 						inherit.aes = FALSE
# 					)
# 			}
# 			p
# 		})
# 	return(plots)
# }

## debugging ggsignif data
# 
# left_join(
# 	diffMethRestRNA %>%
# 		extract(col = comparison,into = c("start","end"),regex = "(\\w+) vs. (\\w+) ") %>%
# 		dplyr::select(tRNAname=name,start,end,comb.p.val,comb.p.adj.fdr,num.sites),
# 	
# 	plotDataBytRNA %>%
# 		filter(type == "meth") %>%
# 		group_by(tRNAname) %>%
# 		dplyr::summarise(
# 			y = max(values) * 1.05#,
# 			#max(values) - min(values)
# 		),
# 	by = "tRNAname"
# ) %>% na.omit() %>%
# 	nest(-tRNAname) %>%
# 	mutate(
# 		y = purrr::map(
# 			data,
# 			~(function(q){q + ((sin(.$y) + 0.15) + seq_along(q))})(.$y)
# 		)
# 	) %>%
# 	unnest() %>%
# 	dplyr::select(-y1)

#group_by(tRNAname) %>%
#dplyr::summarise(min=min(y),max=max(y))

######################################################################################################################
# tRNAMethCovSitesPlotter <- function(
# 	plotdata,
# 	tRNAs = NULL,
# 	uniformLims = FALSE,
# 	RnBeadsPairwiseComparisons = NULL
# ) {
# 	# compLabDat <- NULL
# 	# if (! is.null(RnBeadsPairwiseComparisons)) {
# 	# 	compLabDat <- 
# 	# 		left_join(
# 	# 			RnBeadsPairwiseComparisons %>%
# 	# 				extract(col = comparison,into = c("start","end"),regex = "(\\w+) vs. (\\w+) ") %>%
# 	# 				dplyr::select(tRNAname=name,start,end,comb.p.val,comb.p.adj.fdr,num.sites,mean.mean.diff) %>%
# 	# 				mutate(
# 	# 					#comb.p.val = sprintf("p-value: %4.3g",comb.p.val),#mean.mean.diff
# 	# 					comb.p.val = sprintf("p-value: %4.3g, diff: %3.2f",comb.p.val,(mean.mean.diff*100)),#mean.mean.diff
# 	# 					comb.p.adj.fdr = sprintf("p-value (BH): %4.3g",comb.p.adj.fdr)
# 	# 				),
# 	# 			
# 	# 			plotdata %>% 
# 	# 				filter(type == "meth") %>%
# 	# 				group_by(tRNAname) %>%
# 	# 				dplyr::summarise(y = max(values) * 1.05),
# 	# 			by = "tRNAname"
# 	# 		) %>% 
# 	# 		na.omit() %>%
# 	# 		mutate(type = "meth") %>% 
# 	# 		nest(-tRNAname) %>%
# 	# 		mutate(
# 	# 			y = purrr::map(
# 	# 				data,
# 	# 				~(function(q){
# 	# 					#q + (0.4 + seq_along(q))
# 	# 					#q + ((sin(.$y) + 0.15) + seq_along(q))
# 	# 					q + ( (0.1*.$y) * seq_along(q) )
# 	# 				})(.$y)
# 	# 			)
# 	# 		) %>%
# 	# 		unnest() %>%
# 	# 		dplyr::select(-y1)
# 	# }
# 	
# 	# plotData columns
# 	cols <- c(
# 		"type", "tRNAname", "values", "tRNAge", "ageBand", "age",
# 		"tStrand", "chr", "tStart", "tEnd", "pseudo"
# 	)
# 	if (all(cols %in% colnames(plotdata))) {
# 		if (! all(c("meth","coverage") %in% unique(plotdata$type))) {
# 			stop("type must contain 'meth' & 'coverage'")
# 		}
# 	} else {
# 		missingCols <- paste0(
# 			cols[!(cols %in% colnames(plotdata))], collapse = ", "
# 		)
# 		stop(paste0("Some required Columns are missing: \n", missingCols))
# 	}
# 	
# 	if (is.null(tRNAs)) {
# 		tRNAs <- unique(plotdata$tRNAname)
# 	}
# 	
# 	if (uniformLims == TRUE) {
# 		methLimsData <- data.frame( ##!! NB conflict with ggsignif bracket max values - need to increase to accomodate
# 			type = "meth",
# 			values = range(
# 				plotdata %>%
# 					filter(tRNAname %in% tRNAs) %>%
# 					filter(type == "meth") %>%
# 					pull(values)
# 			)
# 		)
# 		covLimsData <- data.frame(
# 			type = "coverage",
# 			values = range(
# 				plotdata %>%
# 					filter(tRNAname %in% tRNAs) %>%
# 					filter(type == "coverage") %>%
# 					pull(values)
# 			)
# 		)
# 	}
# 	
# 	plots <- plotdata %>%
# 		filter(tRNAname %in% tRNAs) %>%
# 		mutate(age = as.double(as.character(age))) %>%
# 		arrange(age) %>%
# 		mutate(age=as.character(round(age,2))) %>%
# 		mutate(age=factor(age,levels=unique(age))) %>%
# 		mutate(
# 			ageBand=factor(ageBand,levels=unique(ageBand),ordered = TRUE)
# 		) %>%
# 		group_by(tRNAname) %>%
# 		do(plot = {
# 			tRNA <- .$tRNAname
# 			p <- ggplot(., aes(ageBand,values)) + 
# 				geom_boxplot(
# 					data = as.data.frame(.) %>% 
# 						filter(type =="coverage"),
# 					position = "dodge",
# 					show.legend = FALSE,
# 					#position = position_dodge(width = 0.1),
# 					aes(ageBand, values,
# 						fill = age
# 					)
# 				) + 
# 				# geom_violin(
# 				# 	data = as.data.frame(.) %>% 
# 				# 		filter(type == "meth"),
# 				# 	position = "dodge",
# 				# 	aes(ageBand, values, fill = age),
# 				# 	draw_quantiles = c(0.25, 0.5, 0.75),
# 				# 	alpha = 0.8
# 				# ) +
# 				#geom_jitter(
# 				geom_point(
# 					data = as.data.frame(.) %>% 
# 						filter(type == "meth"),
# 					#position = position_jitterdodge(),
# 					show.legend = FALSE,
# 					aes(ageBand, values,fill = age),
# 					size = 0.6
# 				) +
# 				scale_fill_manual(
# 					values = rep(
# 						each = 2,
# 						c("#ffffcc","#a1dab4","#41b6c4","#225ea8") # Names?
# 					)
# 				) +
# 				scale_colour_manual(
# 					values = c("TRUE" = "red","FALSE" = "darkgrey"),
# 					limits = c("TRUE", "FALSE"),
# 					drop = FALSE
# 				) +
# 				#scale_shape_manual(values = 21:24) + 
# 				# stat_summary(
# 				# 	data = as.data.frame(.) %>% 
# 				# 		filter(type == "meth"),
# 				# 	aes(colour = tRNAge, group = 1),#y=value,
# 				# 	fun.y = mean,
# 				# 	geom = "line",
# 				# 	size = 1
# 				# ) +
# 				ylim(0,NA) + 
# 				labs(
# 					title = paste0(
# 						na.omit(.$tRNAname)," (",
# 						na.omit(.$tStrand),") ",
# 						tools::toTitleCase(na.omit(.$chr))," : ",
# 						format(
# 							na.omit(.$tStart),big.mark=",",trim = TRUE,
# 							scientific = FALSE
# 						)
# 						," - ",
# 						format(
# 							na.omit(.$tEnd),big.mark=",",trim = TRUE,
# 							scientific = FALSE
# 						)
# 						," [",
# 						na.omit(.$pseudo),"]"
# 					),
# 					subtitle = "Targeted BS-seq, 4 replicates from 2 pools of 25 people at 4 timepoints",
# 					y = "					meth / % methylated			coverage / N reads",
# 					x = "Age Group"
# 				) +
# 				guides(fill = guide_legend(title = "Pool Mean\nAge / yrs")) + 
# 				facet_wrap(~type, scales = "free_y", nrow = 2)
# 			
# 			if (uniformLims == TRUE) {
# 				p <- p + 
# 					geom_blank(
# 						data = methLimsData, aes(y = values),
# 						inherit.aes = FALSE
# 					) + 
# 					geom_blank(
# 						data = covLimsData, aes(y = values),
# 						inherit.aes = FALSE
# 					)
# 			}
# 			
# 			# if (! is.null(compLabDat)) {
# 			# 	p <- p + 
# 			# 		geom_signif(
# 			# 			inherit.aes = FALSE,
# 			# 			data = compLabDat %>% filter(tRNAname %in% tRNA),
# 			# 			aes(
# 			# 				xmin = start,
# 			# 				xmax = end,
# 			# 				annotations = comb.p.val,
# 			# 				y_position = y
# 			# 			),
# 			# 			#step_increase = c(5,10,15,20,20,30),
# 			# 			textsize = 3,
# 			# 			vjust = -0.2,
# 			# 			manual=TRUE
# 			# 		) + 
# 			# 		labs(caption = "p-values from RnBeads - limma")
# 			# }
# 			
# 			p
# 		})
# 	return(plots)
# }

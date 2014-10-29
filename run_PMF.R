###
# Functions to run PMF using the multilinear engine
#


# folder where running these functions from must have
# 1. iniParams.txt initial parameters
# 2. PMF_bs2.ini initial file for windows
# 3. me2wopt.exe ME2 executable
# and must have write permissions to write .dat files

#requires adjust function from handles package
library(handles)

#####
# runPMF takes in complete datafile and uncertainty file
# and runs PMF using ME2 to get source apportioned results.
#
# datafile is (path) + (file name as csv)
# 	(ndays X ncons (with PM, without date))
# ndays is number of days of observation
# ncons is number of constituents (including PM if used)
# nfact is number of sources
# adjust is censoring adjustment method
# mdlfile is csv file of mdls corresponding to data
# uncfile is csv file of uncertainties corresponding to data
# total indicates whether first column of data is a total
# seeds are optional seeds to pass to ME2
runPMF <- function(datafile, ndays, ncons, nfact, adjust = NULL,
	mdlfile = NULL, uncfile = NULL, total = F, seeds = NULL, pref = "testout") {

    #read in data
	dat <- read.csv(datafile)



	#add uncertainty
	if(!is.null(uncfile)) {
		unc <- read.csv(uncfile)
	}else{
		unc <- matrix(rep(1, nrow(dat) * ncol(dat)), nrow = nrow(dat))
	}
	#make total a weak variable according to PMF rules
	if(total) {
		unc[, 1] <- unc[, 1] * 3
	}


	#adjust below the MDL if needed (substitute, snrat)
	if(!is.null(adjust)) {
		mdl <- read.csv(mdlfile)

		if(class(dat[, 1]) != "Date") {
			dates <- as.Date(rep(1, nrow(dat)), origin = "1970-01-01")
			dat <- data.frame(dates, dat)
		}

		#adjust data removing dates
		adj <- adjust(dat = dat, mdl = mdl, adjust = adjust, unc = unc)
		dat <- adj$dat[, -1]

		if(adjust == "snrat") {
			unc <- adj$unc
		}
	}




	#write concatenated data
	write.table(dat, file = "dat.txt", row.names = F, col.names = F)
	write.table(unc, file = "dat.txt", row.names = F, col.names = F,
		append = T)

	#read current initial parameters using PMF default
	iniParams <- readLines("iniParams.txt")

	#if want to specify seeds directly
	if(!is.null(seeds)) {
		iniParams[6] <- paste(seeds, collapse = " ")
	}

	update number of days, constituents, and sources needed
	iniParams[7] <- paste(c(ndays, ncons, nfact, 0, 0, -12), collapse = " ")
	writeLines(iniParams, "iniParams2.txt")


	#call ME2 (PMF_bs2 references iniParams.txt for params)
	system("me2wopt.exe PMF_bs2", show.output.on.console = F, wait = T)

	#interpret results of PMF/ME2
	solns1 <- readME2(ndays, ncons, nfact, nruns, total, pref)

	return(solns1)
}









###
# Function to apply PMF/ME2
#
# ndays is number of days
# ncons is number of constituents (including PM25)
# nfact is number of sources
# nruns is number of PMF runs to complete
# pref is prefix for data
readME2 <- function(ndays, ncons, nfact, nruns, total = T, pref = "testout"){

	#read in PMF results
	all <- read.table(paste(pref, ".dat", sep= ""),	skip = 2, sep = "",
		strip.white = T)

  	#check output has correct number of lines
  	neach <- ncons + ndays + 1
  	if((nruns * neach) != nrow(all)) {
		print(neach)
  		print(nrow(all))

	}



	#separate by PMF run, name columns as factors
	all <- data.frame(rep(seq(1 : nruns), each = neach), all)
	colnames(all) <- c("nrun", paste("factor", seq(1 : nfact), sep = ""))

	#find best PMF run
	diag <- read.table(paste(pref, ".txt", sep = ""), skip = 2, sep = "",
		fill = T, colClasses = c("character"))
	tasks <- which(diag[, 1] == "task") + 1
	wrun <- which.min(as.numeric(diag[tasks, 2]))

	#Select best PMF run
	all <- all[which(all[, 1] == wrun), ]



	#Get concentrations
	contr <- all[1 : ndays, -1]

	#Get profiles: remove spacer row
	prof <- all[(ndays + 1) : (neach - 1), -1]
	profo <- prof


	#remove PM25 and rescale profile
	prof <- prof[-1, ]
	prof <- sweep(prof, 2, colSums(prof), "/")


	#if use PM2.5, need to rescale concentrations
	if(total == T){
		contr <- sweep(contr, 2, as.numeric(profo[1, ]), "*")
	}


	list(contr = contr, prof = prof)
}

###
# Functions to run PMF using the multilinear engine
#




###
# Function to run PMF
# 
# dats is list with first element data (nDays X nConstituents)
# 	and second element uncertainties (same dimension as data)
#
# ndays is number of days
# ncons is number of chemical constituents (including total)
outerPMF <- function(dats, nameprev = "nothing.dat", ndays, 
	ncons, nfact, total = FALSE, 
	robust = 1, posoutdist = 4,
	negoutdist = 4, precmode = 20,
    nruns = 20, numoldsol = 0, bsinitfact = 1,
    bsmode = 11, simu = 0, contrun = 0,
    dobspull = 0, pullc1 = 1.5, readbscnts = 0,
    samplevari = 1, alowlim = -0.2, normc1 = 0.0,
    acbmodel=0, seed1 = 95, seed2 = 75, seed3 = 75,
    seed4 = 69, seed5 = 72, c1 = 0, c3 = 0, em = -12,
    numpf = 0, numynpf = 0, maxpfdim = 20, modelc1 = 1.2,
    modelc3 = 0.35, pfpullc1 = 1.1, pfsmooc1 = 0.5,
    doresort = 1, dofpeak = 1){
  
  	#create input txt file and run PMF/ME2
	createiniparams(dats = dats, nameprev = nameprev, 
		nameout = "testout", iv = "261107",
        n1 = ndays, n2 = ncons, np = nfact, 
        robust, posoutdist, negoutdist, precmode,
        nruns, numoldsol, bsinitfact, bsmode, simu,
        contrun, dobspull, pullc1, readbscnts, samplevari,
        alowlim, normc1, acbmodel, seed1, seed2, seed3,
        seed4, seed5, c1, c3, em,
        numpf, numynpf, maxpfdim, modelc1, modelc3,
        pfpullc1, pfsmooc1, doresort, dofpeak)
  
	#interpret results of PMF/ME2
	solns1 <- readME2("testout", ndays, nfact, ncons, nruns)
  
  	#create output
 	solns <- list(contr = solns1$contr, prof = solns1$prof)
  	
  	#if use PM2.5, need to rescale concentrations
  	if(total == T){
  		contr <- sweep(solns1$contr, 2, 
  			as.numeric(solns1$profo[1, ]), "*")
		solns <- list(contr, solns1$prof)
	}  

	solns
  
}










###
# Create input text file and run PMF
# 
# dats is list with first element data (nDays X nConstituents)
# 	and second element uncertainties (same dimension as data)
# nameprev is previous run (generally ignore) 
# nameout is name of output file
# iv is key (generally ignore)
# n1 is number of days
# n2 is number of constituents (plus PM25)
# np is number of sources
# 
# There are many other possible inputs.  
# See ME2 documentation for details
createiniparams <- function(dats,nameprev="nothing.dat",
	nameout, iv, n1, n2, np,
	robust = 1, posoutdist = 4,
	negoutdist = 4, precmode = 20,
    numtasks = 20, numoldsol = 0,
    bsinitfact = 1, bsmode = 11, simu = 0,
    contrun = 0, dobspull = 0, pullc1 = 1.5,
    readbscnts = 0, samplevari = 1,
    alowlim = -0.2, normc1 = 0.0, acbmodel = 0,
    seed1 = 24, seed2 = 47, seed3 = 54,
    seed4 = 56, seed5 = 76, c1 = 0, c3 = 0, em = -12,
    numpf = 0, numynpf = 0, maxpfdim = 20,
    modelc1 = 1.2, modelc3 = 0.35,
    pfpullc1 = 1.1, pfsmooc1 = 0.5, doresort = 1,
    dofpeak = 1){
  
	#create lines for text file
	l1 <- paste(iv)
	l2 <- paste(robust, posoutdist, negoutdist, 
		precmode, numtasks, numoldsol, collapse = "\t")
	l3 <- paste(bsinitfact, bsmode, simu, collapse = "\t")
	l4 <- paste(contrun, dobspull, pullc1, readbscnts,
		samplevari, collapse = "\t")
	l5 <- paste(alowlim, normc1, acbmodel, collapse = "\t")
	l6 <- paste(seed1, seed2, seed3, seed4, 
		seed5, collapse = "\t")
	l7 <- paste(n1, n2, np, c1, c3, em, collapse = "\t")
	l8 <- paste(numpf, numynpf, maxpfdim, modelc1, modelc3, 
		pfpullc1, pfsmooc1, collapse = "\t")
	l10 <- paste(doresort, dofpeak)
  
 
  	#create datafile
  	write.table(dats[[1]], file = "dat.txt",
  		row.names = F, col.names = F)
  	#add uncertainty
  	write.table(dats[[2]], file = "dat.txt",
  		row.names = F, col.names = F, append = T)

  
	nameprev2 <- paste("'", nameprev, "'", sep = "")
	nameout2 <- paste("'", nameout, "'", sep = "")
  	l9 <- paste("'dat.txt'", nameprev2, nameout2, 
  		collapse = "\t")
  
  	#create txt file
	fileo<-file("iniParams.txt")
  
  	#write all information
	writeLines(c(l1, l2, l3, l4, l5, l6, 
		l7, l8, l9, l10, 0), fileo)
  
	close(fileo)
		
	#run ME2
    system("me2wopt.exe PMF_bs2",
    	show.output.on.console = F, wait = T)
	
}









###
# Function to apply PMF/ME2
# 
# pref is prefix for data
# ndays is number of days
# nfact is number of sources
# ncons is number of constituents (including PM25)
# nruns is number of PMF runs to complete
readME2 <- function(pref, ndays, nfact, ncons, nruns){
  
	#read in data
	all <- read.table(paste(pref, ".dat", sep= ""),
		skip = 2, sep = "", strip.white = T)
  
  	#check output has correct number of lines
  	neach <- ncons + ndays + 1
  	if((nruns * neach) != nrow(all)){
		print(neach)
  		print(nrow(all))

		browser()
	}
  
  
  
	#separate by PMF run, name columns as factors
	all <- data.frame(rep(seq(1 : nruns),
		each = neach), all)
	colnames(all) <- c("nrun", 
		paste("factor", seq(1 : nfact), sep = ""))
  
	#find best PMF run
	diag <- read.table(paste(pref, ".txt", sep = ""),
		skip = 2, sep = "", fill = T, 
		colClasses = c("character"))
	tasks <- which(diag[, 1] == "task") + 1
	wrun <- which.min(as.numeric(diag[tasks, 2]))
  
	#Select best PMF run
	all <- all[which(all[, 1] == wrun), ]
	#Get concentrations
	contr <- all[1 : ndays, -1]
  
	#Get profiles: remove spacer row
	prof <- all[(ndays + 1) : (neach - 1), -1]
	profo <- prof
  
  
	#remove PM25
	prof <- prof[-1, ]
	#get column sums
	matcs <- matrix(rep(colSums(prof),
		nrow(prof)), nrow = nrow(prof),
		byrow = T)
	#make sum to 1
	prof <- prof / matcs


	list(contr = contr, prof = prof, 
  		profo = profo)
}


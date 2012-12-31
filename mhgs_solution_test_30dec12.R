####
# function to test R gibbs
####

dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/data/"
dir1 <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/"


namef <- "woodiedus"
namef <- "woo"	

timedat <- function(namef) {
	sampl <- sample(seq(1, 300), 1)
	ns <- 2
	ps <- 20
	
	dat <- read.csv(file.path(dir, namef, 
		paste(namef,"_" ,
		sampl, ".csv", sep = "") ))
		
	mdls <- read.csv(file.path(dir, namef, 
		paste("mdls_", namef, "_", sampl,
		"_n", ns, "_p", ps, ".csv", sep = "")))
		
	test <- getimpdat(dat, mdls, N = 5, 
		burnin = 1)	
	}	

for ( i in 1 : 9) {
	namef <- names[i]
	
	x[i, 1] <- proc.time()
	timedat(namef)
	x[i, 2] <- proc.time()
}


Rprof()
test <- getimpdat(dat[1 : 300, ], mdls[1 : 300, ], N = 5, 
	burnin = 1)
summaryRprof()	

	
unlist?
t(yobs - mnobs) 	
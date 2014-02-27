####
# function to test R gibbs
####

dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/data/"
dir1 <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/"

library(ProjectTemplate)
setwd(dir1)
load.project()
run.project()

source("/Users/jennakrall/Dropbox/MDL_sourceapp/git-gibbs-truncnorm/mhgs_solution_8aug12.R")


namef <- "woodiedus"
namef <- "woo"	

timedat <- function(namef) {
	sampl <- sample(seq(1, 300), 1)
	ns <- 11
	ps <- 50
	
	dat <- read.csv(file.path(dir, namef, 
		paste(namef,"_" ,
		sampl, ".csv", sep = "") ))
		
	mdls <- read.csv(file.path(dir, namef, 
		paste("mdls_", namef, "_", sampl,
		"_n", ns, "_p", ps, ".csv", sep = "")))
		
	test <- getimpdat(dat, mdls, N = 10, 
		burnin = 3)	
	}	


namesSOURCE <- sapply(trues, function(x) 
	paste(substr(x, 1, 3), collapse = ""))


#### Time for 20 runs (*2500 for 50000)

x <- matrix(nrow = 9, ncol = 6)
for ( i in 1 : 9) {
	print(i)
	namef <- namesSOURCE[i]
	
	x[i, 1:3] <- proc.time()[1:3]
	timedat(namef)
	x[i, 4:6] <- proc.time()[1:3]
}

times <- x
save(times, file = file.path(dir,"timed_20_31dec12.RData"))

#hours to complete 1 job
hours <- (x[,4] - x[,1]) * 2500/(60*60)  #~ 110 hours per MCMC (50000 iter)
hours <- (x[,4] - x[,1]) * 500/(60*60)  #~ 20 hours per MCMC (10000 iter)
#years to complete all
hours * 300 * 12 / 24 / 365


Rprof()
test <- getimpdat(dat[1 : 300, ], mdls[1 : 300, ], N = 5, 
	burnin = 1)
summaryRprof()	

	
unlist?
t(yobs - mnobs) 	



#put stuff on enigma
sample(c(1: 100000), 1)
set.seed(39045)
x <- sample(c(1: 300), 6, replace = T)
# > x
# [1] 264  42   2 193 128   9
#use first three, woodiedust_n11_p50
#use second three, dusvehcoadie_n11_p80



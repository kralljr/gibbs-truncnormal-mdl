#####
# file to create BDL data sets
#
# 12/15/12
#####


#home
home.dir1 <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/"
home.dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/data/"
library(ProjectTemplate)

#load project
setwd(home.dir1)
load.project()
run.project()
source(file.path( home.dir1, "src/MDL_replace_5mar12.R"))
setwd(home.dir)


#set percs and nums                                                                                          
percs <-  c(.2, .5, .8, .9)
nums <- c(2, 11, 23)
pnames <- c(20, 50, 80, 90)

namessource <- sapply(trues, function(x) {
	paste(substr(x, 1, 3), collapse = "")
})               




mdl.out.createdat <- function(datall, perc) {
	
	ncols <- ncol(datall) * 2
	seqs <- seq(1, ncols, by = 2)
	
	out <- as.matrix(mdl.perc.below.out(datall, perc))

	out2 <- out[,seqs]

	list(out2,out[,-seqs]) 
}


#for each source combination
for ( i in 1 : 9) {
	#setwd
	setwd(file.path(home.dir, namessource[i]))
	
	#set up mdl matrix
	#n1 <- paste(namessource[i], "_", 1, ".csv", sep = "")
	#dat1 <- read.csv(n1)
	# mdlsout <- matrix(nrow = 300 * length(nums) * length(percs),
		# ncol = ncol(dat1))
		
	#for each simulation
	for ( l in 1 : 300 ) {
		
		#read data
		n1 <- paste(namessource[i], "_", l, ".csv", sep = "")
		dat <- read.csv(n1)
	
		#for each number of missing
		for ( k in 1 : length(nums)) {
			#for each percent
			for (j in 1 : length(percs)) {
				
				#read data
				
				
				
				#get correct percent	
				perces <- c( rep(percs[j], nums[k]),
						rep(0, (ncol(dat) - nums[k] )))
				#randomly sample which constitunts BDL to match nums
				persample <- sample(perces, replace = FALSE)
	
	
	
				#get MDLs and eliminate BDL
				temp <- mdl.out.createdat(dat, persample)
				#don't need, dat is unchanged
				#newdat <- temp[[1]]
				mdls <- temp[[2]]
				# mdlsout[m, ] <-  colMeans(mdls, na.rm = T)
				
				#save data
				fn <- paste( namessource[i], "_", l, 
					"_n", nums[k], "_p", pnames[j], ".csv", 
					sep = "")
				#dat has not changed!  use regular dat
				#write.csv(dat, file = fn , row.names = F)
				
				fn1 <- paste("mdls_", fn, sep = "")
				write.csv(mdls, file = fn1 , row.names = F)
				
				#increase row for mdls
				# m <- m + 1
				}
			}
		}
	# fn <- paste("MDLS_", namessource[i], ".csv", sep = "")
	# write.csv(mdlsout, file = fn, row.names = F)
	
}
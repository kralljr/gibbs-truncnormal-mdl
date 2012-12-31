#####
# File to create datasets 
####
# 12/14/12
####

library(ProjectTemplate)
#home.dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012"
home.dir <- "C:/Users/jkrall1/Dropbox/MDL_sourceapp/MDL_project_oct2012"
#home.dir <- "C:/Users/Jenna Krall/Dropbox/MDL_sourceapp/MDL_project_oct2012/"
setwd(home.dir)
load.project()
run.project()

source(file.path(home.dir, "src/PCA_out_pmf_9oct12.R"))
source(file.path(home.dir, "src/PCA_out_pmf_inner_10jul12.R"))
source(file.path(home.dir, "src/absolutepca_5mar12.R"))
source(file.path(home.dir, "src/useME2_11jun12.R"))


#copied from avgdat
createdataset <- function(sourcedat, nsource, trues, 
	k, sds = 0.01,  gwd = gwd, unc = F, name = F) {
	
	#find data, and see if classifies as true without MDL
	td1 <- datafun(sourcedat, nsource = nsource, sds = sds)
	td <- td1[[1]]
	totals <- rowSums(td)
	outT <- SAtrue(td, trues = trues, unc, nsource, 
		sameth = "APCA",
		k2 = k, totals = totals,cat=TRUE,rule=FALSE,gwd=gwd)
	outT2 <- SAtrue(td,trues=trues,unc,nsource,sameth = "PMF",
		k2=k,totals=totals,cat=TRUE,rule=FALSE,gwd=gwd)
 
        #length of correct matches to truth
	 lenw <- length(which(unique(outT[[2]])%in% trues))
     lenw2 <- length(which(unique(outT2[[2]])%in% trues))
     lenw <- min(lenw, lenw2)
                        
     niter <- 1
     
     #repeat until data classifies as truth         
     while (((lenw < nsource) & (niter < 300))) {
		td1 <- datafun(sourcedat,nsource=nsource,sds=sds)
        td <- td1[[1]]
        totals <- rowSums(td)
        outT <- SAtrue(td,trues=trues,unc,nsource,sameth = "APCA",k2=k,
              	totals=totals,cat=TRUE,rule=FALSE,gwd=gwd)
		outT2 <- SAtrue(td,trues=trues,unc,nsource,sameth = "PMF",
			k2=k,totals=totals,cat=TRUE,rule=FALSE,gwd=gwd)
		
		
        #length of correct matches to truth
        lenw1 <- length(which(unique(outT[[2]])%in% trues))
        lenw2 <- length(which(unique(outT2[[2]])%in% trues))
        lenw <- round(min(lenw1, lenw2))
		#print(lenw < nsource)

        niter <- niter+1
        if ( niter == 300 ) {
        	 print("error")
        	 print(c("get rid of ", name))
        	 }    
       }
                 
               
    td           
}               
               
               
namessource <- sapply(trues, function(x) {
	paste(substr(x, 1, 3), collapse = "")
})               



#create 300 complete data files for each source combination
for ( i in 1 : 9) {
	print(i)
	sourcedat <- getsourcedat( trues[[i]], profs, mnvar)
	
	setwd(file.path(home.dir, "data", namessource[i]))
	
	ks <- ksvec[i]
	set.seed(10)
	for ( j in 1 : 300) {
		# print(j)
		dat <- createdataset(sourcedat, nsources[i], trues[[i]], 
			 k = ks, name = j,
      	#gwd = "C:/Users/Jenna Krall/Dropbox/MDL_sourceapp/PMF/tryme/")
			gwd = "C:/Users/jkrall1/Documents/me2" )
			 
		filename <- paste(namessource[i],"_", j, ".csv", sep = "")	 
		write.csv(dat, file = filename, row.names = F) 
	}
	
}	




#####
# TEST
i <- 2
j <- 2
setwd(file.path(home.dir, "data", namessource[i]))
filename <- paste(namessource[i],"_", j, ".csv", sep = "")	
test <- read.csv(filename)
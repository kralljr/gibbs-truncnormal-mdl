# For each simulated dataset, each censoring adjustment method, run APCA

#load library
library(handles)
library(ncdf4)


#print and save arguments passed from BASH/python?
print(commandArgs(TRUE))
args <- commandArgs(TRUE)

# THIS IS WHAT I NEED TO PASS (with file paths) TO THIS FILE
# modelfile <- "models_woodiedus_1_n2_p10.nc"
# datafile <- "woodiedus_1.csv"
# savefile <- "APCA_woodiedus_model_1_n2_p10.RData"




filen <- strsplit(modelfile, "_")[[1]]

#simulation iteration 1-300
iter <- as.numeric(filen[[3]])
#source combo: woodiedus, dusvehdie, woodiedusvehcoa
source <- filen[[2]]
#number of censored constituents
ncons <- filen[[4]]
ncons <- as.numeric(substr(filen[[4]], 2, nchar(filen[[4]])))
#percent censored
censor <- strsplit(filen[[5]], "\\.")[[1]][1]
censor <- as.numeric(substr(censor, 2, 3))




#assign source details
if(source == "woodiedus" | source == "dusvehdie") {
	nsources <- 3
	k <- 1
	type1 <- ifelse(source == "woodiedus", 1, 2)
	if(source == "woodiedus") {
		sourceall <- c("wood", "diesel", "dust")
	}else{
		sourceall <- c("dust", "vehicle", "diesel")
		}
}else{
	nsources <- 5	
	k <- 20
	type1 <- 3
	sourceall <- c("wood", "diesel", "dust", "vehicle", "coal")
}






#read in data
# NEEDS TO BE MODIFIED FOR MODEL DATA
ncdat <- nc_open(modelfile)
data <- ncvar_get(ncdat, "data")
data <- aperm(data, c(2, 1, 3))
nc_close(ncdat)


#set seed
set.seed(15)
seeds <- sample(seq(1, 1000000), 300 * 3 * 5)
seedmat <- array(seeds, dim = c(300, 3, 5))
ncons1 <- c(2, 11)
censor1 <- c(20, 50, 80)
numper <- paste0(rep(ncons1, each = 3), rep(censor1, 2))[-6]
whcen <- which(numper == paste0(ncons, censor))
seed1 <- seedmat[as.numeric(iter), type1, whcen]



#get totals
consdat <- read.csv(datafile) 
tots <- rowSums(consdat)


#run APCA for each of 100 files
set.seed(seed1)

l1 <- vector()
class <- matrix(nrow = 100, ncol = nsources)
means <- class
sds <- class
for(i in 1 : 100) {
	
	#run APCA
	apca1 <- SIMapca(data = data, tots = tots, 
		nsources = nsources, 
		adjust = NULL, k = k)
		
	#get classification	
	class1 <- apca1$class
	l1[i] <- length(which(class1 %in% sourceall))
	#match classification to truth
	match1 <- match(sourceall, class1)
	
	#reorder and save results 
	means[i, ] <-apca1$means[match1]
	sds[i, ] <- apca1$sd[match1]
	class[i, ] <- apca1$class	
		
}


#combine output
l1 <- mean(l1, na.rm = T, trim = 0.2)	
	

apca.res <- list(means = means, sd = sd, 
	class = class, l1 = l1)



#get output format
save(apca.res, file = savefile)

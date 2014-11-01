# For each simulated dataset, each censoring adjustment method, run APCA



#print and save arguments passed from BASH
# list of length 6
# args[[1]] is source abbreviation
# args[[2]] is simulation iteration 1-300
# args[[3]] is number of constituents censored
# args[[4]] is percent censored
# args[[5]] is filepath to data
# args[[6]] is model iteration 1-100
print(commandArgs(TRUE))
args <- commandArgs(TRUE)

#args should include file name of data

#load library
library(handles)


#source combo e.g. woodiedus
sabb <- args[[1]]
#number of sources
if(sabb == "wdd" | sabb == "dvd") {
	source <- ifelse(sabb == "wdd", 
		"woodiedus", "dusvehdie")
	nsources <- 3
	k <- 1
	type1 <- ifelse(sabb == "wdd", 1, 2)
}else{
	nsources <- 5	
	k <- 20
	type1 <- 3
	source <- "woodiedusvehcoa"
}
#simulation number e.g. 124
iter <- args[[2]]
#how many constituents censored e.g. 2
ncons <- args[[3]]
#what percent censored e.g. 50
censor <- args[[4]]
#adjust: NULL 
fp <- args[[5]]
modelsim <- as.numeric(args[[6]])






#read in data
# NEEDS TO BE MODIFIED FOR MODEL DATA
filename <- paste0(source, "_", iter, ".csv")
fp <- file.path(fp, source)
data <- read.csv(file.path(fp, filename))




#set seed
set.seed(15)
seeds <- sample(seq(1, 1000000), 300 * 3 * 100 * 5)
seedmat <- array(seeds, dim = c(300, 3, 100, 5))
ncons1 <- c(2, 11)
censor1 <- c(20, 50, 80)
numper <- paste0(rep(ncons1, each = 3), rep(censor1, 2))[-6]
whcen <- which(numper == paste0(ncons, censor))
seed1 <- seedmat[as.numeric(iter), type1, modelsim, whcen]



#run APCA
set.seed(seed1)
tots <- rowSums(data)
apca.res <- SIMapca(data = data, tots = tots, 
	nsources = nsources, 
	adjust = NULL, k = k)


#save output
name1 <- paste0("APCA_", source, "_model_", iter,  
	"_n", ncons, "_p", censor)
name1 <- paste0(name1, ".RData")



#get output format
fp1 <- file.path(getwd(), "rdata")
fp1 <- file.path(fp1, name1)
print(fp1)

save(apca.res, file = fp1)

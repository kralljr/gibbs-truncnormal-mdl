# For each simulated dataset, each censoring adjustment method, run APCA

#load library
library(handles)
library(ncdf)


#print and save arguments passed from BASH/python?
print(commandArgs(TRUE))
args <- commandArgs(TRUE)

# THIS IS WHAT I NEED TO PASS (with file paths) TO THIS FILE
# datafile <- "woodiedus_1.csv"
# modelfile <- "models_woodiedus_1_n2_p10.nc"
# savefile <- "APCA_woodiedus_model_1_n2_p10.RData"
datafile <- args[[1]]
modelfile <- args[[2]]
savefile <- args[[3]]

filen <- strsplit(modelfile, "_")[[1]]

#simulation iteration 1-300
iter <- as.numeric(filen[[3]])
#source combo: woodiedus, dusvehdie, woodiedusvehcoa
source_combination <- filen[[2]]
#number of censored constituents
ncons <- filen[[4]]
ncons <- as.numeric(substr(filen[[4]], 2, nchar(filen[[4]])))
#percent censored
censor <- strsplit(filen[[5]], "\\.")[[1]][1]
censor <- as.numeric(substr(censor, 2, 3))




#assign source details
if (source_combination == "woodiedus" | source_combination == "dusvehdie") {
    k <- 1
    type1 <- ifelse(source_combination == "woodiedus", 1, 2)
    if (source_combination == "woodiedus") {
        sourceall <- c("wood", "diesel", "dust")
    } else {
        sourceall <- c("dust", "vehicle", "diesel")
    }
} else {
    k <- 20
    type1 <- 3
    sourceall <- c("wood", "diesel", "dust", "vehicle", "coal")
}

nsources <- length(sourceall)






#read in data
# NEEDS TO BE MODIFIED FOR MODEL DATA
ncdat <- open.ncdf(modelfile)
data <- aperm(get.var.ncdf(ncdat, "data"), c(2, 1, 3))
close.ncdf(ncdat)


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


#run APCA for each of draw
set.seed(seed1)

ndraws <- dim(data)[3]
l1 <- vector(length=ndraws)
class <- matrix(nrow = ndraws, ncol = nsources)
means <- class
sds <- class
for(i in 1 : ndraws) {
    #run APCA
    apca1 <- SIMapca(data = data[,, i], tots = tots, nsources = nsources,
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
apca.res <- list(means = means, sds = sds, class = class, l1 = l1)



#get output format
save(apca.res, file = savefile)

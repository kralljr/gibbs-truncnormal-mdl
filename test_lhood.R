#####################################
# Test file for likelihood method 
######################################
	
# source likelihood code
source("lhood.R")

# load simulated test data
# dat is a matrix of 1000 days by 23 chemical constituents of pm2.5
dat <- read.csv("test_data.csv")
# mdls is a matrix of 1000 X 23, which gives 
# the MDL on that day for that constituent
mdls <- read.csv("test_mdls.csv")

#set up initial guess
guess <- list()
guess[[1]] <- dat
guess[[2]] <- colMeans(dat)
guess[[3]] <- cov(dat)


#run likelihood based method for multiple imputation
test <- lhood(dat, mdls, guess, burnin = 1, N = 10)


#####
# or run in C (creates .nc file with output)
runGIBBSc(dat = "test_data.csv", mdls = "test_mdls.csv", 
	outfile = "testdat.nc", seed = 10)



#analyze output
#find posterior means of theta
apply(test[[2]], 1, mean)
#find posterior means of covariance sigma
apply(test[[3]], c(1, 2), mean)
#find posterior mean of data
datmean <- apply(test[[1]], c(1, 2), mean)
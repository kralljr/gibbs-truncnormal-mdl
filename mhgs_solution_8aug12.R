#####################################
# MCMCfunctions for imputing BDL values
#
# Edited 12/12/12 for speed 
#
#
######################################
	
	
	
###############################
# MCMC setup
###############################

# library(truncnorm)
library(msm)
library(tmvtnorm)
library(MCMCpack)





	
##########################################
#get new guess for missing data: trunc normal
#arguments 
	#datk: logged data vector 1 X N_cons
	#gthet: guess vector for mean (N_cons)
	#gsig: guess matrix for covar (N_cons X N_cons)
	#wh: which column is missing (1:N_cons)
##########################################
impmisssingle <- function(datk, gthet, gsig, wh) {
	#for each day

	cov11 <- gsig[-wh, -wh]
	cov00 <- gsig[wh, wh]
	cov01 <- gsig[wh, -wh]
	
	mnobs <- gthet[-wh]
	mnmiss <- gthet[wh]
	yobs <- datk[-wh]
	
	mn <- mnmiss + cov01 %*% solve(cov11, yobs - mnobs)
	ss <- backsolve(chol(cov11), cov01, transpose = T)
	ss <- crossprod(ss)
	var <- cov00 - ss

	list(mn,var)
}
	
	
	
	
		
	
	
	
	
#########	#########		
#########	
#####
#proposal functions
#####
#########	
#########	#########	
	

##########################################
#get new guess for missing data: trunc normal
#arguments 
	#dat: logged data matrix N_days X N_cons
	#mdls: logged data matrix N_days X N_cons of mdls
	#nbdlmat: binary matrix N_days X N_cons, 
		#1 = below DL, 0 above DL
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
##########################################
ymissfun <- function(dat, mdls, nbdls, guessvec,
	wh_miss) {
	
	#set old estimates
	gsig <- guessvec[[3]]		
	gthet <- guessvec[[2]]
	gdat <- guessvec[[1]]
	
	
	#set up for univariate
	m <- 1
	# for each day with BDL
	for (k in 1 : nrow(dat)) {

		if(length(wh_miss[[k]]) > 0) {
		
			#for each missing value on day i
			for (j in 1 : length(wh_miss[[k]])) {
				
				#find conditional mean/var
				mnv <- impmisssingle(gdat[k, ], gthet, 
					gsig, wh_miss[[k]][j])
	
				#propose truncated normal (0 to mdl)	 
					#log scale with cond mean and var		
				newymiss1 <- rtnorm(1, 
					upper = mdls[k, wh_miss[[k]][j]],
					mean = mnv[[1]], sd = sqrt(mnv[[2]]))
				if( is.na(newymiss1)) {browser()}	
		
				#update guess
				guessvec[[1]][k, wh_miss[[k]][j]] <- newymiss1
		
				
			
			}#end loop over col (j)
		
		} # ELSE DO NOTHING 
		
	}#end loop over row (k)
	
	guessvec
}











##########################################
#get new guess for mean: normal
#arguments 
	#dat: logged data matrix N_days X N_cons
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#prv: optional prior for degrees of freedom
	#prS: optional prior for scale matrix
		# (N_cons X N_cons, pd matrix)
##########################################
thetfun<-function(dat, guessvec,
	prmean = rep(0, ncol(dat)), 
	prsig = rep(10^5,ncol(dat))){
		
	#set old estimates
	gsig <- guessvec[[3]]
	gthet <- guessvec[[2]]
	gdat <- guessvec[[1]]

	#Prior: assume independence
	sprsigmat <- diag(1/prsig)

	#Invert current guess
	ssig <- chol(gsig)
	ssig <- chol2inv(ssig)
	
	#Get variance: mu | X, sigma ~ MVN(mns, vars)
	vars <- chol(sprsigmat + nrow(gdat) * ssig)
	vars <- chol2inv(vars)
	
	#Get mean
	mns <- vars %*% (sprsigmat %*% prmean + 
		nrow(gdat)*ssig %*% colMeans(gdat, na.rm = T))
	
	#draw new guess	
	newthet <- rmvnorm(1, mean = mns, sigma = vars)

	#update guess vector
	guessvec[[2]] <- newthet

	
	guessvec
}










	
##########################################
#get new guess for covariance:inverse wishart
#arguments 
	#dat: logged data matrix N_days X N_cons
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#prv: optional prior for degrees of freedom
	#prS: optional prior for scale matrix
		# (N_cons X N_cons, pd matrix)
##########################################
sigfun <- function(dat, guessvec, prv="no prior", prS="no prior") {
	#print("in sigfun")
	
	#set old estimates
	gphi <- guessvec[[3]]
	gthet <- guessvec[[2]]
	gdat <- guessvec[[1]]
	

	prv <- ncol(dat) + 1

	#if gdat is MVN, covar ~ inv wish (S2 = A + Psi , v2 = v + N_days)
	v2 <- nrow(gdat) + prv
	
	#S2 calculation
	swp <- sweep(gdat, 2, gthet)
	arswp <- crossprod(swp)
	Psi <- diag(ncol(gdat))
	S2 <- Psi + arswp
	
	
	newsig <- riwish(v = v2, S = S2)

	#update guess for covariance
	guessvec[[3]] <- newsig

	guessvec
}











	
	#########	#########	
#########		
#######
# outer functions
######
#########	
#########	#########	


##########################################
#update to new guess
	#arguments dat: logged data matrix N_days X N_cons
	#nbdlmat: binary matrix N_days X N_cons, 
			#1 = below DL, 0 above DL
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#wh which fill in 
		#1: missing data
		#2: mean
		#3: covariance
	#wh_miss list of length N_days with 
		#column numbers of missing
##########################################

gibbsfun <- function(dat, nbdlmat, guessvec,
	wh, mdls, wh_miss) {
	#we have initial guesses 
	
	if (wh == 1) {
		#update missing
		guessvec <- ymissfun(dat, mdls, nbdlmat, guessvec,
			wh_miss)
	} else if (wh == 2) {
		#update mean
		guessvec <- thetfun(dat, guessvec)
		
	} else if (wh == 3) {
		#update covariance
		guessvec <- sigfun(dat, guessvec)
		
	}
	
	guessvec
	}
	
	
	
	
	
	
	
	
	
##########################################
#run whole MCMC
#arguments 
	#dat: logged data matrix N_days X N_cons
	#mdls: logged data matrix N_days X N_cons of mdls
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#N is number of iterations
##########################################
mhwithings <- function(dat, mdls, nbdlsmat, 
	wh_miss, guessvec, N){
	

	#create arrays for output
	gymiss <- array(dim = c(nrow(dat), ncol(dat), N + 1))
	gthet <-  array(dim = c(ncol(dat), N + 1))
	gsig <- array(dim = c(ncol(dat), ncol(dat), N + 1))
	
	
	#set first guesses
	gymiss[, , 1] <- guessvec[[1]]
	gthet[, 1] <- guessvec[[2]]
	gsig[, , 1] <- guessvec[[3]]
	
	
	drops <- 0
	#for each iteration (N large)
	for(i in 1:N){
		# print(i)
		j <- 1
		while(j<4){
			#update f/lam

			gibbs<-gibbsfun(dat, nbdlmat, guessvec, j, 
				mdls, wh_miss)
		
			# if(j==3){
				
				# if(gibbs[[2]]==FALSE){

					# drops <- c(drops, i)
				# }else{
					# guessvec <- gibbs[[1]]
				# }
			# }else{
				guessvec <- gibbs
				# }
	
	
			#update guesses, new last
			gymiss[, , i + 1] <- guessvec[[1]]
			gthet[, i + 1] <- guessvec[[2]]
			gsig[, , i + 1] <- guessvec[[3]]
			
			j <- j + 1

		}
	}
	drops <- drops[ -1 ]
	
	#create output
	listout <- list(gymiss, gthet, gsig, drops)
	names(listout) <- c("gymiss", "gthet", "gsig", "drops")
	
	
	listout
}

















#######################################
#######################################
#######################################
#######################################
#######################################

###NOTE: DO NOT NEED THIS TO BE FASTER, THIS
###IS WRAPPER FOR MCMC
##########################################
#######################################
#########################################
##########################################
##################################
#file to wrap MCMC and get imputed
#arguments 
	#dats: unlogged data matrix N_days X N_cons
	#mdls: unlogged data matrix N_days X N_cons of mdls
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#N is number of iterations
	#skips is the number of skips
	#burnin is the length of the burn in period
##########################################
getimpdat <- function(datscomplete, mdls, N = 3000, 
	skips = 5, burnin = 1000) {
		
	#get logged data	
	logdat <- log(datscomplete)
	logmdl <- log(mdls)
	
	
	#create binary matrix of which BDL
	nbdlmat <- 1 * (logdat >= logmdl)
	for ( i in 1 : ncol(logdat)) {
		#find which BDL
		whBDL <- which(nbdlmat[, i] == 0)
		
		#replace BDL with 1/2MDL
		if (length(whBDL) > 0) {
			logdat[whBDL, i] <- 1/2 * logmdl[whBDL, i]
		}
	}
	
	
	#set initial guesses
	guessvec <- vector(mode = "list" , length = 3)
	guessvec[[1]] <- logdat
	guessvec[[2]] <- colMeans(logdat)
	guessvec[[3]] <- cov(logdat)
	
	#get indices of which BDL and which above DL for each day
	wh_miss <- vector(mode = "list", length = nrow(logdat))
	for (i in 1 : nrow(logdat)) {
		wh_miss[[i]] <- which(nbdlmat[i, ] == 0)
	}

	#do MCMC
	mhtest <- mhwithings(logdat, logmdl, nbdlmat, wh_miss, 
		guessvec, N)
	
	#get skips
	# seqs <- seq(1, N-burnin, by = skips)
	
	datMed <- logdat
	datMean <- logdat
	datDraw <- array( dim = c(nrow(logdat), ncol(logdat), 10))
	#for each day
	for (k in 1 : nrow(logdat)) {
		
		#if at least one missing
		if( length(wh_miss[[k]]) > 0 ) {
	
			#if more than 1 missing
			if ( length(wh_miss[[k]]) > 1 ) {
				
				
				#get all data for that day, missing cons
				dat <- mhtest[[1]][k, wh_miss[[k]], ]
				#eliminate burnin
				dat <- dat[, burnin : (N + 1)]
				#undo correlated samples
				#dat <- dat[, seqs]
				

				datMed[k, wh_miss[[k]]] <- exp(apply(dat, 1, median, na.rm = T))
				datMean[k, wh_miss[[k]]] <- exp(apply(dat, 1, mean, na.rm = T))
				
				#get random draws for multiple imputation
				ls <- sample(seq(1, length(dat)), 10, replace = T)
				for( l in 1 : 10) {
					datDraw[k, wh_miss[[k]], l] <- exp(dat[, ls[l]])
					}
					
			#if only 1 missing that day		
			} else {
				dat <- mhtest[[1]][k, wh_miss[[k]], ]
				dat <- dat[ burnin : (N + 1)]
				
				
				#get median and mean of chain
				datMed[k, wh_miss[[k]]] <- exp(median(dat, na.rm = T))
				datMean[k, wh_miss[[k]]] <- exp(mean(dat, na.rm = T))
				
				#get random draws for multiple imputation
				ls <- sample(seq(1, length(dat)), 10, replace = T)
				for( l in 1 : 10) {
					datDraw[k, wh_miss[[k]], l] <- exp(dat[ls[l]])
					}
			}
			

		
		
		} #ELSE DO NOTHING
	}
	
	#exponeniate data for output
	out <- list(datMed, datMean, datDraw)
	names(out) <- c("median", "mean", "draws")
	
	out
}


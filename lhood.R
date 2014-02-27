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
impmisssingle <- function(datk, gthet, gsig, gsig.misscol, wh,
	mnzero) {
	#for each day

#	cov11 <- gsig[-wh, -wh]
	cov00 <- gsig[wh, wh]
	cov01 <- gsig[wh, -wh]
	
	mnobs <- gthet[-wh]
	mnmiss <- gthet[wh]
	yobs <- datk[-wh]

	
	#mn <- mnmiss + cov01 %*% solve(cov11, t(yobs - mnobs))
	mn <- mnmiss + cov01 %*% gsig.misscol[, , wh] %*% 
		mnzero[-wh]
	#ss <- backsolve(chol(cov11), cov01, transpose = T)
	#ss <- crossprod(ss)
	#var <- cov00 - ss
	var <- cov00 - cov01 %*% gsig.misscol[, , wh] %*% cov01

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
	wh_miss, minmdls) {
	
	#set old estimates
	gsig <- guessvec[[3]]		
	gthet <- guessvec[[2]]
	gdat <- guessvec[[1]]
	
	
	#set up for univariate
	m <- 1
	# for each day with BDL
	
	#get inverse of gsig
	gsiginv <- chol2inv(chol(gsig))
	
	#which columns have missing
	wcolmiss <- which(colSums(nbdls) < nrow(nbdls))
	gsig.misscol <- array(dim = c(nrow(gsig) - 1, 
		ncol(gsig) - 1, ncol(gsig)))
		
	#for each column with missingness	
	for( i in 1 : length(wcolmiss)) {
		num <- wcolmiss[i]
		sigcut <- gsiginv[-num, -num] - gsiginv[num, -num] %*% 
			t(gsiginv[-num, num]) / gsiginv[num, num]
		#update array	
		gsig.misscol[, , num] <- sigcut
	}
	
	
	
	for (k in 1 : nrow(dat)) {
		
		mnzero <- t(gdat[k, ] - gthet)

		if(length(wh_miss[[k]]) > 0) {
		
			#for each missing value on day i
			for (j in 1 : length(wh_miss[[k]])) {
		
				#find conditional mean/var
				mnv <- impmisssingle(datk = gdat[k, ], 
					gthet = gthet, 
					gsig = gsig, gsig.misscol = gsig.misscol, 
					wh = wh_miss[[k]][j], mnzero = mnzero)
	
				#propose truncated normal (0 to mdl)	 
					#log scale with cond mean and var		
				newymiss1 <- rtnorm(1, lower = minmdls - 10, 
					upper = mdls[k, wh_miss[[k]][j]],
					mean = mnv[[1]], sd = sqrt(mnv[[2]]))
				if( is.na(newymiss1)) {browser()}	
		
				#update guess
				guessvec[[1]][k, wh_miss[[k]][j]] <- newymiss1
		
				
			
			}#end loop over col (j)
		
		} # ELSE DO NOTHING 
		
	}#end loop over row (k)
	
	list(guessvec, gsiginv)
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
		
	#get inverse of sigma
	ssig <- guessvec[[2]]
	
	#set old estimates
	guessvec <- guessvec[[1]]
	gsig <- guessvec[[3]]
	gthet <- guessvec[[2]]
	gdat <- guessvec[[1]]


	#Prior: assume independence
	sprsigmat <- diag(1/prsig)

	#Invert current guess
	# ssig <- chol(gsig)
	# ssig <- chol2inv(ssig)
	
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
	arswp <- crossprod(as.matrix(swp))
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
	wh, mdls, wh_miss, minmdls) {
	#we have initial guesses 
	if (wh == 1) {
		#update missing
		guessvec <- ymissfun(dat = dat, mdls = mdls, 
			nbdls = nbdlmat, guessvec = guessvec,
			wh_miss = wh_miss, minmdls = minmdls)
		# guessvec <- guessvec1[[1]]
		# gsiginv <- guessvec1[[2]]
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
	wh_miss, guessvec, burnin, N){


	#create arrays for output
	gymiss <- array(dim = c(nrow(dat), ncol(dat), N - burnin ))
	gthet <-  array(dim = c(ncol(dat), N - burnin))
	gsig <- array(dim = c(ncol(dat), ncol(dat), N - burnin))
	
	
	#set first guesses
	gymiss1 <- as.matrix(guessvec[[1]])
	gthet1 <- as.vector(guessvec[[2]])
	gsig1 <- as.matrix(guessvec[[3]])
	
	minmdls <- min(mdls)
	
	#for each iteration (N large)
	for(i in 1:N){
		print(i)
		j <- 1
		while(j<4){
			#update f/lam
			gibbs<-gibbsfun(dat = dat, nbdlmat = nbdlsmat, 
				guessvec = guessvec, wh = j, 
				mdls = mdls, wh_miss = wh_miss, minmdls = minmdls)
		

				if(j == 1) {
					guessvec <- gibbs[[1]]
				}else {
					
					guessvec <- gibbs
					}
	
	
			#update guesses, new last
			#if we are past the burnin period, save
			if (i > burnin) {
				gymiss[, , i - burnin] <- as.matrix(guessvec[[1]])
				gthet[, i - burnin] <- as.vector(guessvec[[2]])
				gsig[, , i - burnin] <- as.matrix(guessvec[[3]])
			}
			
			
			if(j == 1) {
				guessvec <- gibbs
				}
			j <- j + 1

		}
	}
	
	#create output
	listout <- list(gymiss, gthet, gsig)
	names(listout) <- c("gymiss", "gthet", "gsig")
	
	
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
	mhtest <- mhwithings(dat = logdat, mdls = logmdl, 
		nbdlsmat = nbdlmat, wh_miss = wh_miss, 
		guessvec = guessvec, burnin = burnin, N = N)
	
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
				
				#Do THIS INSIDE MHGS
				#eliminate burnin
				#dat <- dat[, burnin : (N + 1)]
				#undo correlated samples
				#dat <- dat[, seqs]
				

				datMed[k, wh_miss[[k]]] <- exp(apply(dat, 1, median, na.rm = T))
				datMean[k, wh_miss[[k]]] <- exp(apply(dat, 1, mean, na.rm = T))
				
				#get random draws for multiple imputation
				ls <- sample(seq(1, ncol(dat)), 10, replace = F)
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
				ls <- sample(seq(1, length(dat)), 10, replace = F)
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




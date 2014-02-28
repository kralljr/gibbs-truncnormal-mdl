#####################################
# Likelihood method for multiply imputing 
# censored observations below the minimum 
# detection limit (MDL)
######################################
	





	
##########################################
#run likelihood method MCMC
#arguments 
	#dat: logged data matrix N_days by N_cons
	#mdls: logged data matrix N_days by N_cons of mdls
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#N is number of iterations
##########################################
lhood <- function(dat, mdls,  
	guessvec, burnin = 10000, N = 50000){

	libmsm <- require(msm)
	libnorm <- require(tmvtnorm)
	libmcmc <- require(MCMCpack)
	if(!(libmsm & libnorm & libmcmc)) {
		stop("Need to install libraries 'msm' 'tmvtnorm' 'MCMCpack'")
	}
	
	
	#create binary matrix of which BDL
	nbdlsmat <- 1 * (dat >= mdls)
	
	#get indices of which BDL and which above DL for each day
	wh_miss <- vector(mode = "list", length = nrow(dat))
	for (i in 1 : nrow(dat)) {
		wh_miss[[i]] <- which(nbdlsmat[i, ] == 0)
	}

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











##########################################
#Function to run in C (C code must be compiled)
##########################################






	
##########################################
#run gibbs sampler in C
#arguments 
	#dat: logged data matrix N_days X N_cons
	#mdls: logged data matrix N_days by N_cons of mdls
	#outfilename: name of output file
	#niter: number of iterations for MCMC
	#burnin: number of burnin samples
	#ndraws: number of draws from the censored data to output
	#seed: random seed
##########################################

runGIBBSc <- function(dat, mdls, outfilename, niter = 100, 
	burnin = 50, ndraws = 1, seed) {
		
	comm <- paste("./gibbs -n", niter, "-b", burnin,
		"-d", ndraws, "-r", seed, dat, mdls, outfilename)
		
	system(comm)	 
	
}



















##########################################
#Internal functions
##########################################








##########################################
#update to new guess
	#arguments 
	#dat: logged data matrix N_days by N_cons
	#nbdlmat: binary matrix N_days by N_cons, 
			#1 = below DL, 0 above DL
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#wh which fill in 
		#1: missing data
		#2: mean
		#3: covariance
	#mdls: logged data matrix N_days by N_cons of mdls 
	#wh_miss list of length N_days with 
		#column numbers of missing
	#minmdls smallest mdl value	
##########################################
gibbsfun <- function(dat, nbdlmat, guessvec,
	wh, mdls, wh_miss, minmdls) {
		
	#we have initial guesses 
	if (wh == 1) {
		#update missing
		guessvec <- ymissfun(dat = dat, mdls = mdls, 
			nbdls = nbdlmat, guessvec = guessvec,
			wh_miss = wh_miss, minmdls = minmdls)
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
#get new guess for missing data: trunc normal
#arguments 
	#dat: logged data matrix N_days X N_cons
	#mdls: logged data matrix N_days X N_cons of mdls
	#nbdls: binary matrix N_days X N_cons, 
		#1 = below DL, 0 above DL
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#wh_miss list of length N_days with 
		#column numbers of missing
	#minmdls smallest mdl value	
##########################################
ymissfun <- function(dat, mdls, nbdls, guessvec,
	wh_miss, minmdls) {
	
	#set old estimates
	gsig <- guessvec[[3]]		
	gthet <- guessvec[[2]]
	gdat <- guessvec[[1]]
	
	
	#set up for univariate
	m <- 1

	
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
	
	
	# for each day with an observation below the MDL
	for (k in 1 : nrow(dat)) {
		
		#set up mean
		mnzero <- t(gdat[k, ] - gthet)

		#if which missing greater than 0
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
#get new guess for missing data: trunc normal
#arguments 
	#datk: logged data vector 1 X N_cons
	#gthet: guess vector for mean (N_cons)
	#gsig: guess matrix for covar (N_cons X N_cons)
	#gsig.misscol:
	#wh: which column is missing (1:N_cons)
	#mnzero: specified by ymissfun
##########################################
impmisssingle <- function(datk, gthet, gsig, gsig.misscol, wh,
	mnzero) {

	#set up current guesses
	cov00 <- gsig[wh, wh]
	cov01 <- gsig[wh, -wh]
	
	mnobs <- gthet[-wh]
	mnmiss <- gthet[wh]
	yobs <- datk[-wh]

	#get conditional mean and variance
	mn <- mnmiss + cov01 %*% gsig.misscol[, , wh] %*% 
		mnzero[-wh]

	var <- cov00 - cov01 %*% gsig.misscol[, , wh] %*% cov01

	list(mn,var)
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

##########################################
#get new guess for mean: normal
#arguments 
	#dat: logged data matrix N_days X N_cons
	#guessvec list of initial guesses
		#guessvec[[1]] missing
		#guessvec[[2]] mean
		#guessvec[[3]] covariance
	#prmean: optional prior for mean
	#prsig: optional prior for variance
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




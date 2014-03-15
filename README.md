gibbs-truncnormal-mdl
=====================

Likelihood-based method for imputing values below the detection limit.  The method assumes the complete, unobserved data are multivariate normal and uses an MCMC approach to estimate the parameters.  We use conjugate priors for the mean and covariance of the multivariate normal distribution and directly sample from the full conditional distributions using a Gibbs sampler.  The full conditional for the censored data is truncated normal, where the data are truncated above by the MDL.  The full model is listed in the files truncated_likelihood_model.  



Running the code in R
---------------------

Necessary packages (with their dependencies) include:
* MCMCpack (for riwish function)
* mvtnorm (for rmvnorm function)
* msm (for rtnorm function)

To run the MCMC you will need:
* data matrix of pollutant concentrations with number of days (rows) and number of chemical constituents (columns)
* MDL matrix in same format as the data giving the MDL for each day and constituent
* list of starting values for the mean, covariance, and censored concentrations

R code files:
* lhood.R includes all functions to run our likelihood based method
* test_lhood.R gives test code for an example dataset



Running the code in C
---------------------

The C code runs the same model as the R code, but has much faster run times.  

Steps for C code
* Compile C code
* Run from R using the function lhoodC in R in lhood.R

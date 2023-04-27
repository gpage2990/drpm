#-----------------------------------------------------------------------------------------------------


# R-wrapper that fits the global model with a number of options.
# 1 - temporal dependence in the likelihood
# 2 - temporal dependence in the latent level
# 3 - temporal dependenc in the partition model
# 4 - Spatial information in the partition model.

drpm_fit <- function(y,s_coords=NULL,
                     M=1,
                     initial_partition=NULL,
					           starting_alpha=0.5,
					           unit_specific_alpha = FALSE, # FALSE implies one alpha for all units
					           time_specific_alpha = FALSE, # FALSE implies one alpha for all time
					           alpha_0=FALSE, # TRUE means alpha=starting_alpha not updated
					           eta1_0=FALSE, # TRUE means that eta1 = 0 not updated
					           phi1_0=FALSE, # TRUE means that phi1 = 0 not updated
					           modelPriors=c(0,100^2,1,1,1,1),
					           alphaPriors=rbind(c(1,1)),  # I need to add this to the DRPM AR1 SPPM function
					           simpleModel = 0, # if 1 then simple model in Sally's simulation.
					           theta_tau2 =c(0, 2), # this is only used if simpleModel = 1
					           SpatialCohesion=4,
					           cParms=c(0, 1, 2, 1),
					           mh=c(0.5, 1, 0.1, 0.1, 0.1),
					           verbose=FALSE,
					           draws=1100,burn=100,thin=1){

	nout <- (draws-burn)/thin

	nsubject = nrow(y)
	if(!is.null(initial_partition)){
	  y <- cbind(rep(0, nsubject), y)
	  if(length(initial_partition) != nsubject) stop("initial partition is not the correct length")
	}
	ntime = ncol(y)

	s1 <- s_coords[,1]; s2 <- s_coords[,2]

	Si <- llike <- fitted <- matrix(1, nrow=nout, ncol=ntime*nsubject)
	gamma <- matrix(0, nrow=nout, ncol=ntime*nsubject)
	mu <- sig2 <- matrix(1, nrow=nout, ncol=ntime*nsubject)
	theta <- tau2 <- alpha_out <- matrix(0.5, nrow=nout, ncol=ntime);
	eta1 <- matrix(0, nrow=nout, ncol=nsubject);
	phi0 <- phi1 <- lam2 <- rep(0, nout)
	lpml <- waic <- rep(0,1)

	sPPM <- ifelse(is.null(s1), FALSE, TRUE)
	if(unit_specific_alpha) alpha_out <- matrix(0.5, nrow=nout, ncol=ntime*nsubject)
  if(unit_specific_alpha & nrow(alphaPriors) != nsubject){
    stop("you have selected unit-specific-alpha option but did not supply correct alpha prior")
  }
	cat("sPPM = ", sPPM, "\n")

#if(!sPPM){
#	C.out <- .C("mcmc_drpm_ar1",
#	C.out <- .C("drpm_ar1_sppm",
#	C.out <- .C("drpm_ar1_sppm2",
#      as.integer(draws), as.integer(burn),as.integer(thin),
#      as.integer(nsubject),as.integer(ntime),
#      as.double(t(y)), as.double(s1), as.double(s2),
#      as.double(M),
#      as.double(modelPriors),as.integer(time_specific_alpha),
#      as.integer(alpha_0), as.integer(eta1_0), as.integer(phi1_0),
#      as.integer(sPPM), as.integer(SpatialCohesion), as.double(cParms),
#      as.double(mh), as.integer(verbose),
#      Si.draws = as.integer(Si),mu.draws = as.double(mu),
#      sig2.draws = as.double(sig2), eta1.draws = as.double(eta1),
#      theta.draws = as.double(theta), tau2.draws = as.double(tau2),
#      phi0.draws = as.double(phi0), phi1.draws = as.double(phi1),
#      lam2.draws = as.double(lam2), gamma.draws=as.integer(gamma),
#      alpha.draws = as.double(alpha_out), fitted.draws = as.double(fitted),
#      llike.draws=as.double(llike), lpml.out = as.double(lpml),
#      waic.out = as.double(waic))
#}


#if(sPPM){
	space_1 <- FALSE;
	alpha <- starting_alpha
	update_alpha <- ifelse(alpha_0==TRUE, 0, 1)
	update_eta1 <- ifelse(eta1_0==TRUE, 0, 1)
	update_phi1 <- ifelse(phi1_0==TRUE, 0, 1)

	ntime_out <- ntime
  cat("ntime = ", ntime, "\n")
  cat("ntime_out = ", ntime_out, "\n")
  if(!is.null(initial_partition)){
	  #cat(initial_partition, "\n")
	  C.out <- .C("informed_ar1_sppm",
	              as.integer(draws), as.integer(burn),as.integer(thin),
	              as.integer(nsubject),as.integer(ntime),
	              as.double(t(y)), as.double(s1), as.double(s2),
	              as.double(M), as.integer(initial_partition), as.double(alpha),
	              as.double(modelPriors),as.double(t(alphaPriors)),
	              as.integer(time_specific_alpha), as.integer(unit_specific_alpha),
	              as.integer(update_alpha), as.integer(update_eta1), as.integer(update_phi1),
	              as.integer(sPPM), as.integer(SpatialCohesion), as.double(cParms),
	              as.double(mh), as.integer(space_1),
	              as.integer(simpleModel), as.double(theta_tau2),
	              Si.draws = as.integer(Si),mu.draws = as.double(mu),
	              sig2.draws = as.double(sig2), eta1.draws = as.double(eta1),
	              theta.draws = as.double(theta), tau2.draws = as.double(tau2),
	              phi0.draws = as.double(phi0), phi1.draws = as.double(phi1),
	              lam2.draws = as.double(lam2), gamma.draws=as.integer(gamma),
	              alpha.draws = as.double(alpha_out),fitted.draws = as.double(fitted),
	              llike.draws=as.double(llike),lpml.out = as.double(lpml),
	              waic.out = as.double(waic))

	  ntime_out <- ntime - 1
	}

	if(is.null(initial_partition)){
	  initial_partition <- 0
	  C.out <- .C("drpm_ar1_sppm",
	              as.integer(draws), as.integer(burn),as.integer(thin),
	              as.integer(nsubject),as.integer(ntime),
	              as.double(t(y)), as.double(s1), as.double(s2),
	              as.double(M), as.double(alpha),
	              as.double(modelPriors),as.double(t(alphaPriors)),
	              as.integer(time_specific_alpha),
	              as.integer(update_alpha), as.integer(update_eta1), as.integer(update_phi1),
	              as.integer(sPPM), as.integer(SpatialCohesion), as.double(cParms),
	              as.double(mh), as.integer(space_1),
	              as.integer(simpleModel), as.double(theta_tau2),
	              Si.draws = as.integer(Si),mu.draws = as.double(mu),
	              sig2.draws = as.double(sig2), eta1.draws = as.double(eta1),
	              theta.draws = as.double(theta), tau2.draws = as.double(tau2),
	              phi0.draws = as.double(phi0), phi1.draws = as.double(phi1),
	              lam2.draws = as.double(lam2), gamma.draws=as.integer(gamma),
	              alpha.draws = as.double(alpha_out),fitted.draws = as.double(fitted),
	              llike.draws=as.double(llike),lpml.out = as.double(lpml),
	              waic.out = as.double(waic))
	}

#}

  cat("ntime = ", ntime, "\n")
  cat("ntime_out = ", ntime_out, "\n")

	out <- NULL
	out$Si <- array(C.out$Si.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  out$gamma <- array(C.out$gamma.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  out$mu <- array(C.out$mu.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  out$sig2 <- array(C.out$sig2.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  if(unit_specific_alpha) out$alpha <- array(C.out$alpha, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  if(!unit_specific_alpha) out$alpha <- matrix(C.out$alpha.draws,nrow=nout, byrow=TRUE)[,1:ntime_out,drop=FALSE]
  out$theta <- matrix(C.out$theta.draws, nrow=nout, byrow=TRUE)[,1:ntime_out,drop=FALSE]
  out$tau2 <- matrix(C.out$tau2.draws,nrow=nout, byrow=TRUE)[,1:ntime_out,drop=FALSE]

  out$eta1 <- matrix(C.out$eta1.draws,nrow=nout, byrow=TRUE)

	out$phi0 <- matrix(C.out$phi0.draws,nrow=nout, byrow=TRUE)
  out$phi1 <- matrix(C.out$phi1.draws,nrow=nout, byrow=TRUE)
  out$lam2 <- matrix(C.out$lam2.draws,nrow=nout, byrow=TRUE)

  out$llike <- array(C.out$llike.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  out$fitted <- array(C.out$fitted.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
  out$lpml <- C.out$lpml.out
  out$waic <- C.out$waic.out

  out$initial_partition = initial_partition

  if(simpleModel==1){
    out <- NULL
    out$Si <- array(C.out$Si.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
    out$gamma <- array(C.out$gamma.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
    if(unit_specific_alpha) out$alpha <- array(C.out$alpha, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
    if(!unit_specific_alpha) out$alpha <- matrix(C.out$alpha.draws,nrow=nout, byrow=TRUE)[,1:ntime_out,drop=FALSE]
    out$mu <- array(C.out$mu.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
    out$sig2 <- 1
    out$theta <- theta_tau2[1]
    out$tau2 <- theta_tau2[2]
    out$llike <- array(C.out$llike.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
    out$fitted <- array(C.out$fitted.draws, c(ntime,nsubject,nout))[1:ntime_out,,,drop=FALSE]
    out$lpml <- C.out$lpml.out
    out$waic <- C.out$waic.out
    out$initial_partition = initial_partition

  }

  out

}





#-----------------------------------------------------------------------------------------------------




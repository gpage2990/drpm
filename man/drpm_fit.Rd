\name{drpm_fit}
\alias{drpm_fit}
\title{Function that fits spatio-temporal hierarchical Gaussian model with a dependent random partition model}
\usage{
drpm_fit(y,
         s_coords=NULL,
         M=1,
         initial_partition = NULL,
         starting_alpha = 0.5,
         unit_specific_alpha = FALSE,
         time_specific_alpha = FALSE,
         alpha_0=FALSE,
         eta1_0=FALSE,
         phi1_0=FALSE,
         modelPriors=c(0,100^2,1,1,1,1),
         alphaPriors=rbind(c(1,1)),
         simpleModel = 0,
         theta_tau2 =c(0, 2),
         SpatialCohesion=4,
         cParms=c(0, 1, 2, 1),
         mh=c(0.5, 1, 0.1, 0.1, 0.1),
         verbose=FALSE,
         draws=1100,burn=100,thin=1)
}

\description{
\code{drpm_fit} is the main function used to fit model with Gaussian likelihood and spatio-temporal dependent random partition model.  An initial ``guess'' on the partition can be provided as well as prior information on the ``certainty'' associated with the initial guess.   In the model, temporal dependence could possibly appear in three places a) AR(1) process in data model, b) AR(1) process associated with mean model of atoms, c) temporal dependence in the partition.

}

\arguments{
\item{y}{An m x T numeric matrix containing the response values measured over time.}

\item{s_coords}{Two-column matrix containing spatial locations (i.e., longitude and lattitude).}

\item{M}{Parameter related to Dirichlet process scale or dispersion parameter.}

\item{initial_partition}{The initial partition that is illicited from an expert.}

\item{starting_alpha}{Initial value for alpha.  If alpha_0 = TRUE (see below), then alpha is fixed to this value.}

\item{unit_specific_alpha}{Logical that indicates if a unique alpha is to be included for each unit.}

\item{time_specific_alpha}{Logical that indicates if a unique alpha is to be included for each time point.}

\item{alpha_0}{Logical where TRUE indicates alpha_t = starting_alpha (i.e., alpha is not updated)}

\item{eta1_0}{Logical where TRUE indicates eta1 = 0 (i.e., observations are conditionally independent).}

\item{phi1_0}{Logical where TRUE indicates phi1 = 0 (i.e., iid model for atoms)}

\item{modelPriors}{Vector containing the following hierarchical model prior values in order
\itemize{
  \item {m0 - mean of phi0}
	\item {s20 - variance of phi0}
	\item {A - upper bound on sigma }
	\item {At - upper bound on tau }
	\item {Al - upper bound on lambda }
	\item {be - scale prior for xi}
  }
}
\item{alphaPriors}{Matrix containing beta prior values for the alpha parameters.  The columns correspond to the 1st and 2nd shape beta shape parameters.  The number of rows corresponds to the number of units.  If no initial partition is supplied then only the first row is used.
}
\item{simpleModel}{If 1 then the simple model in Sally's simulation is employed.}

\item{theta_tau2}{Only used if simpleModel=1.  This is a two-dimensional vector with the first entry corresponding to the value that is fixed to theta and the second to tau2}

\item{SpatialCohesion}{Scalar indicating what cohesion to use if spatial coordinates are provided
\itemize{
\item {3 - auxiliary similarity }
\item {4 - double dipper similarity }
  }
}

\item{cParms}{Vector containing spatial cohesion parameter values (see below for more details)
\itemize{
\item {mu0 - center of NNIG }
\item {k0 - scale parameter of gaussian part of NNIG }
\item {v0 - degrees of freedom for inverse-gamma part of NNIG}
\item {L0 - scale parameter for inverse-gamma part of NNIG}
  }
}
\item{mh}{Vector of gaussian standard deviations for metropolis updates of sigma2, tau, lambda, eta1, phi1.}

\item{draws}{Number of MCMC samples to collect}

\item{burn}{Number of the MCMC samples discarded in the burn-in phase of the sampler}

\item{thin}{The amount of thinning desired for the chain}

\item{verbose}{Logical indicating if MCMC progression should be print to screen}
}

\value{
This function returns, in a list, all MCMC interates for each model parameter.  In addition the LPML and WAIC model fit metric is provided.
}
\examples{


\dontrun{


   Still needs to be done.

  }


}

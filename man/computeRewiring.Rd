% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rmimp.R
\name{computeRewiring}
\alias{computeRewiring}
\title{Score wt and mt sequences for a pwm}
\usage{
computeRewiring(obj, mut_ps, prob.thresh = 0.5, log2.thresh = 1,
  include.cent = F, degenerate.pwms = F, .degenerate.groups = c("DE",
  "KR", "ILMV"))
}
\arguments{
\item{obj}{MIMP kinase object containing PWM, auc, GMM parameters, family name, etc.}

\item{mut_ps}{psnvs data frame containing wt and mt sequences computed from pSNVs function}

\item{prob.thresh}{Probability threshold of gains and losses. This value should be between 0.5 and 1.}

\item{log2.thresh}{Threshold for the absolute value of log ratio between wild type and mutant scores. Anything less than this value is discarded (default: 1).}

\item{include.cent}{If TRUE, gains and losses caused by mutation in the central STY residue are kept}
}
\description{
Score wt and mt sequences for a pwm
}
\keyword{internal}
\keyword{mut}
\keyword{psites}
\keyword{score}
\keyword{snp}
\keyword{snv}


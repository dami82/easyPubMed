% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epm_all_fx.R
\name{EPM_esearch_basic_q}
\alias{EPM_esearch_basic_q}
\title{Submit a Query to the NCBI ESearch Server.}
\usage{
EPM_esearch_basic_q(params)
}
\arguments{
\item{params}{List including the information
for querying the NCBI ESearch Server.}
}
\value{
Character vector including the response from the server.
}
\description{
Submit a Query to the NCBI ESearch Server and 
capture the response.
}
\details{
The \code{params} list must include the
elements listed below.
\itemize{
  \item `q`. String corresponding to the Query to be submitted to the server.
  \item `api_key`. (Optional) String corresponding to the NCBI API key.
}
}
\examples{
# Note: a time limit can be set in order to kill the operation when/if 
# the NCBI/Entrez server becomes unresponsive.
setTimeLimit(elapsed = 4.9)
try({
  my_q <- 'easyPubMed'
  my_params <- list(q = my_q)
  easyPubMed:::EPM_esearch_basic_q(params = my_params)
}, silent = TRUE)
setTimeLimit(elapsed = Inf)

}
\references{
\url{https://www.data-pulse.com/dev_site/easypubmed/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
\keyword{internal}

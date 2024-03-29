% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epm_all_fx.R
\name{EPM_efetch_basic_q}
\alias{EPM_efetch_basic_q}
\title{Submit a Query to the NCBI EFetch Server.}
\usage{
EPM_efetch_basic_q(params)
}
\arguments{
\item{params}{List including the information
for querying the NCBI EFetch Server.}
}
\value{
Character vector including the response from the server.
}
\description{
Submit a Query to the NCBI EFetch Server and 
capture the response.
}
\details{
The input list must include the 
elements listed below.
\itemize{
  \item `web_env`. String, unique value returned 
  by the NCBI ESearch server.
  \item `format`. String corresponding to the desired 
  response data format (e.g., "xml").
  \item `query_key`. Integer, key value returned by the 
  NCBI ESearch server.
  \item `retstart`. Integer, numeric index of the first 
  record to be request.
  \item `retmax`. Integer, maximum number of records to be retrieved 
  from the server.
  \item `encoding`. String, encoding of the data (e.g., "UTF-8").
}
}
\examples{
# Note: a time limit can be set in order to kill the operation when/if 
# the NCBI/Entrez server becomes unresponsive.
setTimeLimit(elapsed = 4.9)
try({
  x <- easyPubMed:::EPM_esearch_basic_q(params = list(q = "easyPubMed"))
  x <- easyPubMed:::EPM_esearch_parse(x)
  my_params <- list(web_env = x$web_env, 
                    query_key = x$query_key, 
                    format = "uilist")
  easyPubMed:::EPM_efetch_basic_q(params = my_params)
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

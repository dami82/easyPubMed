% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epm_surrogate_fx.R
\name{fetch_pubmed_data}
\alias{fetch_pubmed_data}
\title{Retrieve PubMed Data in XML or TXT Format}
\usage{
fetch_pubmed_data(
  pubmed_id_list,
  retstart = 0,
  retmax = 500,
  format = "xml",
  encoding = "UTF8",
  api_key = NULL,
  verbose = TRUE
)
}
\arguments{
\item{pubmed_id_list}{An easyPubMed object.}

\item{retstart}{Integer (>=0): this argument is ignored.}

\item{retmax}{Integer (>=1): this argument is ignored.}

\item{format}{String: element specifying the output format. The following values are allowed: 
c("xml", "medline", "uilist").}

\item{encoding}{String, the encoding of the records retrieved from Pubmed. 
This argument is ignored and set to 'UTF-8'.}

\item{api_key}{String, corresponding to the NCBI API token (if available). 
NCBI token strings can be requested from NCBI. Record download will be 
faster if a valid NCBI token is used. This argument can be NULL.}

\item{verbose}{Logical, shall details about the 
progress of the operation be printed to console.}
}
\value{
Character vector of length >= 1.
If format is set to "xml" (default), a single String including all 
PubMed records (decorated with XML tags) is returned. If a different format 
is selected, a vector of strings 
is returned, where each element corresponds to a line of the output document.
}
\description{
Retrieve PubMed records from Entrez following a search performed via the 
get_pubmed_ids() function. Data are downloaded in the XML or TXT format and are 
retrieved in batches of up to 5000 records.
}
\details{
The `fetch_pubmed_data()` function is now obsolete. 
You should use the `epm_fetch()` function instead.
Please, have a look at the manual or the vignette.
The `fetch_pubmed_data()` function will be retired in
the second half of 2024.
}
\examples{
## Example 01: retrieve PubMed record Unique Identifiers (uilist)
# Note: a time limit can be set in order to kill the operation when/if 
# the NCBI/Entrez server becomes unresponsive.
setTimeLimit(elapsed = 4.9)
try({ 
  q <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
  x <- get_pubmed_ids(pubmed_query_string = q)
  y <- fetch_pubmed_data(x, format = "uilist")
  y
}, silent = TRUE)
setTimeLimit(elapsed = Inf)

\dontrun{
## Example 02: retrieve data in XML format
q <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
x <- epm_query(query_string = q)
y <- fetch_pubmed_data(x, format = "xml")
y
}

}
\references{
\url{https://www.data-pulse.com/dev_site/easypubmed/}
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/}
}
\author{
Damiano Fantini \email{damiano.fantini@gmail.com}
}

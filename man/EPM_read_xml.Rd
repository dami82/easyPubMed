% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epm_all_fx.R
\name{EPM_read_xml}
\alias{EPM_read_xml}
\title{Import PubMed Records Saved Locally in XML Format.}
\usage{
EPM_read_xml(x)
}
\arguments{
\item{x}{Path to an XML file on the local machine.}
}
\value{
List including four elements: `guide` (data.frame), `meta` (list),      
`job_info` (data.frame) and `contents` (named list).
}
\description{
Read the contents of an XML file and import Metadata and
PubMed records for use by easyPubMed. 
The XML file must be generated by easyPubMed (ver >= 3) via 
the `epm_fetch()` function or via the `fetchEPMData()` method. 
XML files downloaded from the Web or using other software are 
currently unsupported. This function can only process one file.
}
\examples{
\dontrun{
  x <- epm_query(query_string = 'easyPubMed', verbose = TRUE)
  x <- epm_fetch(x = x, write_to_file = TRUE, store_contents = FALSE, 
                 outfile_prefix = 'qpm_qry_', verbose = TRUE)
  y <- EPM_read_xml(x = 'qpm_qry__batch_01.txt')
  try(unlink('qpm_qry__batch_01.txt'), silent = TRUE)
  y
}



}
\references{
\url{https://www.data-pulse.com/dev_site/easypubmed/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
\keyword{internal}

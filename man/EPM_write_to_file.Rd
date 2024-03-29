% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epm_all_fx.R
\name{EPM_write_to_file}
\alias{EPM_write_to_file}
\title{Write PubMed Records to Local Files.}
\usage{
EPM_write_to_file(x, to, format, addon = NULL, verbose = FALSE)
}
\arguments{
\item{x}{List including raw PubMed records.}

\item{to}{Path to the destination file on the local disc.}

\item{format}{String, format of the raw PubMed records that 
will be saved to the destination file (e.g., 'xml').}

\item{addon}{String, optional chunk of text in XML format to be written
to the destination file (header). This argument is only used when
`format` is set to 'xml'. It can be NULL.}

\item{verbose}{Logical, shall details about the 
progress of the operation be printed to console.}
}
\value{
integer in the range c(0, 1). A result of 0 indicates that an error
occurred while writing the file. A result of 1 indicates that the operation 
was completed successfully.
}
\description{
Write a list of PubMed records to a local file. If already existing, 
the destination file will be over-written. Original formatting of the
PubMed records should be declared and will be 
preserved in the output file. Format conversion is NOT supported.
}
\examples{
test <- list('Record #1', 'Record #2')
outfile = './test_file.txt'
file.exists(outfile)
easyPubMed:::EPM_write_to_file(x = test, to = './test_file.txt', format = 'xml')
file.exists(outfile)
readLines(outfile)
unlink(outfile)


}
\references{
\url{https://www.data-pulse.com/dev_site/easypubmed/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
\keyword{internal}

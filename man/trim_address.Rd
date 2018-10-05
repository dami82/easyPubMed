\name{trim_address}
\alias{trim_address}
\title{
Trim and Format Address Information
}
\description{
Set of rules for trimming and standardizing the format of address information retrieved from PubMed records. Affiliations including more than one address will be trimmend and only the first address will be returned. 
}
\usage{
trim_address(addr)
}
\arguments{
  \item{addr}{
Character string including an address as extracted from PubMed records.
}
}
\value{
Character string including a formatted and trimmed address (if available).
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
addr_string <- " 2 Dept of Urology, Feinberg School of Medicine," 
addr_string <- paste(addr_string, "Chicago, US; Dept of Mol Bio as well...")
addr_string
trim_address(addr = addr_string)
}

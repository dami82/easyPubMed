\name{table_articles_byAuth}
\alias{table_articles_byAuth}
\title{
Extract Publication and Affiliation Data from PubMed Records 
}
\description{
Extract Publication Info from PubMed records and cast data into a data.frame where each row corresponds to a different author. It is possible to retrieve data from first authors or last authors only as well as information from all authors of each PubMed record.
}
\usage{
table_articles_byAuth(pubmed_data, 
                      included_authors = "all", 
                      max_chars = 500, 
                      autofill = TRUE, 
                      dest_file = NULL, 
                      getKeywords = TRUE)
}
\arguments{
  \item{pubmed_data}{
PubMed Data in XML format: typically, an XML file resulting from a batch_pubmed_download() call or an XML object, result of a fetch_pubmed_data() call.
}
  \item{included_authors}{
Character: c("first", "last", "all"). Only includes information from the first, the last or all authors of a PubMed record.
}
  \item{max_chars}{
Numeric: maximum number of chars to extract from the AbstractText field.
}
  \item{autofill}{
Logical. If TRUE, missing affiliations are imputed according to the available values (from the same article).
}
  \item{dest_file}{
String (character of length 1). Name of the file that will be written for storing the output. If NULL, no file will be saved.
}
  \item{getKeywords}{
Logical. If TRUE, the operation will attempt to extract PubMed record keywords (MESH topics, keywords).
}
}
\details{
Retrieve publiccation and author information from PubMed data in the form of a data frame. 
}
\value{
Data frame including the following fields: c("article.title","article.abstract", "date.year","date.month","date.day", "journal.abbrv","journal.title", "keywords", "auth.last","auth.fore","auth.address", "auth.email").
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
\dontrun{
#
#
dami_query <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query)
dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed)
table_articles_byAuth(pubmed_data = dami_abstracts_xml, 
                      included_authors = "first", 
                      max_chars = 100, 
                      autofill = TRUE)[1:2,]
#
#
dami_query <- "Damiano Fantini[AU]"
curr.file <- batch_pubmed_download(dami_query, dest_file_prefix = "test_bpd_")
table_articles_byAuth(pubmed_data = curr.file[1], 
                      included_authors = "all", 
                      max_chars = 20, 
                      autofill =FALSE)
}
}

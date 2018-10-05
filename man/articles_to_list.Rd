\name{articles_to_list}
\alias{articles_to_list}
\title{
Cast PubMed Data into a List of Articles
}
\description{
Convert an XML object of PubMed records into a list of strings (character vector of length 1) corresponding to individual PubMed articles. PubMed records are identified by a "/PubmedArticle" XML tag. This automatically cast all the content of each PubMed record to a character-class object without removing XML tags.
}
\usage{
articles_to_list(pubmed_data)
}
\arguments{
  \item{pubmed_data}{
String corresponding to the name of an XML file (typically, the result of a batch_pubmed_download() call). Alternatively, an XML Object, such as the result of a fetch_pubmed_data() call.
}
}
\details{
The input is an XML object or an XML file, typically the result of a fetch_pubmed_data() call or a batch_pubmed_download() call. The function returns a list where each element is a different PubMed record.
}
\value{
List of character elements including the all records from the original XML object. Elements in the list are not named and are only accessible via the numeric index.
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
#
# retrieve PubMed data and return a list ot articles
dami_query <- "Damiano Fantini[AU]"
outfile <- batch_pubmed_download(dami_query, dest_file_prefix = "easyPM_ex001_")
listed_articles <- articles_to_list(pubmed_data = outfile)
listed_articles[[3]]
}

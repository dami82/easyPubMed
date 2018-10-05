\name{article_to_df}
\alias{article_to_df}
\title{
Extract Data from a PubMed Record
}
\description{
Extract publication-specific information from a PubMed record driven by XML tags. The input record is a string (character-class vector of length 1) and includes PubMed-specific XML tags. Data are returned as a data frame where each row corresponds to one of the authors of the PubMed article.
}
\usage{
article_to_df(pubmedArticle, autofill = FALSE, 
              max_chars = 500, getKeywords = FALSE, 
              getAuthors = TRUE)
}
\arguments{
  \item{pubmedArticle}{
String including one PubMed record.
}
  \item{autofill}{
Logical. If TRUE, missing affiliations are automatically imputed based on other non-NA addresses from the same record. 
}
  \item{max_chars}{
Numeric (integer). Maximum number of characters to be extracted from the Article Abstract field. Set max_chars to -1 for extracting the full-length abstract. Set max_chars to 0 to extract no abstract
}
  \item{getKeywords}{
Logical. If TRUE, an attempt to extract article Keywords will be made.
}
  \item{getAuthors}{
Logical. If FALSE, author information won't be extracted. This will considerably speed up the operation.
}
}
\details{
Given one Pubmed Article record, this function will automatically extract a set of features. Extracted information include: PMID, DOI, article title, article abstract, publication date (year, month, day), journal name (title, abbreviation), keywords, and a set of author-specific info (names, affiliation, email address). Each row of the output data frame corresponds to one of the authors of the PubMed record. Author-independent info (publication ID, title, journal, date) are identical across all rows. If information about authors are not required, set 'getAuthors' = TRUE. 
}
\value{
Data frame including the extracted features. Each row correspond a different author.
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
#
# Query PubMed, retrieve a selected citation and format it as a data frame
dami_query <- "Damiano Fantini[AU] AND 2017[PDAT]"
dami_on_pubmed <- get_pubmed_ids(dami_query)
dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed)
dami_abstracts_list <- articles_to_list(dami_abstracts_xml)
article_to_df(pubmedArticle = dami_abstracts_list[[1]], autofill = FALSE)
article_to_df(pubmedArticle = dami_abstracts_list[[2]], autofill = TRUE, max_chars = 300)[1:2,]
}

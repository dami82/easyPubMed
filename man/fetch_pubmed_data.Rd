\name{fetch_pubmed_data}
\alias{fetch_pubmed_data}
\title{
Retrieve PubMed Data in XML or TXT Format
}
\description{
Retrieve PubMed records from Entrez following a search performed via the get_pubmed_ids() function. Data are downloaded in the XML or TXT format and are retrieved in batches of up to 5000 records.
}
\usage{
fetch_pubmed_data(pubmed_id_list, retstart = 0, retmax = 500, 
                  api_key = NULL, format = "xml")
}
\arguments{
  \item{pubmed_id_list}{List: the result of a get_pubmed_ids() call.}
  \item{retstart}{Integer (>=0): index of the first UID in the retrieved PubMed Search Result set to be included in the output (default=0, corresponding to the first record of the entire set).}
  \item{retmax}{Integer (>=1): size of the batch of PubMed records to be retrieved at one time.}
  \item{api_key}{
is a string and corresponds to the personal API key linked to the NCBI/E-utilities account. If provided, E-services can queried at faster pace (10 queries per second). Can be NULL.
}
  \item{format}{Character: element specifying the output format. The following values are allowed: c("asn.1", "xml", "medline", "uilist", "abstract").}
}
\details{
Retrieve PubMed records based on the results of a get_pubmed_ids() query. Records are retrieved from Entrez via the PubMed API efetch function. The first entry to be retrieved may be adjusted via the retastart parameter (this allows the user to download large batches of PubMed data). The maximum number of entries to be retrieved can also be set adjusting the retmax parameter (1 < retmax < 5000). Data will be downloaded on the fly (no files are saved locally as a result of a fetch_pubmed_data() call).
}
\value{
If format == "xml": a XMLInternalDocument-class object. For accessing these data, use a XML parser.
If format != "xml": a "character" vector. Each element corresponds to one line of data.
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
\dontrun{
# Example 01: retrieve data in XML format
dami_query_string <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query_string)
dami_papers <- fetch_pubmed_data(dami_on_pubmed)
titles <- unlist(xpathApply(dami_papers, "//ArticleTitle", saveXML))
title_pos <- regexpr("<ArticleTitle>.*<\\\\/ArticleTitle>", titles)
titles <- substr(titles, title_pos + 14, title_pos + attributes(title_pos)$match.length - 16)
print(titles)
#
}
# Example 02: retrieve data in TXT format
dami_query_string <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query_string)
dami_papers <- fetch_pubmed_data(dami_on_pubmed, format = "abstract")
dami_papers[dami_papers == ""] <- "\n"
cat(paste(dami_papers[1:65], collapse = ""))
}

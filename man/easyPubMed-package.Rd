\name{easyPubMed-package}
\alias{easyPubMed-package}
\alias{easyPubMed}
\docType{package}
\title{
Retrieve and Process Scientific Publication Records from Pubmed
}
\description{
Query NCBI Entrez and retrieve PubMed records in XML or TXT format. PubMed records can be downloaded and saved as XML or text files. Data integrity is enforced during data download, allowing to retrieve and save very large number of records effortlessly. PubMed records can be processed to extract publication- and author-specific information.
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\references{
 \url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\examples{
# Example 01: retrieve data in XML format
dami_query_string <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query_string)
dami_papers <- fetch_pubmed_data(dami_on_pubmed)
titles <- unlist(xpathApply(dami_papers, "//ArticleTitle", saveXML))
title_pos <- regexpr('<ArticleTitle>.*<\\\\/ArticleTitle>', titles)
titles <- substr(titles, title_pos + 14, title_pos + attributes(title_pos)$match.length - 16)
print(titles)
#
# Example 02: retrieve data in TXT format
dami_query_string <- "Damiano Fantini[AU]"
dami_on_pubmed <- get_pubmed_ids(dami_query_string)
dami_papers <- fetch_pubmed_data(dami_on_pubmed, format = "abstract")
dami_papers[dami_papers == ""] <- "\n"
cat(paste(dami_papers[1:65], collapse = ""))
#
\dontrun{
# Example 03: retrieve data from PubMed and save as XML file
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out1 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
XML::xmlParse(out1[1])
#
# Example 04: retrieve data from PubMed and save as TXT file
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out2 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180, format = "medline")
readLines(out2[1])[1:30]
#
# Example 05: extract information from a single PubMed record 
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out3 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
PM_data <- articles_to_list(out3[1])
PM_record_df <- article_to_df(PM_data[[100]])
print(PM_record_df[1,])
print(PM_record_df[,"address"])
#
# Example 06: query PubMed and extract information from multiple records in one step 
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out4 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
PM_tab <- table_articles_byAuth(out4[1], autofill = TRUE, included_authors = "last")
PM_tab$address <- substr(PM_tab$address, 1, 15)
PM_tab[50:70,c("pmid", "jabbrv", "year", "lastname", "address")]
#
}
}

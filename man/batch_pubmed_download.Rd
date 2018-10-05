\name{batch_pubmed_download}
\alias{batch_pubmed_download}
\title{
Download PubMed Records in XML or TXT Format
}
\description{
Performs a PubMed Query (via the get_pubmed_ids() function), downloads the resulting data (via multiple fetch_pubmed_data() calls) and then saves data in a series of xml or txt files on the local drive. The function is suitable for downloading a very large number of records.
}
\usage{
batch_pubmed_download(pubmed_query_string, dest_dir = NULL, 
                      dest_file_prefix = "easyPubMed_data_", 
                      format = "xml", api_key = NULL, 
                      batch_size = 400, res_cn = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pubmed_query_string}{
String (character-vector of length 1): this is the string used for querying PubMed (the standard PubMed Query synthax applies).
}
  \item{dest_dir}{
String (character-vector of length 1): this string corresponds to the name of the existing folder where files will be saved. Existing files will be overwritten. If NULL, the current working directory will be used.
}
  \item{dest_file_prefix}{
String (character-vector of length 1): this string is used as prefix for the files that are written locally.
}
  \item{format}{
String (character-vector of length 1): data will be requested from Entrez in this format. Acceptable values are: c("medline","uilist","abstract","asn.1", "xml") When format != "xml", data will be saved as text files (txt).
}
  \item{api_key}{
is a string and corresponds to the personal API key linked to the NCBI/E-utilities account. If provided, E-services can queried at faster pace (10 queries per second). Can be NULL.
}
  \item{batch_size}{
Integer (1 < batch_size < 5000): maximum number of records to be saved in a single xml or txt file.
}
  \item{res_cn}{
Integer (> 0): numeric index of the data batch to start downloading from. This parameter is useful to resume an incomplete download job after a system crash.
}
}
\details{
Download large number of PubMed records as a set of xml or txt files that are saved in the folder specified by the user. This function enforces data integrity. If a batch of downloaded data is corrupted, it is discarded and downloaded again. Each download cycle is monitored until the download job is successfully completed. This function should allow to download a whole copy of PubMed, if desired. The function informs the user about the current progress by constantly printing to console the number of batches still in queue for download. pubmed_query_string accepts standard PubMed synthax. The function will query PubMed multiple times using the same query string. Therefore, it is recommended to use a [EDAT] or a [PDAT] filter in the query if you want to ensure reproducible results.
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
\dontrun{
# Example 01: retrieve data from PubMed and save as XML file
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out1 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
XML::xmlParse(out1[1])
#
# Example 02: retrieve data from PubMed and save as TXT file
ml_query <- "Machine Learning[TI] AND 2016[PD]"
out2 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180, format = "medline")
readLines(out2[1])[1:30]
}
}


#  easyPubMed, ver 3.1.3
#  Retrieve and Process Scientific Publication Records from Pubmed

#  Copyright (C) 2024, Damiano Fantini

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.



# ---



#' @title Preprocessed PubMed Records and Data
#'
#' @description This dataset includes a collection of sample 
#' data obtained from PubMed records and saved in different formats.  
#' This dataset is used to demonstrate specific functionalities of 
#' the `easyPubMed` R library. Each element in the `epm_samples` list 
#' corresponds to a different input or intermediate object. 
#'
#' @usage data("epm_samples")
#' 
#' @format The dataset is formatted as a list including 4 elements: 
#' 
#' * `bladder_cancer_2018`: List of 4
#' 
#' * `bladder_cancer_40y`: List of 1
#' 
#' * `fx`: List of 5
#' 
#'   
#'
#' @examples 
#' ## Display some contents
#' data("epm_samples")
#' # Display Query String used for collecting the data
#' print(epm_samples$bladder_cancer_2018$demo_data_01)
"epm_samples"


#' @title PubMed Query Stopwords 
#'
#' @description Collection of 133 Stopwords that can be removed from query strings to improve the 
#'              accuracy of exact-match PubMed queries.
#'
#' @usage data("epm_stopwords")
#' 
#' @format A character vector including all PubMed stopwords tat are typically filtered out from queries. 
#'   
#' @details Number of stopwords included, n=133.
#'
#' @examples 
#' ## Display some contents
#' data("epm_stopwords")
#' head(epm_stopwords)
"epm_stopwords"



#' @title Retrieve and Process Scientific Publication Records from Pubmed
#' 
#' @description Query NCBI Entrez and retrieve PubMed records in XML or TXT format. PubMed records 
#' can be downloaded and saved as XML or text files. Data integrity is enforced during data download, 
#' allowing to retrieve and save very large number of records effortlessly. PubMed records can be processed 
#' to extract publication- and author-specific information.
#' 
#' @details 
#' This software is based on the information 
#' included in the Entrez Programming Utilities Help manual 
#' authored by Eric Sayers, PhD and available on the NCBI Bookshelf (NBK25500).
#' This R library is NOT endorsed, supported, maintained NOR affiliated with NCBI.
#' 
#' 
#' 
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references 
#' \itemize{
#'   \item Tutorials and Help Webpage: \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'   \item NCBI PubMed Help Manual: \url{https://pubmed.ncbi.nlm.nih.gov/help/}
#'   \item Entrez Programming Utilities Help (NBK25500): \url{https://www.ncbi.nlm.nih.gov/books/NBK25500/}
#' }
#' 
#' 
#' @examples
#' ## Example 01: retrieve data in XML format, extract info, show
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   my_query_string <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   epm <- epm_query(my_query_string)
#'   epm <- epm_fetch(epm)
#'   epm <- epm_parse(epm, max_authors = 5, max_references = 10)
#'   processed_data <- get_epm_data(epm)
#'   utils::head(processed_data)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' \dontrun{
#' ## Example 02: retrieve data in medline format
#' my_query_string <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#' epm <- epm_query(my_query_string)
#' epm <- epm_fetch(epm, format = 'medline')
#' medline_data <- get_epm_raw(epm)
#' first_record <- medline_data[[1]] 
#' cat(first_record, sep = '\n')
#' 
#' 
#' ## Additional Examples: show easyPubMed Vignette
#' library(easyPubMed)
#' vignette("easyPubMed_demo")
#' 
#' }
#' 
#' @keywords internal
"_PACKAGE"


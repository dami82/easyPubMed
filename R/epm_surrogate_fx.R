
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



## ---------------------------
## --- Surrogate Functions ---
## ---------------------------


#' Simple PubMed Record Search
#' 
#' Query PubMed (Entrez) in a simple way via the PubMed API eSearch function. 
#' Calling this function results in posting the query results on the PubMed 
#' History Server. This allows later access to the resulting data via the 
#' fetch_pubmed_data() function, or other easyPubMed functions. 
#' NOTE: this function has become obsolete. You should use the epm_query() 
#' function instead. Please, have a look at the manual or the vignette.
#' The \code{get_pubmed_ids()} function will be retired in 2024.
#' 
#' 
#' 
#' @details 
#' This function will use the String provided as argument for 
#' querying PubMed via the eSearch function of the PubMed API. 
#' The Query Term can include one or multiple words, as well as the 
#' standard PubMed operators (AND, OR, NOT) and 
#' tags (i.e., [AU], [PDAT], [Affiliation], and so on). ESearch will post 
#' the UIDs resulting from the search operation onto the History server 
#' so that they can be used directly in a subsequent fetchPubmedData() call. 
#' 
#' 
#' 
#'
#' @param pubmed_query_string String (character vector of length 1), 
#' corresponding to the query string used for querying PubMed.
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be NULL.
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   qry <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   get_pubmed_ids(pubmed_query_string = qry)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return An easyPubMed object which includes no PubMed records.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'  
#' @importFrom rlang warn
#'
#' 
#' @export
get_pubmed_ids <- function(pubmed_query_string, api_key = NULL) {
  
  warn_msg <- paste(
    'The get_pubmed_ids() function has become', 
    'obsolete. You should use the EPM_query() function instead.',
    'Please, have a look at the manual or the vignette.',
    'The get_pubmed_ids() function will be retired in the second half of 2024.')
  
  rlang::warn(message = warn_msg, .frequency = 'once', 
              .frequency_id='epm_wrn001_gpi')  
  
  y <- epm_query(query_string = pubmed_query_string, 
                 api_key = api_key, verbose = FALSE)
  
  return(y)
}




#' @title Retrieve PubMed Data in XML or TXT Format
#'
#' @description Retrieve PubMed records from Entrez following a search performed via the 
#' get_pubmed_ids() function. Data are downloaded in the XML or TXT format and are 
#' retrieved in batches of up to 5000 records.
#' 
#' 
#' @param pubmed_id_list An easyPubMed object.
#' @param retstart Integer (>=0): this argument is ignored.
#' @param retmax Integer (>=1): this argument is ignored.
#' @param format String: element specifying the output format. The following values are allowed: 
#' c("xml", "medline", "uilist").
#' @param encoding String, the encoding of the records retrieved from Pubmed. 
#' This argument is ignored and set to 'UTF-8'.
#' @param api_key String, corresponding to the NCBI API token (if available). 
#' NCBI token strings can be requested from NCBI. Record download will be 
#' faster if a valid NCBI token is used. This argument can be NULL. 
#' @param verbose Logical, shall details about the 
#' progress of the operation be printed to console.
#' 
#' 
#' @details 
#' The `fetch_pubmed_data()` function is now obsolete. 
#' You should use the `epm_fetch()` function instead.
#' Please, have a look at the manual or the vignette.
#' The `fetch_pubmed_data()` function will be retired in
#' the second half of 2024.
#' 
#' @return 
#' Character vector of length >= 1.
#' If format is set to "xml" (default), a single String including all 
#' PubMed records (decorated with XML tags) is returned. If a different format 
#' is selected, a vector of strings 
#' is returned, where each element corresponds to a line of the output document.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/}
#'
#' @examples
#' ## Example 01: retrieve PubMed record Unique Identifiers (uilist)
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({ 
#'   q <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   x <- get_pubmed_ids(pubmed_query_string = q)
#'   y <- fetch_pubmed_data(x, format = "uilist")
#'   y
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' \dontrun{
#' ## Example 02: retrieve data in XML format
#' q <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#' x <- epm_query(query_string = q)
#' y <- fetch_pubmed_data(x, format = "xml")
#' y
#' }
#' 
#' @importFrom utils head
#' @export
fetch_pubmed_data <- function (pubmed_id_list,
                               retstart = 0,
                               retmax = 500,
                               format = "xml", 
                               encoding = "UTF8", 
                               api_key = NULL, 
                               verbose = TRUE) 
{
  
  warn_msg <- paste(
    'The fetch_pubmed_data() function has become', 
    'obsolete. You should use the epm_fetch() function instead.',
    'Please, have a look at the manual or the vignette.',
    'The fetch_pubmed_data() function will be retired in', 
    'the second half of 2024.')
  
  rlang::warn(message = warn_msg, .frequency = 'once', 
              .frequency_id='epm_wrn002_gpi') 
  
  # Ignore these args
  zz <- retstart; zz <- retmax; zz <- encoding
  
  # check other args
  stopifnot(inherits(pubmed_id_list, 'easyPubMed'), 
            is.character(format), length(format)==1, !is.na(format))
  
  
  # Fetch the data
  x <- epm_fetch(x = pubmed_id_list, format = format, 
                 api_key = api_key, write_to_file = FALSE, 
                 outfile_path = NULL, outfile_prefix = NULL, 
                 store_contents = TRUE, encoding = 'UTF-8', 
                 verbose = verbose)
  x <- getEPMRaw(x)
  
  # take care of XML records
  if (format == 'xml') {
    li <- list(
      '<?xml version="1.0" ?>', 
      '<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, ', 
      '1st January 2023//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/', 
      'out/pubmed_230101.dtd">', 
      '<PubmedArticleSet>')
    
    # write each element of x
    for (j in seq(1, length(x), by = 1)) {
      if (nchar(x[[j]]) > 0) {
        li[[length(li) + 1]] <- 
          tryCatch({paste0('<PubmedArticle>', x[[j]], '</PubmedArticle>') }, 
                   error = function(e) { NULL })
      }
    }
    li[[length(li) + 1]] <- '</PubmedArticleSet>'
    x <- do.call(paste0, li)
    
    
  } else if (format == "medline") {
  
    x <- lapply(x, function(i) { c(i, ' ', ' ', ' ')})
    names(x) <- NULL
    x <- do.call(c, x)
    
    
  } else if (format == 'uilist') {
    
    x <- do.call(c, x)
    names(x) <- NULL
  }
  
  return(x)
}








#' @title Extract Publication and Affiliation Data from PubMed Records
#'
#' @description Extract Publication Info from PubMed records and cast data into a 
#' data.frame where each row corresponds to a different author. It is possible to limit
#' data extraction to first authors or last authors only, or get information about 
#' all authors of each PubMed record.
#' 
#' 
#' @param pubmed_data PubMed Data in XML format: typically, an XML file resulting from a 
#' batch_pubmed_download() call or an XML object, result of a fetch_pubmed_data() call.
#' @param included_authors Character: c("first", "last", "all"). Only includes information 
#' from the first, the last or all authors of a PubMed record.
#' @param max_chars This argument is ignored. In this version of the function, 
#' the whole Abstract Text is returned.
#' @param autofill Logical. If TRUE, missing affiliations are imputed according to the available 
#' values (from the same article).
#' @param dest_file String (character of length 1). Name of the file that will be written for 
#' storing the output. If NULL, no file will be saved.
#' @param getKeywords This argument is ignored. In this version of the function 
#' MeSH terms and codes (i.e., keywords) are parsed by default.
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the function 
#' base::iconv(). See iconv() help page for how to find out more about encodings that can be 
#' used on your platform. Here, we recommend using "UTF-8".
#' 
#' @details 
#' The `table_articles_byAuth()` function is now obsolete. 
#' You should use the `epm_parse()` function instead.
#' Please, have a look at the manual or the vignette.
#' The `table_articles_byAuth()` function will be retired in
#' the second half of 2024.
#' 
#' @return 
#' Data frame including the following fields: `c("pmid", "doi", 
#' "title", "abstract", "year", "month", "day", "jabbrv", "journal", 
#' "keywords", "mesh", "lastname", "firstname", "address", "email")`.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   q0 <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   q1 <- easyPubMed::get_pubmed_ids(pubmed_query_string = q0)
#'   q2 <- fetch_pubmed_data(pubmed_id_list = q1)
#'   df <- table_articles_byAuth(q2, included_authors = 'first')
#'   df[, c('pmid', 'lastname', 'jabbrv', 'year', 'month', 'day')]
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' 
#' @importFrom rlang warn
#' @importFrom utils write.table
#' 
#' @export
table_articles_byAuth <- function (pubmed_data, 
                                   included_authors = "all", 
                                   max_chars = 500, 
                                   autofill = TRUE, 
                                   dest_file = NULL, 
                                   getKeywords = TRUE, 
                                   encoding = "UTF8") {
  
  # inner import f(x)
  inner_fx <- function(x) {
    
    pmart_pat <- '<PubmedArticle(>|([[:space:]]+[^>]*>))'
    y <- tryCatch({strsplit(x, split = pmart_pat)[[1]][-1]}, 
                  error = function(e) { NULL })
    y <- tryCatch({lapply(y, function(k) {
      tryCatch({sub('</PubmedArticle>.*$', '', k)}, 
               error = function(e) { NULL })})}, 
      error = function(e) { NULL })
    
    return(y)
  }
  
  # These are ignored
  zz <- max_chars
  zz <- getKeywords
  zz <- NULL
  
  # Warn that this f(x) is becoming obsolete
  warn_msg <- paste(
    'The fetch_pubmed_data() function has become', 
    'obsolete. You should use the epm_fetch() function instead.',
    'Please, have a look at the manual or the vignette.',
    'The fetch_pubmed_data() function will be retired in', 
    'the second half of 2024.')
  
  rlang::warn(message = warn_msg, .frequency = 'once', 
              .frequency_id='epm_wrn003_gpi') 
  
  stopifnot(is.character(included_authors), length(included_authors) == 1, 
            !is.na(included_authors), 
            included_authors %in% c("all", "first", "last")) 

  incl_auth <- c('all'=1000, 'first'=1, 'last'=-1)[included_authors]
  
  # init collector
  x_raw <- list()
  
  if (is.character(pubmed_data) && length(pubmed_data) == 1 && 
      !is.na(pubmed_data) && grepl('^<\\?xml', substr(pubmed_data, 1, 20))) {
    
    # read from single string
    x_raw <- inner_fx(pubmed_data)

  } else if (is.character(pubmed_data) && length(pubmed_data) > 0 && 
             sum(file.exists(pubmed_data)) > 0) {
    
    pubmed_data <- pubmed_data[file.exists(pubmed_data)]
    pubmed_data <- unique(pubmed_data)

    for (fi in pubmed_data) {
      
      TMP_i <- readLines(fi)
      TMP_i <- paste(TMP_i, collapse = ' ')
      TMP_i <- inner_fx(TMP_i)
      
      for (ti in TMP_i) {
        x_raw[[length(x_raw) + 1]] <- ti
      }
    }
  } else {
    stop('Unsupported Input Type.')
  }

    # prep for messaging (if requested)            
  if (length(x_raw) > 20) {
    prc_seq <- unique(as.integer(
      seq(1, length(x_raw), length.out = 19)))
  } else {
    prc_seq <- length(x_raw)
  }
  
  if (TRUE) {
    vrb_li1 <- '|----+----+----+----| 100%'
    vrb_li2 <- '...................|'
    message(vrb_li1)
    message('|', appendLF = FALSE)
  }
  
  # collector (init) and loop
  y <- list()
  for (i in seq_len(length(x_raw))) {
    
    # message
    if (TRUE && length(prc_seq) > 1 && i %in% prc_seq){
      message('.', appendLF = FALSE)  
    }
    
    # parse
    y[[length(y) + 1 ]] <- tryCatch({
      epm_parse_record(x_raw[[i]], 
                       max_authors = incl_auth, 
                       autofill_address = autofill, 
                       compact_output = FALSE,
                       include_abstract = TRUE,
                       max_references = 0,
                       ref_id_type = 'doi')
    }, error = function(e) { NULL })
  }
  
  # message
  if (TRUE) {
    if (length(prc_seq) > 1) {
      message('|', appendLF = TRUE)  
    } else {
      message(vrb_li2, appendLF = TRUE)
    }
  }
  
  # collapse
  tot_processed_n <- length(y)
  y <- tryCatch({do.call(rbind, y)}, error = function(e) { NULL })
  stopifnot(rownames(y) > 0)
  rownames(y) <- NULL
  
  yy <- data.frame(pmid=y$pmid, doi=y$doi, title=y$title, abstract=y$abstract, 
                   year=y$year, month=y$month, day=y$day, 
                   jabbrv=y$jabbrv, journal=y$journal, 
                   keywords=y$mesh_terms, mesh=y$mesh_codes, 
                   lastname=y$last_name, firstname=y$first_name, 
                   address=y$affiliation, email=y$email, 
                   stringsAsFactors = FALSE)
  rownames(yy) <- NULL
  
  # shall we write to file?
  if (!is.null(dest_file)) {
    if (is.character(dest_file) && length(dest_file) == 1) {
      tryCatch(utils::write.table(yy, dest_file, fileEncoding = encoding), 
               error = function(e) {
                 NULL
               })
    }
  }
  
  # return
  return(yy)
}


  
 



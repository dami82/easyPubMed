
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



## ----------------
## --- Internal ---
## ----------------


#' Custom XML Tag Matching.
#'
#' Extract text form a string containing XML or HTML tags. 
#' Text included between tags of interest will be returned. 
#' If multiple tagged substrings are found, they will be returned as different 
#' elements of a list or character vector.
#'
#' @param xml_data String (character vector of length 1), this is a string
#' including PubMed records or string including XML/HTML tags.
#' 
#' @param tag String (character vector of length 1), the tag of 
#' interest (e.g., "Title") (should NOT include < > chars).
#' 
#' @param xclass String (character vector of length 1), a tag decorator 
#' of interest (e.g., "EIdType=\"doi\""). Can be NULL.
#' 
#' @param format String. Must be a value in c("list", "char"). Indicates the 
#' type of output. Defaults to "list".
#' 
#' @examples 
#' x <- "This string includes <Ti>an XML Tag</Ti>."
#' easyPubMed:::EPM_custom_grep(x, tag = "Ti")
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return List or vector where each element corresponds to an in-tag substring.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#' @keywords internal
EPM_custom_grep <- function (xml_data, tag, xclass = NULL, format = "list") 
{
  x <- xml_data[[1]]
  
  if (is.null(xclass)) {
    #tag.op <- paste0("\\<", tag, "((\\>)|([[:space:]]([^[<]]*)\\>))")
    tag.op <- paste0("\\<", tag, "((\\>)|([[:space:]]([^[<]]*)>))")
  } else {
    tag.op <- paste0("\\<", tag, "[[:space:]]([^>])*(", xclass, ")([^>]*)>")
  }
 
  tag.cl <- paste("(<\\/)", tag, "(\\>)", sep = "")
  out.result <- list()
  i = 1
  while (!is.null(x) && !is.na(x) && x != "" && 
         nchar(x) > 0 && regexpr(tag.op, x) > 0 && 
         regexpr(tag.cl, x) > 0) {
    
    tag.op.pos <- regexpr(tag.op, x)
    nu.x <- substr(x, (tag.op.pos - 1), nchar(x))
    inner.trim <- regexpr(">", nu.x, fixed = TRUE)
    nu.x <- substr(nu.x, (inner.trim + 1), nchar(nu.x))
    tag.cl.pos <- regexpr(tag.cl, nu.x)
    tag.cl.full <- tag.cl.pos + attributes(tag.cl.pos)$match.length +  1
    x <- substr(nu.x, tag.cl.full, nchar(x))
    nu.x <- substr(nu.x, 1, (tag.cl.pos - 1))
    out.result[[i]] <- nu.x
    i <- i + 1
  }
  if (format != "list") {
    out.result <- do.call(c, out.result)
  }
  return(out.result)
}


#' Submit a Query and Read the Response from the Server.
#' 
#' Submit a request to a server (typically, the Entrez Eutils server) and 
#' capture the response. 
#'
#' @param qurl String (character vector of length 1), corresponding to the 
#' query URL to the remote server.
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   qry <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
#'                 "esearch.fcgi?db=pubmed&term=easyPubMed")
#'   easyPubMed:::EPM_submit_q(qry)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Character vector including the response from the server.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' @importFrom utils head
#' 
#' @keywords internal
EPM_submit_q <- function(qurl) { 
  
  stopifnot(is.character(qurl), 
            length(qurl) == 1, 
            nchar(qurl) > 0)
  
  # establish connection
  IDconnect <- tryCatch({ suppressWarnings(
    url(qurl, open = "rb", encoding = "UTF8"))}, 
    error = function(e) { NULL })
  
  y <- tryCatch({
    
    # read remote content    
    idXML <- suppressWarnings(
      readLines(IDconnect, warn = FALSE, encoding = "UTF8"))
    
    # If error, return NULL
    # Default behavior in case of an error -> TRUE (return NULL)
    CHK1 <- tryCatch({
      tmp_chk <- lapply(utils::head(idXML), substr, start = 1, stop = 50)
      tmp_chk <- paste(do.call(c, tmp_chk), collapse = '')
      as.logical(grepl("<ERROR>", tmp_chk))
    }, error = function(e) { TRUE })
    
    if (CHK1) { idXML <- NULL }
    
    # Output
    idXML
    
  }, error = function(e) {
    NULL
  }, finally = {
    try(suppressWarnings(close(IDconnect)), silent = TRUE)
  })
  
  return(y)
}


#' Submit a Query to the NCBI ESearch Server.
#' 
#' Submit a Query to the NCBI ESearch Server and 
#' capture the response. 
#'
#' @param params List including the information
#' for querying the NCBI ESearch Server. 
#' 
#' @details 
#' The \code{params} list must include the
#' elements listed below.
#' \itemize{
#'   \item `q`. String corresponding to the Query to be submitted to the server.
#'   \item `api_key`. (Optional) String corresponding to the NCBI API key.
#' }
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   my_q <- 'easyPubMed'
#'   my_params <- list(q = my_q)
#'   easyPubMed:::EPM_esearch_basic_q(params = my_params)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Character vector including the response from the server.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_esearch_basic_q <- function(params) {
  
  # Check
  stopifnot(is.list(params), length(params) > 0,
            !is.null(params$q), is.character(params$q), 
            length(params$q) == 1, nchar(params$q) > 0)
  
  # remove spaces from query
  q <- tryCatch({gsub('[[:space:]]', '+', params$q)}, 
                error = function(e) { params$q })
  
  myPubmedURL <- paste(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
    "db=pubmed&term=", q, "&usehistory=y", sep = "")
  
  api_key <- tryCatch({params$api_key}, error = function(e) { NULL })
  
  if (!is.null(api_key)) {
    if (is.character(api_key) && length(api_key) == 1 && nchar(api_key) > 0) {
      myPubmedURL <- paste(myPubmedURL, "&api_key=", api_key, sep = "")
    }
  }
  
  # Exec, collapse and return
  y <- tryCatch(EPM_submit_q(myPubmedURL), error = function(e) { NULL })
  y <- tryCatch({paste(y, collapse = '')}, error = function(e) { NULL })
  return(y)
}



#' Parse Responses from the NCBI ESearch Server.
#' 
#' Parse Responses from the NCBI ESearch Server and 
#' return a list of information that can be used for retrieving 
#' PubMed records from the NCBI EFetch Server. 
#'
#' @param x String (character vector of length 1), this is the 
#' xml string returned by the NCBI ESearch Server. 
#' 
#' 
#' @examples
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   my_q <- 'easyPubMed'
#'   my_params <- list(q = my_q)
#'   x <- easyPubMed:::EPM_esearch_basic_q(params = my_params)
#'   easyPubMed:::EPM_esearch_parse(x)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return List including information extracted from the NCBI ESearch 
#' Server response. 
#' 
#' @details 
#' The output list includes the following items.
#' \itemize{
#'   \item `web_env`. String, unique identifier for fetching PubMed records 
#'   corresponding to the current query.
#'   \item `query_key`. Integer, unique numeric key for fetching PubMed records
#'   corresponding to the current query.
#'   \item `count`. Integer, expected number of records returned 
#'   by the current query.
#'   \item `query_translation`. String, translation of the Query string 
#'   provided by the user.
#' }
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_esearch_parse <- function(x) {
  
  # check x
  stopifnot(!is.null(x), is.character(x), 
            length(x) == 1, !is.na(x))
  
  # init collector
  y <- list()
  
  # attempt extraction
  y$web_env <- tryCatch({
    suppressWarnings(
      as.character(
        EPM_custom_grep(xml_data = x, tag = 'WebEnv', format = 'char')))
  }, error = function(e) { NA })
  
  y$query_key <- tryCatch({
    suppressWarnings(
      as.numeric(
        EPM_custom_grep(xml_data = x, tag = 'QueryKey', format = 'char')))
  }, error = function(e) { NA })
  
  y$count <- tryCatch({
    suppressWarnings(
      as.numeric(
        EPM_custom_grep(xml_data = x, tag = 'Count', format = 'char')))
  }, error = function(e) { NA })
  
  y$query_translation <- tryCatch({
    suppressWarnings(
      as.character(
        EPM_custom_grep(xml_data = x, tag = 'QueryTranslation', 
                        format = 'char')))
  }, error = function(e) { NA })
  
  return(y)
}


#' Submit a Query to the NCBI EFetch Server.
#' 
#' Submit a Query to the NCBI EFetch Server and 
#' capture the response. 
#'
#' @param params List including the information
#' for querying the NCBI EFetch Server. 
#' 
#' @details The input list must include the 
#' elements listed below.
#' \itemize{
#'   \item `web_env`. String, unique value returned 
#'   by the NCBI ESearch server.
#'   \item `format`. String corresponding to the desired 
#'   response data format (e.g., "xml").
#'   \item `query_key`. Integer, key value returned by the 
#'   NCBI ESearch server.
#'   \item `retstart`. Integer, numeric index of the first 
#'   record to be request.
#'   \item `retmax`. Integer, maximum number of records to be retrieved 
#'   from the server.
#'   \item `encoding`. String, encoding of the data (e.g., "UTF-8").
#' }
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- easyPubMed:::EPM_esearch_basic_q(params = list(q = "easyPubMed"))
#'   x <- easyPubMed:::EPM_esearch_parse(x)
#'   my_params <- list(web_env = x$web_env, 
#'                     query_key = x$query_key, 
#'                     format = "uilist")
#'   easyPubMed:::EPM_efetch_basic_q(params = my_params)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Character vector including the response from the server.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_efetch_basic_q <- function(params) {
  
  # check
  stopifnot(is.list(params), length(params) > 0, 
            !is.null(params$web_env), 
            is.character(params$web_env), 
            length((params$web_env)) == 1, 
            nchar(params$web_env) > 0)
  
  # retrieve params
  enco0 <- tryCatch({params$encoding}, error = function(e) {'UTF-8'})
  wenv0 <- tryCatch({params$web_env}, error = function(e) { NULL })
  qkey0 <- tryCatch({params$query_key}, error = function(e) { NULL })
  rsta0 <- tryCatch({params$retstart}, error = function(e) { 0 })
  rmax0 <- tryCatch({params$retmax}, error = function(e) { 500 })
  frmt0 <- tryCatch({params$format}, error = function(e) {'xml'})
  
  # Retrieval params  
  ret_li <- list(
    "xml" = list(ret_type = 'null', ret_mode = 'xml'), 
    "medline" = list(ret_type = 'medline', ret_mode = 'text'), 
    "uilist" = list(ret_type = 'uilist', ret_mode = 'text'))
  
  # Compose the query string
  efetch_url <- tryCatch({
    paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", 
          "db=pubmed&WebEnv=", wenv0, "&query_key=", qkey0, 
          "&retstart=", rsta0, "&retmax=", rmax0, 
          "&rettype=", ret_li[[frmt0]]$ret_type, 
          "&retmode=", ret_li[[frmt0]]$ret_mode, sep = "")
  }, error = function (e) {NULL}) 
  
  api_key <- params$api_key
  if (!is.null(api_key)) {
    efetch_url <- tryCatch({
      paste(efetch_url, "&api_key=", api_key, sep = "")
    }, error = function(e) { NULL })
  }
  
  # Exec and return
  # No collapsing nor handling
  y <- tryCatch(EPM_submit_q(efetch_url), error = function(e) { NULL })
  
  # Deal with conditional collapsing & splitting
  if (!is.null(y)) {
    
    if (params$format == 'xml') {
      y <- tryCatch(paste(y, collapse = ''), error = function(e) { NULL })
      
      # Split
      # pmart_pat <- '<PubmedArticle(>|([[:space:]]+[^>]*>))'
      pmart_pat <- paste0('(<PubmedArticle(>|([[:space:]]+[^>]*>)))|', 
                          '(<PubmedBookArticle(>|([[:space:]]+[^>]*>)))')

      y <- tryCatch({strsplit(y, split = pmart_pat)[[1]][-1]}, 
                     error = function(e) { NULL })

      # y <- tryCatch({lapply(y, function(k) {
      #   tryCatch({sub('</PubmedArticle>.*$', '', k)}, 
      #            error = function(e) { NULL })})}, 
      #   error = function(e) { NULL })
      y <- tryCatch({lapply(y, function(k) {
        tryCatch({sub('</PubmedBookArticle>.*$', '', sub('</PubmedArticle>.*$', '', k))}, 
                 error = function(e) { NULL })})}, 
        error = function(e) { NULL })
      
    } else if (params$format == 'medline') {
      y <- tryCatch({
        idx1 <- which(nchar(y) == 0) + 1
        idx2 <- which(grepl(pattern = '^PMID', y))
        idx1 <- idx1[idx1 %in% idx2]
        ct_guide <- data.frame (init = idx1, 
                                stop = c((idx1[-1])-1, length(y)))
        
        lapply(seq(1, nrow(ct_guide), by=1), function(i) {
          keep_i <- seq(ct_guide$init[i], ct_guide$stop[i], by = 1)
          tmp <- y[keep_i]
          tmp <- tryCatch({tmp[nchar(tmp) > 0]}, error = function(e) { tmp })
          tmp
        })
      }, error = function(e) { NULL })
    } else if (params$format == 'uilist') {
      y <- lapply(y, function(i) { i })
      
    }
  }
  
  # Return
  return(y)
}


#' Retrieve Results via an Esearch and Efetch sequence.
#' 
#' Submit a Query to the NCBI ESearch Server, 
#' capture the response and retrieve the corresponding PubMed records from 
#' the NCBI EFetch Server. 
#' Up to the first n=10,000 records returned by the query 
#' will be retrieved (as per the NCBI policy). This does not include
#' a timeout limit to complete the operation.
#' 
#' 
#'
#' @param query_string String (character vector of length 1), 
#' corresponding to the query URL to the remote server.
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be NULL.
#' @param batch_size Integer, max number of records to be retrieved 
#' as a batch. This corresponds to the "retmax" NCBI parameter.
#' @param encoding String (character vector of length 1), encoding
#' of the resulting records (e.g., "UTF-8").
#' @param format String (character vector of length 1),
#' desired format of the Pubmed records. This must be one of the values in
#' c("xml", "medline", "uilist"). 
#' @param max_restart_attempts Integer, max number of attempts in case of a
#' failed iteration.
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   qry <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   easyPubMed:::EPM_esearch_efetch_seq(query_string = qry, format = "uilist")
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Character vector including the response from the server.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_esearch_efetch_seq <- function(query_string, api_key = NULL, 
                                   batch_size = 500, encoding = 'UTF-8', 
                                   format = 'xml', max_restart_attempts = 10) {
  
  # check params
  stopifnot(!is.null(query_string), is.character(query_string), 
            length(query_string) == 1, !is.na(query_string), 
            nchar(query_string) > 0, 
            !is.null(batch_size), is.numeric(batch_size), 
            length(batch_size) == 1, !is.na(batch_size), 
            !is.null(max_restart_attempts), is.numeric(max_restart_attempts), 
            length(max_restart_attempts) == 1, !is.na(max_restart_attempts), 
            !is.null(format), is.character(format), length(format) == 1, 
            !is.null(encoding), is.character(encoding), length(encoding) == 1)
  
  if (!is.null(api_key)) {
    stopifnot(is.character(api_key), 
              length(api_key) == 1, 
              !is.na(api_key), 
              nchar(api_key) > 0)
    
    # pace of query (secs per query) 
    rltime <- 0.15
  } else {
    
    # pace of query w/o API KEY (secs per query)
    rltime <- 0.35
  }
  
  # Handle format arg
  format <- match.arg(arg = format, 
                      choices = c('xml', 'medline', 'uilist'), 
                      several.ok = FALSE)
  
  # Handle params
  batch_size <- round(batch_size)
  batch_size[batch_size < 5] <- 5
  batch_size[batch_size > 2000] <- 2000
  
  max_restart_attempts <- round(max_restart_attempts)
  max_restart_attempts[max_restart_attempts < 5] <- 5
  max_restart_attempts[max_restart_attempts > 100] <- 100
  
  # Init, collect expectations
  q_list <- list(q=query_string, api_key=api_key)
  q_0 <- EPM_esearch_basic_q(q_list)
  q_0 <- EPM_esearch_parse(q_0)
  
  # Limit requests rate
  Sys.sleep(rltime)
  
  if (q_0$count > 10000) {
    w_msg <- paste0('This query returned a large number of records (n=', 
                    q_0$count, '). Only the first n=10,000 records ', 
                    'will be retrieved (as per the NCBI Entrez policy).')
    warning(w_msg)
  }
  
  # prepare guide
  rec_num <- min(c(q_0$count, 10000))
  tmp_retstarts <- seq(0, rec_num, by = batch_size)
  my_guide <- data.frame(idx = seq_along(tmp_retstarts), 
                         retstart = tmp_retstarts, 
                         retmax = batch_size, 
                         completed = FALSE, 
                         stringsAsFactors = FALSE)
  
  # Initialize collector and monitor 
  out <- list()
  attempt_n <- 1
  check_t0 <- Sys.time()
  iter_success <- FALSE
  
  while(sum(!my_guide$completed) > 0 && attempt_n <= max_restart_attempts) {
    
    # If latest cycle returned no result, re-query search
    if (!iter_success) {
      q_0 <- EPM_esearch_basic_q(q_list)
      q_0 <- EPM_esearch_parse(q_0)
    }
    
    # Select next
    tryCatch({
      
      i0 <- min(which(!my_guide$completed))
      rst0 <- my_guide$retstart[i0]
      rmx0 <- my_guide$retmax[i0]
      lst0 <- list(encoding = encoding, 
                   web_env = q_0$web_env, 
                   query_key = q_0$query_key, 
                   retstart = rst0, 
                   retmax = rmx0, 
                   format = format)
      
      #debug
      #message('*', appendLF = FALSE)
      
      tmp0 <- EPM_efetch_basic_q(params = lst0)
      
      # if NULL, retry. Otherwise, write partial result.
      if (!is.null(tmp0)){
        out[[length(out) + 1]] <- tmp0
        iter_success <- TRUE
        my_guide$completed[i0] <- TRUE
        
      } else {
        iter_success <- FALSE
        attempt_n <- attempt_n + 1
        Sys.sleep(time = (3 * rltime))
      }
      
    }, error = function(e) { NULL })
    
    # Pace server requests...
    check_t1 <- Sys.time()
    check_tdiff <- as.numeric(difftime(check_t1, check_t0, units = 'sec'))
    wait_time <- rltime - check_tdiff
    if (wait_time > 0){
      Sys.sleep(wait_time)
    }
    check_t0 <- Sys.time()
  }
  
  # Finalize
  out <- tryCatch({do.call(c, out)}, error = function(e) { NULL })
  tryCatch({names(out) <- NULL}, error = function(e) { NULL })
  
  # Return
  return(out)
}




#' Submit a Query and Retrieve Results from PubMed.
#' 
#' Submit a Query to the NCBI ESearch Server, 
#' capture the response and retrieve the corresponding PubMed records from 
#' the NCBI EFetch Server. 
#' Up to the first n=10,000 records returned by the query 
#' will be retrieved (as per the NCBI policy). The operation must be completed
#' within a user-defined timeout window otherwise it will be killed.
#' 
#' 
#'
#' @param query_string String (character vector of length 1), 
#' corresponding to the query string.
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be NULL.
#' @param batch_size Integer, max number of records to be retrieved 
#' as a batch. This corresponds to the "retmax" NCBI parameter.
#' @param encoding String (character vector of length 1), encoding
#' of the resulting records (e.g., "UTF-8").
#' @param format String (character vector of length 1),
#' desired format of the Pubmed records. This must be one of the values in
#' c("xml", "medline", "uilist"). 
#' @param max_restart_attempts Integer, max number of attempts in case of a
#' failed iteration.
#' @param timeout Integer, time allowed for completing the operation
#' (in seconds).
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   qry <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   easyPubMed:::EPM_retrieve_data(qry, format = "uilist")
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Character vector including the response from the server.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_retrieve_data <- function(query_string, api_key = NULL, format = 'xml', 
                              encoding = 'UTF-8', timeout = 600,
                              batch_size = 500, max_restart_attempts = 10) {
  
  # check arguments
  stopifnot(!is.null(query_string), is.character(query_string), 
            length(query_string) == 1, !is.na(query_string), 
            nchar(query_string) > 0, 
            !is.null(batch_size), is.numeric(batch_size), 
            length(batch_size) == 1, !is.na(batch_size), 
            !is.null(max_restart_attempts), is.numeric(max_restart_attempts), 
            length(max_restart_attempts) == 1, !is.na(max_restart_attempts),
            !is.null(timeout), is.numeric(timeout), 
            length(timeout) == 1, !is.na(timeout),
            !is.null(format), is.character(format), length(format) == 1, 
            !is.null(encoding), is.character(encoding), length(encoding) == 1)
  
  if (!is.null(api_key)) {
    stopifnot(is.character(api_key), 
              length(api_key) == 1, 
              !is.na(api_key), 
              nchar(api_key) > 0)
  }
  
  # inner f(x)
  # This includes a timeout call
  # Args were evaluated before
  inner_fx <- function(query_string, api_key = NULL, format = 'xml', 
                       encoding = 'UTF-8', timeout = 180,
                       batch_size = 500, max_restart_attempts = 10) {
    
    setTimeLimit(elapsed = timeout, transient = TRUE)
    y <- NULL
    while (is.null(y)) {
      y <- tryCatch({
        EPM_esearch_efetch_seq(query_string = query_string,  
                               batch_size = batch_size, encoding = encoding, 
                               format = format, api_key = api_key,
                               max_restart_attempts = max_restart_attempts) }, 
        error = function(e) { NULL })
      if (is.null(y)) {
        Sys.sleep(1)
      }
    }
    
    return(y)
  }
  
  # And execute
  y <- tryCatch({
    inner_fx(query_string = query_string,  
             batch_size = batch_size, encoding = encoding, 
             timeout = timeout,
             format = format, api_key = api_key,
             max_restart_attempts = max_restart_attempts ) }, 
    error = function(e) { NULL })
  
  #setTimeLimit(elapsed = Inf)
  return(y)
}






#' Split A PubMed Retrieval Job into Manageable Batches.
#' 
#' Assess the number of PubMed records expected from a user-provided query
#' and split the job in multiple sub-queries if the number is bigger than 
#' "max_records_per_batch" (typically, n=10,000). 
#' Sub-queries are split according to the "Create Date" of 
#' PubMed records. This does not support splitting jobs returning more than 
#' "max_records_per_batch" (typically, n=10,000) records that have 
#' the same "Create Date" (i.e., "[CRDT]").
#' 
#' 
#'
#' @param query_string String (character vector of length 1), 
#' corresponding to the query string.
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be NULL.
#' @param max_records_per_batch Integer, maximum number of records that should
#' be expected be sub-query. This number should be in the range 1,000 
#' to 10,000 (typicall, max_records_per_batch=10,000). 
#' @param verbose logical, shall progress information be printed to console.
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   qry <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   easyPubMed:::EPM_job_split(query_string = qry, verbose = TRUE)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#'                            
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Character vector including the response from the server.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_job_split <- function(query_string, api_key = NULL, 
                          max_records_per_batch = 9999, verbose = FALSE) {
  
  
  stopifnot(is.character(query_string), length(query_string) == 1, 
            !is.na(query_string), is.numeric(max_records_per_batch),
            length(max_records_per_batch) == 1, is.logical(verbose), 
            length(verbose) == 1)
  
  max_records_per_batch[max_records_per_batch > 9999] <- 9999
  max_records_per_batch[max_records_per_batch < 1000] <- 1000
  
  buil_dtflt <- function(init_date, end_date) {
    return(paste0('"', init_date, '"[CRDT]:"', end_date, '"[CRDT]') ) 
  }
  
  ave_date <- function(init_date, end_date) {
    
    diff_days <- as.numeric(difftime(as.Date(end_date), 
                                     as.Date(init_date), units = 'days'))
    
    stopifnot(diff_days >= 1)
    
    half_diff_days <- as.integer(diff_days/2)
    if (half_diff_days < 1) {
      mid_1 <- init_date
      mid_2 <- end_date
    } else {
      mid_1 <- format(as.Date(init_date) + half_diff_days, format = '%Y/%m/%d')
      mid_2 <- format(as.Date(mid_1) + 1, format = '%Y/%m/%d')
    }
    
    y <- list(
      init = init_date, 
      mid1 = mid_1, 
      mid2 = mid_2, 
      end = end_date)
    
    return(y)
  }
  
  compute_subquery <- function(query, init_date, end_date, api_key = NULL) {
    setTimeLimit(elapsed = 120)
    qy <- paste0('(', query, ') AND (', 
                 buil_dtflt(init_date = init_date, end_date = end_date), ')')
    job0 <- NULL
    while(is.null(job0)) {
      job0 <- tryCatch({
        EPM_esearch_parse(EPM_esearch_basic_q(list(q=qy, api_key=api_key)))}, 
        error = function(e) { NULL })
      if (is.null(job0)){
        Sys.sleep(0.5)
      }
    }
    
    dif0 <- as.numeric(difftime(time1 = as.Date(end_date), 
                                time2 = as.Date(init_date), units = 'days'))       
    y <- data.frame(
      query_string = query,
      init_date = init_date,
      end_date = end_date,
      diff_days = dif0,
      exp_count = tryCatch({job0$count}, error = function(e) { 0 }), 
      stringsAsFactors = FALSE)
    
    return(y)
  }
  
  expand_rows <- function(x, max_n = 9999, api_key = NULL) {
    
    # determine waiting time
    if (!is.null(api_key)) {
      wtime <- 0.15
    } else {
      wtime <- 0.35
    }
    
    # init position and collector
    x$pos <- seq(1, length.out = nrow(x), by = 2)
    y <- list()
    
    for (i in seq(1, nrow(x))) {
      
      pos_i <- as.numeric(x$pos[i])
      diff_di <- as.numeric(x$diff_days[i])
      cnt_i <- as.numeric(x$exp_count[i])
      
      if (diff_di > 1 && cnt_i > max_n) {
        A0 <- tryCatch({
          avedt_i <- ave_date(init_date = x$init_date[i], 
                              end_date = x$end_date[i])
          
          A1 <- compute_subquery(query = x$query_string[i], 
                                 init_date = avedt_i$init, 
                                 end_date = avedt_i$mid1, api_key = api_key)
          A1$pos <- x$pos[i]
          Sys.sleep(wtime)
          
          A2 <- compute_subquery(query = x$query_string[i], 
                                 init_date = avedt_i$mid2, 
                                 end_date = avedt_i$end, api_key = api_key)
          A2$pos <- 1 + x$pos[i]
          Sys.sleep(wtime)
          
          rbind(A1, A2)}, error = function(e) {x[i, ]})
        
        y[[length(y) + 1]] <- A0
        
      } else {
        y[[length(y) + 1]] <- x[i, ]
      }
    }
    y <- do.call(rbind, y)
    y <- y[y$exp_count > 0, ]
    y <- y[order(y$pos), ]
    rownames(y) <- NULL
    
    return(y[, c(1, 2, 3, 4, 5)])
  }
  
  # Declare Init dates for the job
  init_date <- format('1850/01/01', format = '%Y/%m/%d')
  end_date <- format(1e7+Sys.time(), format = '%Y/%m/%d')
  
  # Hardcode the max number of vals
  max_recs <- max_records_per_batch
  
  TMP <- compute_subquery(query = query_string,
                          init_date = init_date, 
                          end_date = end_date, 
                          api_key = api_key)
  
  # Initial check
  tot_expected_cnt <- TMP$exp_count
  
  ## debug (keep commented)
  # message(paste0('tot_expected_cnt: ', tot_expected_cnt))
  # message(paste0('max_recs: ', max_recs))
  
  check_00 <- (tot_expected_cnt > max_recs) & 
    (TMP$diff_days > 1)
  
  # are there zero records?
  if (verbose && tot_expected_cnt < 1) {
    message('The query returned NO records!')
  }
  
  # are there too many records?
  if (sum(check_00) > 0) {
    if (verbose) {
      message(paste0('The query returned a large number of records (n=', 
                     TMP$exp_count, '). \nPlease, wait while the query is ', 
                     'further processed into a list of manegeable ', 
                     'sub-jobs. This may take several minutes.\n'))
    }
    
    # Loop
    while (is.null(TMP) || sum(check_00) > 0) {
      
      # Split
      TMP <- tryCatch({
        expand_rows(x = TMP, max_n = max_recs, api_key = api_key)}, 
        error = function(e) { TMP })
      
      # Check again
      check_00 <- (TMP$exp_count > max_recs) & 
        (TMP$diff_days > 1)
      
      # Debug
      if (verbose)
        message('.', appendLF = FALSE)
    }
  }
  
  if (verbose)
    message('', appendLF = TRUE)
  
  # Return
  r2bm <- (TMP$exp_count - 9999)
  r2bm[r2bm < 0] <- 0
  my_meta <- list(
    query_string = query_string, 
    max_records_per_batch = max_records_per_batch,
    exp_count = tot_expected_cnt,
    exp_num_of_batches = nrow(TMP),
    all_records_covered = sum(TMP$exp_count <= 9999) == nrow(TMP),
    exp_missed_records = sum(r2bm))
  
  return(list(meta = my_meta, 
              query_guide = TMP))
}





#' Validate Parameters of a PubMed Retrieval Job.
#' 
#' Check and correct (if needed) the parameters of an easyPubMed
#' retrieval job.  
#'
#' @param params list of user-provided parameters.
#' 
#' @details 
#' The following elements are expected and/or parsed from the
#' `params` list:
#' \itemize{
#' \item `encoding`. String, e.g. "UTF-8".
#' \item `format`. String, must be one of the
#' following values: `c('uilist', 'medline', 'xml')`.
#' \item `store_contents`. Logical, shall retrieved contents 
#' be stored in the object. If `FALSE`, the `write_to_file`
#' argument must be `TRUE`.
#' \item `write_to_file` Logical, shall retrieved contents
#' be written to a file (or list of files). If `FALSE`, the
#' `store_contents` argument must be `TRUE`.
#' \item `outfile_path`. String, path to the folder where 
#' files will be written. This argument is evaluated only
#' if `write_to_file` is `TRUE`.
#' \item `outfile_prefix`. String, prefix of the 
#' files that will be written locally. This argument 
#' is evaluated only
#' if `write_to_file` is `TRUE`.
#' \item `api_key`. String, NCBI API key. Can be NULL.
#' \item `max_records_per_batch`. Integer scalar (numeric 
#' vector of length 1), this is the maximum number of 
#' records retrieved per batch. It deafualts to 10,000.
#' \item `verbose`. Logical, shall details about the 
#' progress of the operation be printed to console.
#' }
#' 
#' 
#' @examples 
#' prms <- list(
#'   encoding  = 'UTF-8', 
#'   format = 'xml', 
#'   api_key = NULL,
#'   store_contents = TRUE, 
#'   write_to_file = FALSE, 
#'   verbose = TRUE)
#' easyPubMed:::EPM_validate_fetch_params(prms)
#' 
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return list including the vetted parameters.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_validate_fetch_params <- function(params) {
  
  # Get params
  stopifnot(is.list(params), length(params) > 0)
  encoding  <- params[['encoding']]
  format  <- params[['format']]
  api_key  <- params[['api_key']]
  store_contents  <- params[['store_contents']]
  write_to_file  <- params[['write_to_file']]
  outfile_path  <- params[['outfile_path']]
  outfile_prefix  <- params[['outfile_prefix']]
  verbose  <- params[['verbose']]
  
  # Parse arguments & correct
  if (is.null(encoding)) {
    encoding <- 'UTF-8'
  }
  stopifnot(!is.null(encoding), 
            is.character(encoding), 
            length(encoding) == 1, 
            !is.na(encoding), 
            nchar(encoding) > 0)
  
  stopifnot(!is.null(format), 
            is.character(format), 
            length(format) == 1, 
            !is.na(format), 
            nchar(format) > 0, 
            format %in% c('uilist', 'medline', 'xml'))
  
  stopifnot(!is.null(store_contents), 
            is.logical(store_contents), 
            length(store_contents) == 1, 
            !is.na(store_contents))
  
  stopifnot(!is.null(write_to_file), 
            is.logical(write_to_file), 
            length(write_to_file) == 1, 
            !is.na(write_to_file))
  
  stopifnot(!is.null(verbose), 
            is.logical(verbose), 
            length(verbose) == 1, 
            !is.na(verbose))
  
  if (!is.null(api_key)) {
    stopifnot(is.character(api_key), 
              length(api_key) == 1, 
              !is.na(api_key), 
              nchar(api_key) > 0)
  }
  
  if (!is.null(outfile_path)) {
    stopifnot(is.character(outfile_path), 
              length(outfile_path) == 1, 
              !is.na(outfile_path),  
              nchar(outfile_path) > 0,
              dir.exists(outfile_path))
    
  } else if (write_to_file) {
    outfile_path <- '.'
  }
  
  
  if (!is.null(outfile_prefix) && write_to_file) {
    stopifnot(is.character(outfile_prefix), 
              length(outfile_prefix) == 1, 
              !is.na(outfile_prefix), 
              nchar(outfile_prefix) > 0)
  } else if (write_to_file) {
    outfile_prefix <- paste0('easypubmed_job_', 
                             format(Sys.time(), 
                                    format = '%Y%m%d%H%M'))
  }
  
  
  # Checks for addtl requirement(s)
  stopifnot(sum(c(write_to_file,  store_contents)) > 0)
  
  # Provisional output
  y <- list(encoding = encoding, 
            format = format, 
            api_key = api_key,
            store_contents = store_contents, 
            write_to_file = write_to_file, 
            outfile_path = outfile_path, 
            outfile_prefix = outfile_prefix, 
            verbose = verbose)
  
  return(y)
}




#' Validate Parameters of a PubMed Record Parsing Job.
#' 
#' Check and correct (if needed) the parameters of an easyPubMed
#' Record Parsing job.  
#'
#' @param params list of user-provided parameters.
#' 
#' 
#' @details 
#' The following elements are expected and/or parsed from the
#' `params` list:
#' \itemize{
#' \item `max_authors`. Numeric, maximum number of authors to retrieve. If this
#' is set to -1, only the last author is extracted. If this is set to 1, 
#' only the first author is returned. If this is set to 2, the first and the 
#' last authors are extracted. If this is set to any other positive 
#' number (i), up to the leading (i-1) authors are retrieved together with the 
#' last author. If this is set to a number larger than the number of authors in
#' a record, all authors are returned. Note that at least 1 author has to be
#' retrieved, therefore a value of 0 is not accepted (coerced to -1).
#' \item `autofill_address`. Logical, shall author affiliations be 
#' propagated within each record to fill missing values.
#' \item `compact_output`. Logical, shall record data be returned in a 
#' compact format where each row is a single record and author names are
#' collapsed together. If `FALSE`, each row corresponds to a single author of
#' the publication and the record-specific data are recycled for all included
#' authors. 
#' \item `include_abstract`. Logical, shall abstract text be included in the 
#' output `data.frame`.
#' \item `max_references`. Numeric, maximum number of references to return 
#' (from each PubMed record).
#' \item `ref_id_type`. String, must be one of the
#' following values: `c('pmid', 'doi')`.
#' \item `verbose`. Logical, shall details about the 
#' progress of the operation be printed to console.
#' }
#' 
#' 
#' @examples 
#' 
#' prms <- list(
#'   max_authors  = 12, 
#'   autofill_address = TRUE, 
#'   compact_output = FALSE,
#'   include_abstract = TRUE, 
#'   max_references = 100, 
#'   ref_id_type = 'doi',
#'   verbose = TRUE)
#' easyPubMed:::EPM_validate_parse_params(prms)
#' 
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return list including the vetted parameters.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_validate_parse_params <- function(params) {
  
  # Get params
  stopifnot(is.list(params), length(params) > 0)
  
  max_authors  <- params[['max_authors']]
  autofill_address  <- params[['autofill_address']]
  compact_output  <- params[['compact_output']]
  include_abstract  <- params[['include_abstract']]
  max_references <- params[['max_references']]
  ref_id_type  <- params[['ref_id_type']]
  verbose  <- params[['verbose']]
  
  # Parse arguments & correct / set to default if missing
  
  # max_authors
  if (is.null(max_authors)) {
    max_authors <- 15
  } else {
    max_authors <- tryCatch({suppressWarnings(
      as.integer(as.numeric(max_authors)))}, 
      error = function(e) { 15 })
    
    if (!is.numeric(max_authors)) {
      max_authors <- 15
    }
    
    if (max_authors < 1) { max_authors <- (-1)}
  }
  
  # max_references
  if (is.null(max_references)) {
    max_references <- 150
  } else {
    max_references <- tryCatch({suppressWarnings(
      as.integer(as.numeric(max_references)))}, 
      error = function(e) { 150 })
    
    if (!is.numeric(max_references)) {
      max_references <- 150
    }
    
    if (max_references < 1) { max_references <- 0}
  }
  
  # autofill_address
  if (is.null(autofill_address)) {
    autofill_address <- TRUE
  } else {
    autofill_address <- tryCatch({suppressWarnings(
      as.logical(autofill_address))}, 
      error = function(e) { TRUE })
    
    if (!is.logical(autofill_address)) {
      autofill_address <- TRUE
    }
  }
  
  # compact_output
  if (is.null(compact_output)) {
    compact_output <- TRUE
  } else {
    compact_output <- tryCatch({suppressWarnings(
      as.logical(compact_output))}, 
      error = function(e) { TRUE })
    
    if (!is.logical(compact_output)) {
      compact_output <- TRUE
    }
  }
  
  # include_abstract
  if (is.null(include_abstract)) {
    include_abstract <- TRUE
  } else {
    include_abstract <- tryCatch({suppressWarnings(
      as.logical(include_abstract))}, 
      error = function(e) { TRUE })
    
    if (!is.logical(include_abstract)) {
      include_abstract <- TRUE
    }
  }
  
  # ref_id_type
  stopifnot(!is.null(ref_id_type), 
            is.character(ref_id_type), 
            length(ref_id_type) == 1, 
            !is.na(ref_id_type))
  id_choices <- c('pmid', 'doi')
  ref_id_type <- match.arg(arg = ref_id_type, choices = id_choices, 
                           several.ok = FALSE)
  
  # verbose
  if (is.null(verbose)) {
    verbose <- TRUE
  } else {
    verbose <- tryCatch({suppressWarnings(
      as.logical(verbose))}, 
      error = function(e) { TRUE })
    
    if (!is.logical(verbose)) {
      verbose <- TRUE
    }
  }
  
  # Package parsed params
  y <- list(
    max_authors = max_authors,
    autofill_address = autofill_address, 
    compact_output = compact_output, 
    include_abstract = include_abstract, 
    max_references = max_references,
    ref_id_type = ref_id_type, 
    verbose = verbose)
  
  # return
  return(y)
}






#' Detect PubMed Record Identifiers.
#' 
#' Parse a list of pubmed records in XML or Medline format, 
#' extract and return the corresponding PubMed record identifiers (PMID).
#' 
#'
#' @param x list including PubMed record data (either in `xml` or 
#' `abstract` format).
#' 
#' @param format string (character of length 1) indicating the format of
#' each element in x (either `xml` or `medline`).
#' 
#' @param as.list logical (of length 1). Shall results be returned as a list.
#' 
#' 
#' 
#' 
#' @examples 
#' x <- list(A='First record: <PMID>Rec_1A</PMID> Lorem ipsum dolor sit amet', 
#'           B='Another record: <Ti>Title</Ti><PMID>Rec_2</PMID> Lorem ipsum ')
#' easyPubMed:::EPM_detect_pmid(x, format = 'xml')
#' 
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return list of PubMed record identifiers.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_detect_pmid <- function(x, format = 'xml', as.list = TRUE) {
  
  tmp_h1 <- function(x) {
    tryCatch({
      tmp <- EPM_custom_grep(x, tag = "PMID", format = "char"); tmp[1]}, 
      error = function(e) { NA })
  }
  
  det_xml <- function(x, as.list = TRUE) {
    tmp_id <- substr(x, 1, 500)
    tmp_id <- lapply(tmp_id, tmp_h1)
    if (!as.list) {
      tmp_id <- as.character(do.call(c, tmp_id))
    }
    return(tmp_id)
  }
  
  det_mdl <- function(x, as.list = TRUE) {
    y <- lapply(x, function(xi) {
      if (sum(grepl(pattern = '^PMID', xi)) > 0) {
        ti <- grep('^PMID', xi, value = TRUE)
        ti <- sub('^PMID[[:punct:]]*[[:space:]]*([[:digit:]]+)[^[:digit:]]*$', 
                  replacement =  '\\1', x = ti)
      } else {
        ti <- NA
      }
      ti
    })   
    if (as.list) {
      return(y)
    } else {
      return(as.character(do.call(c, y)))
    }
  }
  
  if (format == 'xml') {
    y <- tryCatch({det_xml(x, as.list = as.list)}, 
                  error = function(e) { message (e); NULL})
  } else {
    y <- tryCatch({det_mdl(x, as.list = as.list)}, 
                  error = function(e) { message (e); NULL})
  }
  
  return(y)
}



#' Harmonize the Elements of a Vector by Adding Leading Zeros.
#' 
#' Coerce a vector to character and then harmonize the number of characters 
#' (nchar) of each element by adding a suitable number of leading
#' zeroes (or other user-character). 
#' 
#' @param x vector (numeric or character).
#' 
#' @param fillchar string corresponding to a single character. 
#' This character is going to be added (one or more times) 
#' in front of each element of the input vector. 
#' 
#' 
#' @examples 
#' # Example 1
#' easyPubMed:::EPM_zerofill(c(1, 100, 1000))
#' # Example 2
#' easyPubMed:::EPM_zerofill(c('Hey,', 'hello', 'there!'), '_')
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return character vector whose elements have all the same size 
#' (number of characters).
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_zerofill <- function(x, fillchar = '0') {
  
  stopifnot(is.vector(x), 
            length(x) > 0, 
            sum(is.na(x)) == 0, 
            is.vector(fillchar), length(fillchar) == 1, 
            !is.na(fillchar), nchar(fillchar) == 1)
  
  xx <- as.character(x)
  x_nc <- nchar(xx)
  max_nchar <- max(x_nc) + 1
  char_to_add <- max_nchar - x_nc
  
  xx <- lapply(seq(1, length(xx), by = 1), function(i) {
    paste(c(rep(fillchar, char_to_add[i]), xx[i]), collapse = '')
  })
  
  y <- do.call(c, xx)
  names(y) <- NULL
  return(y)
}





#' Map Job Batches to Filenames.
#' 
#' Build Filenames Matching job sub-tasks. Each filename corresponds to
#' a series of records returned by a specific job batch. The associated
#' filename indicates where the corresponding records will be written 
#' on the local disc (if requested by the user). 
#' 
#' @param job_list data.frame. This is the `job_list` data.frame included in 
#' the `misc` slot of an `easyPubMed` object.  
#' 
#' @param path folder on the local computer where files will be saved. It must 
#' be an existing directory.
#' 
#' @param prefix string used as common prefix for all files written as part of 
#' the same PubMed record download job.
#' 
#' 
#' 
#' @examples 
#' test_df <- data.frame(query_string = c('ANY', 'ANY'), 
#'                       init_date = c('2020/01/01', '2020/01/10'), 
#'                       end_date = c('2020/01/11', '2020/01/20'), 
#'                       diff_days = c(10, 10), 
#'                       exp_count = 100, 100)
#' easyPubMed:::EPM_prep_outfile(test_df, path = '.', prefix = 'my_test_job')                   
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return character vector pointing to the target files where Pubmed
#' records will be written.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_prep_outfile <- function(job_list, path , prefix) {
  
  # Initial checks
  stopifnot(!is.null(job_list), is.data.frame(job_list), 
            nrow(job_list) > 0, !is.null(path), 
            is.character(path), length(path) == 1, 
            !is.na(path), dir.exists(path), 
            !is.null(prefix), 
            is.character(prefix), length(prefix) == 1, 
            nchar(prefix)> 0) 
  
  # adjust folder
  my_path <- sub('/$', '', path)
  
  # expected files
  my_filex <- paste0(prefix, '_batch_', 
                     EPM_zerofill(seq(1, nrow(job_list), by = 1)), '.txt')
  
  # files to be written
  y <- file.path(my_path, my_filex)
  
  # return
  return(y)
}



#' Write PubMed Records to Local Files.
#' 
#' Write a list of PubMed records to a local file. If already existing, 
#' the destination file will be over-written. Original formatting of the
#' PubMed records should be declared and will be 
#' preserved in the output file. Format conversion is NOT supported.
#' 
#' @param x List including raw PubMed records.
#' 
#' @param to Path to the destination file on the local disc.
#' 
#' @param format String, format of the raw PubMed records that 
#' will be saved to the destination file (e.g., 'xml').
#' 
#' @param addon String, optional chunk of text in XML format to be written
#' to the destination file (header). This argument is only used when
#' `format` is set to 'xml'. It can be NULL.
#' 
#' @param verbose Logical, shall details about the 
#' progress of the operation be printed to console.
#' 
#' 
#' @examples 
#' test <- list('Record #1', 'Record #2')
#' outfile = './test_file.txt'
#' file.exists(outfile)
#' easyPubMed:::EPM_write_to_file(x = test, to = './test_file.txt', format = 'xml')
#' file.exists(outfile)
#' readLines(outfile)
#' unlink(outfile)
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return integer in the range c(0, 1). A result of 0 indicates that an error
#' occurred while writing the file. A result of 1 indicates that the operation 
#' was completed successfully.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_write_to_file <- function(x, to, format, addon = NULL, verbose = FALSE) {
  
  stopifnot(is.list(x), length(x) > 0, 
            is.character(to), length(to) == 1, !is.na(to), 
            is.character(format), length(format) == 1, !is.na(format))
  
  # XML chunks
  xml_init <- c(
    '<?xml version="1.0" ?>', 
    paste0('<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, ', 
           '1st January 2023//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/', 
           'out/pubmed_230101.dtd">'))
  
  if (format == 'xml' && !is.null(addon)) {
    xml_init <- c(xml_init, addon)
  }
  
  xml_init <- c(xml_init, '<PubmedArticleSet>')
  
  xml_close <- c('</PubmedArticleSet>') 
  xml_record_tags <- c('<PubmedArticle>', '</PubmedArticle>')
  
  # match format
  my_format <- match.arg(format, several.ok = FALSE, 
                         choices = c('xml', 'medline', 'uilist'))
  
  # if file exists, remove
  if (file.exists(to)) {
    if (verbose)
      message(paste0('File `', to, '` is going to be overwritten!'))
    
    tryCatch({unlink(to, force = TRUE)}, error = function(e) {
      stop('The output file already exists and cannot be overwritten!')
    })
  }
  
  if (my_format == 'xml') {
    
    # write XML head
    zz <- tryCatch({write(xml_init, file=to, append=FALSE)}, 
                   error = function(e) { NULL })
    
    # write each element of x
    for (j in seq(1, length(x), by = 1)) {
      zz <- tryCatch({write(
        paste0(xml_record_tags[1], x[[j]], xml_record_tags[2]), 
        file=to, append=TRUE)},
        error = function(e) { NULL })
    }
    
    # write XML foot
    zz <- tryCatch({write(xml_close, file=to, append=TRUE)}, 
                   error = function(e) { NULL })
    
  } else if (my_format == 'medline') {
    
    # Medline (init)
    write('', file=to, append=FALSE)
    
    for (j in seq(1, length(x), by = 1)) {
      zz <- tryCatch({write(c(x[[j]], ''), file=to, append=TRUE)},
                     error = function(e) { NULL })
    }
    
    
  } else if (my_format == 'uilist') {
    
    tmp_out <- tryCatch({do.call(c, x)}, 
                        error = function(e) { NULL })
    tmp_out <- tryCatch({tmp_out[nchar(tmp_out) > 0]}, 
                        error = function(e) { NULL })
    
    zz <- tryCatch({write(tmp_out, file=to, append=FALSE)},
                   error = function(e) { NULL })
  }
  
  # Check is the file was written and has size > 0
  chk0 <- tryCatch({as.numeric(file.info(to)$size > 0)}, 
                   error = function(e) { 0 })
  
  return(chk0)
}





#' Parse and Format a Pubmed Date Field.
#' 
#' Extract Date Information form a slice of a raw XML PubMed record. 
#' Day, month and year are returned. Months are recoded as numeric
#' if needed (e.g., `Oct` and `October` are converted to 10). If
#' month and/or day information are missing, these are imputed to 1.
#' If the year is missing, NA is returned.
#' 
#' @param x String (character vector of length 1) including an XML date
#' field from a PubMed record.
#' 
#' 
#' @examples 
#' dt0 <- '<Year>2021</Year><Month>03</Month><Day>12</Day>'
#' easyPubMed:::EPM_date_parse(dt0)
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return list including n=3 numeric elements: day, month and year.
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_date_parse <- function(x) {
  
  # Initial check(s)
  stopifnot(is.character(x), length(x) == 1, !is.na(x))
  
  # declare vars
  month_dict <- data.frame(
    id = c("jan", "january", "01", "1", 
           "feb", "february", "02", "2", 
           "mar", "march", "03", "3", 
           "apr", "april", "04", "4", 
           "may", "05", "5", 
           "jun", "june", "06", "6", 
           "jul", "july", "07", "7", 
           "aug", "august", "08", "8", 
           "sep", "september", "09", "9", 
           "oct", "october", "10", 
           "nov", "november", "11", 
           "dec", "december", "12"),
    
    val = c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4), 
            rep(5, 3), rep(6, 4), rep(7, 4), rep(8, 4), 
            rep(9, 4), rep(10, 3), rep(11, 3), rep(12, 3)))
  
  # Extract date elements
  y <- list(
    day = EPM_custom_grep(xml_data = x, tag = 'Day', format = 'char'),
    month = EPM_custom_grep(xml_data = x, tag = 'Month', format = 'char'),
    year = EPM_custom_grep(xml_data = x, tag = 'Year', format = 'char'))
  
  # Adjust Day
  if (is.null(y$day)) {
    y$day <- 1
  } else {
    y$day <- tryCatch({as.numeric(y$day)}, error = function(e) { 1 })
    y$day[is.na(y$day)] <- 1
  }
  
  # Adjust Month
  if (is.null(y$month)) {
    y$month <- 1
  } else {
    y$month <- tryCatch({month_dict$val[month_dict$id == tolower(y$month)]}, 
                        error = function(e) { 1 })
    y$month <- tryCatch({as.numeric(y$month)}, error = function(e) { 1 })
    y$month[is.na(y$month)] <- 1
  }
  
  # Adjust Year
  y$year <- tryCatch({
    ye <- y$year[1]
    if (!is.na(ye)) {
      if (nchar(ye) == 4 && grepl('(19)|(20)[[:digit:]]{2}', ye)) {
        as.numeric(ye)
      } else if (nchar(ye) == 2 && grepl('[[:digit:]]{2}', ye)) {
        as.numeric(paste0('20', ye))
      } else {
        NA
      }
    }}, error = function(e) {
      NA
    })
  
  # Coerce to missing is year could not be parsed
  if (is.na(y$year)) {
    y$day <- NA
    y$month <- NA
  }
  
  # return
  return(y)
}



#' Parse and Format Pubmed MeSH terms.
#' 
#' Extract MeSH Information form a slice of a raw XML PubMed record. 
#' Both MeSH codes and MeSH terms are returned.
#' 
#' @param x String (character vector of length 1) including an XML Mesh
#' term field/section from a PubMed record.
#' 
#' 
#' @examples 
#' msh <- paste0('<MeshHeading><DescriptorName UI=\"D000465\" >',
#'               'Algorithms</DescriptorName></MeshHeading>')
#' easyPubMed:::EPM_mesh_parse(msh)
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return list including n=2 elements (character vectors): mesh_codes and 
#' mesh_terms.
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_mesh_parse <- function(x) {
  
  # check input
  stopifnot(is.character(x), length(x) == 1, !is.na(x))
  
  # split mesh terms
  xx <- tryCatch({strsplit(x, split = '<MeshHeading') }, 
                 error = function(e) { NULL })
  xx <- tryCatch({ xx[[1]][-1] }, error = function(e) { NULL })
  
  mesh_codes <- tryCatch({
    tmp1 <- sub('^.{0,10}DescriptorName UI=\"([[:alnum:]]+)\".*$', '\\1', xx);
    tmp1 <- tmp1[nchar(tmp1) > 5 & nchar(tmp1) < 12]
    tmp1
  }, error = function(e) NA)
  
  mesh_terms <- tryCatch({
    tmp1 <- sub('^.{0,10}DescriptorName [^>]+>([^<]+)</.*$', '\\1', xx)
    tmp1 <- tmp1[nchar(tmp1) > 3 & nchar(tmp1) < 150]
    tmp1
  }, error = function(e) NA)
  
  return(list(
    mesh_codes = mesh_codes, 
    mesh_terms = mesh_terms))
}


#' Parse and Format Author Names and Affiliations.
#' 
#' Extract Author Information form a slice of a raw XML PubMed record. 
#' Last Name, First Name, Address and emails are returned. Only the first
#' address of each author is returned. A collapsed version of the
#' author list is also returned.
#' 
#' @param x String (character vector of length 1) including an XML Author List
#' section from a PubMed record.
#' 
#' @param max_authors Numeric, maximum number of authors to include. See 
#' details for additional information.
#' 
#' @param autofill Logical, shall non-missing address information be propagated 
#' to fill missing address information for other authors 
#' in the same publication. 
#' 
#' 
#' @details 
#' The value of the `max_authors` argument should be tuned to control which 
#' author information to extract from the input. If 
#' `max_authors` is set to `0`, no author information are extracted. If 
#' `max_authors` is set to `-1` (or any negative number), only information
#' corresponding to the last author are extracted. If `max_authors` is set to 
#' `+1`, only the first author information are extracted. If `max_authors`
#' is set to any other positive integer, only information for the 
#' indicated number of authors is extracted. In this case, information for both 
#' the first and the last author will be included.
#' 
#' @examples 
#' aff <- paste0('<Author><LastName>Doe</LastName><ForeName>John</ForeName>', 
#'               '<Affiliation>Univ A</Affiliation></Author>',
#'               '<Author><LastName>Doe</LastName><ForeName>Jane</ForeName>', 
#'               '<Affiliation>jane_doe@@univ_a.edu</Affiliation></Author>',
#'               '<Author><LastName>Foo</LastName><ForeName>Bar</ForeName>', 
#'               '<Affiliation>Univ B</Affiliation></Author>')
#' easyPubMed:::EPM_auth_parse(aff)
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return list including 2 elements: `authors` is a data.frame 
#' including one row for each author and n=4 columns:
#' lastname, forename, address and email; `collapsed` is a list
#' including 2 elements (each element is a string): authors and address. 
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'  
#' @importFrom utils tail  
#'
#' 
#' @keywords internal
EPM_auth_parse <- function(x, max_authors = 15, autofill = TRUE) {
  
  # nested f(x)
  inner_fx <- function(z) {
    
    emapat <- paste0("([[:alnum:]]([[:alnum:]]|\\.|\\-|\\_){2,200}", 
                     "@[[:alnum:]]([[:alnum:]]|\\.|\\-|\\_){1,200}(\\.)", 
                     "([[:alnum:]]){2,6})")
    
    lastnm <- tryCatch({
      EPM_custom_grep(xml_data = z, tag = 'LastName', format = 'char')}, 
      error = function(e) { NA })
    forenm <- tryCatch({
      EPM_custom_grep(xml_data = z, tag = 'ForeName', format = 'char')}, 
      error = function(e) { NA })
    affi <- tryCatch({
      tmp0 <- NULL
      tmp0 <- EPM_custom_grep(xml_data = z, tag = 'Affiliation', 
                              format = 'char')
      if (is.null(tmp0)) { tmp0 <- NA }
      tmp0[1]}, 
      error = function(e) { NA })
    
    addr <- affi
    emai <- NA
    
    emai_pos <- regexpr(emapat, affi)
    if (!is.na(as.numeric(emai_pos))) {
      if (as.numeric(emai_pos) > 0) {
        p01 <- as.numeric(emai_pos)
        p02 <- p01 + attributes(emai_pos)$match.length - 1
        addr <- substr(x = affi, start = 1, stop = (p01-1))
        addr <- gsub('(^[[:space:]]+)|([[:space:]]+$)', '', addr)
        emai <- substr(x = affi, start = p01, stop = p02)
      }
    }
    
    # Check if fields are empty
    lastnm <- ifelse(nchar(lastnm) < 1, NA, lastnm)
    forenm <- ifelse(nchar(forenm) < 1, NA, forenm)
    
    # Address: patch for email declaration
    ema_pat <- '[[:space:]]*Electronic[[:space:]]*address[[:punct:]]*$'
    addr <- tryCatch({sub(ema_pat, '', addr, ignore.case = TRUE)})
    addr <- ifelse(nchar(addr) < 1, NA, addr)
    
    emai <- ifelse(nchar(emai) < 1, NA, emai)
    
    # Package and return
    Y <- tryCatch({data.frame(
      lastname = lastnm, forename = forenm, 
      address = addr, email = emai, stringsAsFactors = FALSE)}, 
      error = function(e) { NULL })
    
    return(Y)
  }
  
  # Dummy outputs
  y1 <- data.frame(
    lastname = NA, forename = NA, 
    address = NA, email = NA, stringsAsFactors = FALSE)
  y2 <- list(authors = NA, address = NA)
  Y0 <- list(authors = y1, collapsed = y2)
  
  # Args check
  stopifnot(!is.null(x), is.character(x), 
            length(x) == 1, !is.na(x), 
            is.numeric(max_authors), length(max_authors) == 1, 
            !is.na(max_authors))
  
  max_authors <- as.integer(max_authors)
  
  if (max_authors == 0) {
    
    # return NULL    
    return(Y0)
    
  } else {
    
    add_dots <- 0
    xx <- tryCatch({strsplit(x = x, split = '</Author>', fixed = TRUE)[[1]]}, 
                   error = function(e) { NULL } )
    
    # if max_aut
    if (max_authors < 0) {
      
      # we also need to save the first affiliation
      # in case something autofill is TRUE and 
      # max_authors is -1.
      tmp_add1 <- tryCatch({inner_fx(xx[1])}, error = function(e) { NA })
      tmp_add1 <- tryCatch({tmp_add1$address[1]}, error = function(e) { NA })
      
      # only get last author
      xx <- xx[length(xx)]
      
    } else if (max_authors == 1) {
      xx <- xx[1]
      
    } else if (length(xx) > max_authors) {
      
      idx <- c(seq(1, (max_authors-1), by = 1), length(xx))
      xx <- xx[idx]
      add_dots <- 1
      
    } 
    
    AU0 <- tryCatch({
      
      # Create a guide data.frame that informs about
      # how to propagate affiliation info for selected authors.
      tmp0 <- NULL
      tmp0 <- do.call(rbind, lapply(xx, inner_fx))
      
      # Handle exceptions. If affiliation is missing for the first author
      # propagate from the last author. If both are NA, NA is used. 
      if (nrow(tmp0) > 1 && autofill &&
          (is.na(tmp0$address[1]) || nchar(tmp0$address[1]) < 1)) {
        
        tmp0$address[1] <- tmp0$address[nrow(tmp0)]
      }
      
      if (autofill) {
        
        if (max_authors < 0) {
          if (is.na(tmp0$address[nrow(tmp0)])) {
            tmp0$address <- tmp_add1
          }
          
        } else {
          aut_2_fill <- data.frame(pos = seq_len(nrow(tmp0)), 
                                   fld = as.numeric(!is.na(tmp0$address)), 
                                   final = 0)
          
          # Update the guide so that position 1 is never fld==0.
          aut_2_fill$fld[1] <- 1

          for (j in seq_len(nrow(aut_2_fill))) {
            
            if (aut_2_fill$fld[j] == 1) {
              aut_2_fill$final[j] <- aut_2_fill$pos[j]
              
            } else {
              
              log_jj <- as.logical(aut_2_fill$fld==1 & aut_2_fill$pos<j)
              log_jj[1] <- TRUE
              jj <- tryCatch({max(which(log_jj))}, error = function(e) { NA })
              
              aut_2_fill$final[j] <- aut_2_fill$pos[jj]
            }
          }
          
          tmp0$address <- tmp0$address[aut_2_fill$final]          
        }
      }
      tmp0  
      
    }, error = function(e) { NULL })
    
    # collapsed version
    au_last <- tryCatch({AU0$lastname}, error = function(e) { '' })
    au_last [is.na(au_last)] <- ''
    
    au_first <- tryCatch({AU0$forename}, error = function(e) { '' })
    au_first [is.na(au_first)] <- ''
    
    if (add_dots == 1) {
      au_cmp <- tryCatch({paste(paste(
        c(au_first[-length(au_first)], '', au_first[length(au_first)]), 
        c(au_last[-length(au_last)], '...', au_last[length(au_last)])), 
        collapse = ', ')})
      
    } else {
      au_cmp <- tryCatch({paste(paste(au_first, au_last), collapse = ', ')})
    }
    
    au_add <- tryCatch({AU0$address}, error = function(e) { '' })
    au_add <- au_add[!is.na(au_add)]
    if (length( au_add ) > 0) {
      au_add <- utils::tail(au_add, 1)
    } else {
      au_add <- NA
    }
    
    # Return
    OU1 <- list(
      authors = AU0, 
      collapsed = list(authors = au_cmp, address = au_add)
    )
    
    return(OU1)
  }   
}


#' Parse and Format References.
#' 
#' Extract Reference Information form a raw XML string, typically 
#' extracted from a PubMed record. Users can select the type of identifier
#' to extract and return, as well as the maximum number of references to 
#' be returned. 
#'  
#' 
#' @param x String (character vector of length 1) including a List
#' of references obtained from a PubMed record.
#' 
#' @param max_references Numeric (of length 1). Maximum number of references 
#' to extract/include. This should be an integer `>=0`.
#' 
#' @param id_type String (character vector of length 1). Type of identifier 
#' to be used for references. One of the following values is expected: 
#' `c('pmid', 'doi', 'pmc')`.
#' 
#' 
#' 
#' @examples 
#' ref <- paste0('<xml><Reference><Citation>',
#'               '<ArticleId IdType=\"pubmed\">25822800</ArticleId>',
#'               '<ArticleId IdType=\"pmc\">PMC4739640</ArticleId>',
#'               '</Citation></Reference></xml>')
#' easyPubMed:::EPM_reference_parse(ref)
#' easyPubMed:::EPM_reference_parse(ref, id_type = 'pmc')
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return data.frame including one row for each author and n=4 columns:
#' lastname, forename, address and email. 
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_reference_parse <- function(x, max_references = 100, id_type = 'pmid') {
  
  stopifnot(is.character(x), length(x) == 1, !is.na(x), nchar(x) > 10, 
            is.numeric(max_references), length(max_references) == 1, 
            !is.na(max_references), is.character(id_type), 
            length(id_type) == 1, !is.na(id_type))
  
  id_type <- match.arg(arg = id_type, choices = c('pmid', 'doi', 'pmc'), 
                       several.ok = FALSE)
  
  # If we do not wish to return references, output NA
  if (max_references < 1) {
    
    return(NA)
  }
    
  if (id_type == 'pmid') {
    xx <- tryCatch({
      unique(EPM_custom_grep(x, tag = 'ArticleId', 
                             xclass = 'pubmed', format = 'char'))}, 
      error = function(e) { NULL } )
    
  } else if (id_type == 'doi') {
    xx <- tryCatch({
      unique(EPM_custom_grep(x, tag = 'ArticleId', 
                             xclass = 'doi', format = 'char'))}, 
      error = function(e) { NULL } )
    
  } else if (id_type == 'pmc') {
    xx <- tryCatch({
      unique(EPM_custom_grep(x, tag = 'ArticleId', 
                             xclass = 'pmc', format = 'char'))}, 
      error = function(e) { NULL } )
  }
  
  if (is.null(xx)) { 
    xx <- NA 
  } else if (length(xx) > 0) {
    keep <- seq_len(min(c(length(xx), max_references)))
    xx <- xx[keep]
  } 
  
  # return
  return(xx)
}  






#' Generate a Unique Query Key.
#' 
#' Generate a pseudo-random key that uniquely identifies easyPubMed 
#' objects. The key is a 46-char string
#' that includes the current date + time and a list of 
#' randomly selected characters, numbers and special characters.
#' The unique key is typically saved in the `meta` slot of an
#' easyPubMed object, and is also written to local files when
#' records are donwloaded and saved in XML format.
#' This function takes NO arguments.
#'  
#' 
#' @examples 
#' easyPubMed:::EPM_init_unique_key()
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return string, a 46-char unique key.
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_init_unique_key <- function() {
  
  pr1 <- 'EPMJ_'
  pl1 <- c(rep(letters,2), 
           rep(seq(0, 9, by = 1),4), 
           c('!', '@', '#', '$', '%', '&'))
  pl1 <- sample(pl1, size = length(pl1), replace = FALSE)
  
  da1 <- format(Sys.time(), '%Y%m%d%H%M%S')
  ch1 <- sample(pl1, size = 8, replace = TRUE)
  ch2 <- sample(pl1, size = 8, replace = TRUE)
  ch3 <- sample(pl1, size = 8, replace = TRUE)
  ch4 <- sample(letters, size = 2, replace = TRUE)
  
  my_uid <- paste0(pr1, da1, '_',
                   paste(ch1, collapse = ''), 
                   paste(ch2, collapse = ''),
                   paste(ch3, collapse = ''), 
                   paste(ch4, collapse = ''))
  
  return(my_uid)
}







#' Encode Metadata to an XML String.
#' 
#' Encode a list of meta information from an easyPubMed object 
#' into an XML string. These meta-information are used to keep track 
#' of easyPubMed query jobs and/or to re-build objects
#' starting from XML files saved on a local disk.
#' 
#' 
#' @param meta List including metadata associated with an easyPubMed query job. 
#' It corresponds to the contents of the `meta` slot of an
#' easyPubMed object.
#' 
#' @param job_list Data.frame that defines the list of sub-queries of an 
#' easyPubMed query job. It corresponds to the
#' 'job_list' data.frame included in the `misc` slot of an
#' easyPubMed object.
#'  
#' 
#' @param i Integer, index of the batch (query sub-job) being written to file.
#' 
#' @param encoding String, this is the Encoding of the contents/text being 
#' retrieved from the Entrez server (typically, 'UTF-8').  
#' 
#' 
#' 
#' @examples 
#' tmp_meta <- list(max_records_per_batch = 1000, 
#'                  exp_count = 10, 
#'                  exp_num_of_batches = 1, 
#'                  all_records_covered = TRUE, 
#'                  exp_missed_records = 0, 
#'                  query_date = "2023-10-16 23:13:29", 
#'                  UID = 'EPMJ_20231017141741_c4das',
#'                  EPM_version = "3.01")
#' tmp_jobs <- data.frame(query_string = 'my test query', 
#'                        init_date = '1990/01/01', 
#'                        end_date = '2023/01/01', 
#'                        diff_days = 12053,
#'                        exp_count = 10, 
#'                        stringsAsFactors = FALSE)
#' easyPubMed:::EPM_encode_meta_to_xml(meta = tmp_meta, job_list = tmp_jobs, 
#'                                     i = 1, encoding = 'UTF-8' )
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return String, chunck of XML-decorated text including meta information.
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_encode_meta_to_xml <- function(meta, job_list, i, encoding) {
  
  stopifnot(is.list(meta), length(meta) >= 8,  
            is.data.frame(job_list),
            nrow(job_list) > 0,
            is.numeric(i), length(i) == 1, 
            i <= nrow(job_list), i>= 1)
  
  y <- tryCatch({
    paste0(
      '<EPMxJobData>',
      '<EPMxJobUniqueId>', meta$UID,'</EPMxJobUniqueId>',
      '<EPMxJobBatchNum>', nrow(job_list),'</EPMxJobBatchNum>',
      '<EPMxJobBatch>', i,'</EPMxJobBatch>',
      
      '<EPMxQuery>', job_list$query_string[i],'</EPMxQuery>',
      '<EPMxQBatchInitDate>', job_list$init_date[i],'</EPMxQBatchInitDate>',
      '<EPMxQBatchEndDate>', job_list$end_date[i],'</EPMxQBatchEndDate>',
      '<EPMxQBatchDiffDays>', job_list$diff_days[i],'</EPMxQBatchDiffDays>',
      '<EPMxQBatchExpCount>', job_list$exp_count[i],'</EPMxQBatchExpCount>',
      
      '<EPMxMaxRecordsPerBatch>', meta$max_records_per_batch ,
      '</EPMxMaxRecordsPerBatch>',
      '<EPMxExpCount>', meta$exp_count,'</EPMxExpCount>',
      '<EPMxExpNumOfBatches>', meta$exp_num_of_batches,
      '</EPMxExpNumOfBatches>',
      '<EPMxAllRecordsCovered>', meta$all_records_covered,
      '</EPMxAllRecordsCovered>',
      '<EPMxExpMissedRecords>', meta$exp_missed_records,
      '</EPMxExpMissedRecords>',
      '<EPMxQueryDate>', meta$query_date,'</EPMxQueryDate>',
      '<EPMxRawFormat>', 'xml','</EPMxRawFormat>',
      '<EPMxRawEncoding>', encoding,'</EPMxRawEncoding>',
      '<EPMxRawDate>', format(Sys.time(), tz = 'GMT'),'</EPMxRawDate>',
      '<EPMxLibVersion>', meta$EPM_version,'</EPMxLibVersion>',
      
      '</EPMxJobData>')}, 
    error = function(e) { NULL })
  
  return(y)
}



#' Decode an XML String into the Corresponding Metadata.
#' 
#' Decode an XML String including a list of meta information 
#' associated to an easyPubMed object whose contents were written to 
#' a text file on a local disk. These meta-information are used to keep track 
#' of easyPubMed query jobs and/or to re-build objects
#' starting from XML files saved on a local disk.
#' 
#' 
#' @param x String corresponding to the XML-decorated text including
#' metadata from an easyPubMed object/query job.
#' 
#' 
#' @examples 
#' xml <- paste0('<EPMxJobData><EPMxJobUniqueId>EPMJ_20231017151112_mi7xvol743', 
#'               'rvz5ry5z3n8qm0ww</EPMxJobUniqueId><EPMxJobBatchNum>4</EPMxJo', 
#'               'bBatchNum><EPMxJobBatch>1</EPMxJobBatch><EPMxQuery>Test_Quer', 
#'               'y</EPMxQuery><EPMxQBatchInitDate>1937/01/22</EPMxQBatchInitD', 
#'               'ate><EPMxQBatchEndDate>1980/08/01</EPMxQBatchEndDate><EPMxQB', 
#'               'atchDiffDays>15897</EPMxQBatchDiffDays><EPMxQBatchExpCount>2', 
#'               '13</EPMxQBatchExpCount><EPMxMaxRecordsPerBatch>1000</EPMxMax', 
#'               'RecordsPerBatch><EPMxExpCount>2083</EPMxExpCount><EPMxExpNum', 
#'               'OfBatches>4</EPMxExpNumOfBatches><EPMxAllRecordsCovered>TRUE', 
#'               '</EPMxAllRecordsCovered><EPMxExpMissedRecords>0</EPMxExpMiss', 
#'               'edRecords><EPMxQueryDate>2023-10-17 15:11:12</EPMxQueryDate>', 
#'               '<EPMxRawFormat>xml</EPMxRawFormat><EPMxRawEncoding>UTF-8</EP', 
#'               'MxRawEncoding><EPMxRawDate>2023-10-17 15:14:12</EPMxRawDate>', 
#'               '<EPMxLibVersion>3.01</EPMxLibVersion></EPMxJobData>')
#' easyPubMed:::EPM_decode_xml_meta(xml)
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return String, chunck of XML-decorated text including meta information.
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_decode_xml_meta <- function(x) {
  
  # Check
  stopifnot(is.character(x), length(x) == 1, !is.na(x), nchar(x) > 0)
  
  # Declare tags of interest
  char_tags <- c('EPMxJobUniqueId', 'EPMxQuery', 'EPMxQBatchInitDate', 
                 'EPMxQBatchEndDate', 'EPMxQueryDate', 'EPMxRawFormat', 
                 'EPMxRawEncoding', 'EPMxRawDate', 'EPMxLibVersion')
  num_tags <- c('EPMxJobBatchNum', 'EPMxJobBatch', 'EPMxQBatchDiffDays', 
                'EPMxQBatchExpCount', 'EPMxMaxRecordsPerBatch', 'EPMxExpCount', 
                'EPMxExpNumOfBatches', 'EPMxExpMissedRecords')
  log_tags <- c('EPMxAllRecordsCovered')
  
  # initialize collector
  y <- list()
  
  # Loop
  for (ci in char_tags) {
    y[[sub('^EPMx', '', ci)]] <- as.character(
      EPM_custom_grep(xml_data = x, tag = ci, format = 'char'))
  }
  for (ci in num_tags) {
    y[[sub('^EPMx', '', ci)]] <- as.numeric(
      EPM_custom_grep(xml_data = x, tag = ci, format = 'char'))
  }
  for (ci in log_tags) {
    y[[sub('^EPMx', '', ci)]] <- as.logical(
      EPM_custom_grep(xml_data = x, tag = ci, format = 'char'))
  }
  
  return(y)
}



#' Import PubMed Records Saved Locally in XML Format.
#' 
#' Read the contents of an XML file and import Metadata and
#' PubMed records for use by easyPubMed. 
#' The XML file must be generated by easyPubMed (ver >= 3) via 
#' the `epm_fetch()` function or via the `fetchEPMData()` method. 
#' XML files downloaded from the Web or using other software are 
#' currently unsupported. This function can only process one file.
#' 
#' 
#' 
#' @param x Path to an XML file on the local machine. 
#' 
#' 
#' @examples 
#' \dontrun{
#'   x <- epm_query(query_string = 'easyPubMed', verbose = TRUE)
#'   x <- epm_fetch(x = x, write_to_file = TRUE, store_contents = FALSE, 
#'                  outfile_prefix = 'qpm_qry_', verbose = TRUE)
#'   y <- EPM_read_xml(x = 'qpm_qry__batch_01.txt')
#'   try(unlink('qpm_qry__batch_01.txt'), silent = TRUE)
#'   y
#' }
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return List including four elements: `guide` (data.frame), `meta` (list),      
#' `job_info` (data.frame) and `contents` (named list). 
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_read_xml <- function(x) {
  
  # nested f(x)
  inner_fx <- function(x, y) {
    
    pats <- list('<?xml', '<!doctype pubmedarticleset', 
                 '<epmxjobdata>', '<pubmedarticleset', 
                 '<pubmedarticle>')
    
    # init inner_fx output
    out <- list(check = y, type = 0, id = NA, data = NA)
    
    # isolate initial string chunck
    sx <- tolower(substr(x, 1, 30))
    
    # check 1
    if (grepl(pats[[1]], sx, fixed = TRUE)) {
      
      out$check$xml <- TRUE
      
      # check 2  
    } else if (grepl(pats[[2]], sx, fixed = TRUE)) {
      
      out$check$doc <- TRUE
      
      # check 3
    } else if (grepl(pats[[3]], sx, fixed = TRUE)) {
      
      out$check$epmx <- TRUE
      out$data <- EPM_decode_xml_meta(x)
      out$type <- 1
      
      # check 4
    } else if (grepl(pats[[4]], sx, fixed = TRUE)) {
      
      out$check$pmas <- TRUE
      
      # check 5
    } else if (sum(as.numeric(do.call(c, out$check))) == 4 & 
               grepl(pats[[5]], sx, fixed = TRUE)) {
      
      out$data <- x
      out$id <- EPM_detect_pmid(x, as.list = FALSE)
      out$type <- 2
      
    }
    
    return(out)
  }
  
  
  # Check if file exists
  # Must be a single file
  stopifnot(is.character(x), length(x) == 1, 
            !is.na(x), file.exists(x))
  
  # Create connection
  fcon <- file(x, "r")
  
  # Init Checklist
  chk <- list(xml = FALSE, doc = FALSE, epmx = FALSE, pmas = FALSE)
  
  # init collector
  my_rids <- list()
  my_cnts <- list()
  my_meta <- list()
  
  # Inint counter
  li_num <- 1
  
  # Start reading 
  # make sure to stop if checks fail within the first 25 lines
  tryCatch({
    while ( TRUE ) {
      
      li <- readLines(fcon, n = 1)
      if ( length(li) == 0 ) {
        break
      }
      if ( li_num > 25 && sum(as.logical(tmp_li$check)) < 4) {
        break
      }
      
      tmp_li <- inner_fx(li, chk)
      chk <- tmp_li$check 
      
      if (tmp_li$type == 1) {
        
        my_meta <- tmp_li$data
        
      } else if (tmp_li$type == 2) {
        
        if (!is.na(tmp_li$id) & nchar(tmp_li$id) >= 1) {
          my_rids[[length(my_rids) + 1]] <- tmp_li$id
          my_cnts[[length(my_cnts) + 1]] <- tmp_li$data
          
        }
      }
      li_num <- li_num + 1
    }
  }, error = function(e) { NULL }, finally = {close(fcon)})
  
  
  # Re-arrange output
  out_meta <- list(
    max_records_per_batch=my_meta$MaxRecordsPerBatch, 
    exp_count=my_meta$ExpCount,
    exp_num_of_batches=my_meta$ExpNumOfBatches,
    all_records_covered=my_meta$AllRecordsCovered,
    exp_missed_records=my_meta$ExpMissedRecords,
    query_date=my_meta$QueryDate,
    raw_format=my_meta$RawFormat,
    raw_encoding=my_meta$RawEncoding,
    raw_data_embedded=TRUE,
    raw_record_num=length(my_cnts),
    raw_date=my_meta$RawDate,
    processed_compact_output=NA,
    processed_record_num=NA,
    processed_ref_type=NA,
    UID=my_meta$JobUniqueId,
    EPM_version=my_meta$LibVersion)
  
  out_job_info <- data.frame(
    query_string=my_meta$Query,  
    init_date=my_meta$QBatchInitDate,   
    end_date=my_meta$QBatchEndDate, 
    diff_days=my_meta$QBatchDiffDays, 
    exp_count=my_meta$QBatchExpCount, 
    stringsAsFactors = FALSE)
  
  out_guide <- data.frame(
    JobUniqueId=my_meta$JobUniqueId,
    JobQuery=my_meta$Query,
    JobBatch=my_meta$JobBatch,
    JobBatchNum=my_meta$JobBatchNum, 
    stringsAsFactors = FALSE)
  
  names(my_cnts) <- as.character(do.call(c, my_rids))
  
  # Package output
  return(list(guide = out_guide, 
              meta = out_meta, 
              job_info = out_job_info, 
              contents = my_cnts))
}



#' Check Metadata from Imported XML Files. 
#' 
#' Analyze the Metadata from different XML files that were imported using
#' easyPubMed and identify which records / files can be merged together and
#' which ones to exclude. Only files with the same unique ID can be merged 
#' together a this step. The goal is to re-build a consistent easyPubMed
#' object. 
#' 
#' 
#' 
#' @param x Data.frame including information from the imported XML files. The
#' following columnnames are expected: `index`, `file`, `JobUniqueId`, 
#' `JobQuery`, `JobBatch`.
#' 
#' 
#' @examples 
#' gx <- data.frame(
#'   index = c(1, 2, 3, 4, 5),
#'   JobUniqueId = rep('xyz0x', 5),
#'   JobQuery = rep('test_query', 1),
#'   JobBatch = c(1, 2, 3, 4, 3),
#'   JobBatchNum = rep(4, 5),
#'   stringsAsFactors = FALSE)
#' easyPubMed:::EPM_check_guide(gx)
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return Data.frame identical to `x` with an additional *numeric) column 
#' (`pass` column). 
#' 
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @keywords internal
EPM_check_guide <- function(x) {
  
  # check input x
  stopifnot(is.data.frame(x), nrow(x) > 0, ncol(x) > 3, 
            'JobUniqueId' %in% colnames(x),
            'JobBatchNum' %in% colnames(x),
            'JobBatch' %in% colnames(x) )
  
  # add PASS column, default to 1
  x$pass <- 1
  
  # what is the majority ID
  epmj_id <- names(sort(table(x$JobUniqueId), decreasing = TRUE))[1]
  x$pass[x$JobUniqueId != epmj_id] <- 0
  
  # expected Batch Number
  epmj_bn <- as.numeric(names(sort(table(
    x$JobBatchNum[x$pass == 1]), decreasing = TRUE)))[1]
  x$pass[x$JobBatchNum != epmj_bn] <- 0
  
  # Check batch Is
  unq_bi <- unique(x$JobBatch[x$pass == 1])
  unq_bi <- unq_bi[unq_bi <= epmj_bn]
  
  lil <- list()
  for (j in seq_len(nrow(x))) {
    if( x$pass[j] == 1 ) {
      if (is.null(lil[[as.character(x$JobBatch[j])]])) {
        lil[[as.character(x$JobBatch[j])]] <- 1  
      } else {
        x$pass[j] <- 0
      }
    }
  }
  
  # Return 
  return(x)
}


## ----------------
## --- Exported ---
## ----------------



#' Search for PubMed Records.
#' 
#' Query PubMed (Entrez) via the PubMed API eSearch utility. 
#' Calling this function results in submitting a query to the NCBI EUtils 
#' server and then capturing and parsing the response.
#' The number of records expected to be returned by the query is 
#' determined. If this number is bigger than n=10,000, the record retrieval job 
#' is automatically split in a list of smaller manageable sub-queries. 
#' This function returns an "easyPubMed" object, which includes all
#' information required to retrieve PubMed records using the \code{epm_fetch()}
#' function. 
#' 
#' 
#' @details 
#' This function will use "query_string" for querying PubMed. 
#' The Query Term can include one or multiple words, as well as the standard 
#' PubMed operators (AND, OR, NOT) and tags (i.e., [AU], [PDAT], 
#' [Affiliation], and so on). 
#' 
#'
#' @param query_string String (character vector of length 1), 
#' corresponding to the query string.
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be `NULL`.
#' @param verbose logical, shall progress information be printed to console.
#' Defaults to `TRUE`.
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   qry <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
#'   epm_query(query_string = qry, verbose = FALSE)
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
#' @importFrom methods new
#'
#' 
#' @export
epm_query <- function (query_string, api_key = NULL, verbose = TRUE) {
  
  # Initial checks
  stopifnot(!is.null(query_string), 
            is.character(query_string), 
            !is.na(query_string),
            length(query_string) == 1, 
            nchar(query_string) > 0)
  
  if (!is.null(api_key)) {
    stopifnot(is.character(api_key), 
              length(api_key)  == 1, 
              !is.na(api_key), 
              nchar(api_key) > 0)  
  }
  # max records
  # set this to 9999 for production
  max_recs <- 9999
  
  # set this to 999 for debugging
  #max_recs <- 999
  
  # Assess the query job at hand
  tmp <- EPM_job_split(query_string = query_string, api_key = api_key, 
                       max_records_per_batch = max_recs, 
                       verbose = verbose)
  
  # Build an EPM object, then return
  y <- tryCatch({
    methods::new("easyPubMed", 
                 query_string = query_string, 
                 job_info = tmp) }, 
    error = function(e) { NULL })
  
  return(y)
}











#' Query PubMed by PMIDs.
#' 
#' Query PubMed using a list of PubMed record identifiers (PMIDs) as input. 
#' The list of identifiers is automatically split into a series of 
#' manageable-sized chunks (max n=50 PMIDs per chunk).
#' 
#' 
#' 
#' @param pmids Vector (character or numeric), list of 
#' Pubmed record identifiers (PMIDs). Values will be coerced to
#' character.
#' 
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be `NULL`.
#' 
#' @param verbose Logical, shall details about the 
#' progress of the operation be printed to console.
#' 
#' 
#' 
#' @examples
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   my_pmids <- c(34097668, 34097669, 34097670)
#'   epm_query_by_pmid(my_pmids)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return an easyPubMed object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @export
epm_query_by_pmid <- function(pmids, api_key = NULL, verbose = TRUE) {
  
  # Initial arg checks
  stopifnot(is.vector(pmids), length(pmids) > 0, 
            is.logical(verbose), length(verbose) == 1, 
            !is.na(verbose))
  
  if (!is.null(api_key)) {
    stopifnot(is.character(api_key), length(api_key) == 1, 
              !is.na(api_key), nchar(api_key) > 0)  
    waitime <- 0.15
  } else {
    waitime <- 0.35
  }
  
  # Handle pmids
  pmids <- as.character(pmids)
  pmids <- pmids[!is.na(pmids)]
  pmids <- pmids[nchar(pmids) > 0]
  stopifnot(is.vector(pmids), length(pmids) > 0)
  
  # split pmids by batches (of n=50)  
  ba_guide <- seq(1, length(pmids), by = 50)
  ba_guide <- data.frame(start=ba_guide, stop=(ba_guide+50-1))
  ba_guide$stop[ba_guide$stop > length(pmids)] <- length(pmids)
  
  # build subqueries
  init_date <- format('1850/01/01', format = '%Y/%m/%d')
  end_date <- format(1e7+Sys.time(), format = '%Y/%m/%d')
  dif0 <- as.numeric(difftime(time1 = as.Date(end_date), 
                              time2 = as.Date(init_date), units = 'days'))       
  
  y <- lapply(seq_along(nrow(ba_guide)), function(i) {
    tmp_ids <- pmids[seq(ba_guide$start[i], ba_guide$stop[i], by = 1)]
    tmp_qry <- paste0(tmp_ids, "[PMID]")
    tmp_qry <- paste(tmp_qry, collapse = " OR ")
    tmp_n <- NULL
    cnt_i <- 0
    while(is.null(tmp_n) && cnt_i < 10) {
      Sys.sleep(waitime)
      tmp_n <- tryCatch({
        EPM_esearch_basic_q(params = list(q=tmp_qry, api_key=api_key))}, 
        error = function(e) { NULL })
      if(is.null(tmp_n)) {
        Sys.sleep(1)
      }
      cnt_i <- cnt_i + 1
    }
    tmp_n <- tryCatch({EPM_esearch_parse(tmp_n)$count}, 
                      error = function(e) {NA})
    
    # return
    data.frame(
      query_string = tmp_qry, 
      init_date = init_date, 
      end_date = end_date, 
      diff_days = dif0, 
      exp_count = tmp_n, 
      stringsAsFactors = FALSE)
  })
  
  # Prepare elements
  qry_text <- 'Custom query (epm_query_by_pmid)'
  job_list <- do.call(rbind, y)
  miss_recs <- length(pmids) - sum(job_list$exp_count, na.rm = TRUE)
  miss_recs[miss_recs < 0] <- 0
  
  job_info <- list(
    meta = list(
      query_string = qry_text, 
      max_records_per_batch = 50, 
      exp_count = sum(job_list$exp_count, na.rm = TRUE), 
      exp_num_of_batches = nrow(job_list), 
      all_records_covered = as.logical(miss_recs > 0), 
      exp_missed_records = miss_recs),
    
    query_guide = job_list)
  
  if (miss_recs > 0 && verbose) {
    
    message(paste0('There is a mismatch between the number of PMIDs provided ', 
                   'by the user (n=', length(pmids),
                   ') and the number of expected ', 
                   'results (n=', 
                   sum(job_list$exp_count, na.rm = TRUE), ')'))
  }
  
  # Return object
  out <- tryCatch({
    methods::new("easyPubMed", 
                 query_string = qry_text, 
                 job_info = job_info) }, 
    error = function(e) { NULL })
  
  return(out)
}


#' Query PubMed by Full-length Title.
#' 
#' Execute a PubMed query using a full-length 
#' publication title as query string. Tokenization and 
#' stopword removal is automatically performed. The goal is to mimic 
#' a Pubmed citation matching search. Because of this approach, 
#' it is possible that a query by full-length title may return 
#' more than one record. 
#' 
#' 
#' @param fulltitle String (character vector of length 1) that corresponds 
#' to the full-length publication title used for querying PubMed 
#' (titles should be used as is, without adding trailing filter tags).
#' 
#' @param field String (character vector of length 1). This indicates the 
#' PubMed record field where the full-length string (fulltitle) should be 
#' searched in. By default, this points to the 'Title' field. 
#' However, the field can be changed (always use fields supported by PubMed) as 
#' required by the user (for example, to attempt an exact-match query 
#' using a specific sentence included in the abstract of a record).
#' 
#' @param api_key String (character vector of length 1), 
#' corresponding to the NCBI API key. Can be `NULL`.
#' 
#' @param verbose Logical, shall details about the 
#' progress of the operation be printed to console.
#' 
#' 
#' @examples
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   q <- 'Analysis of Mutational Signatures Using the mutSignatures R Library.'
#'   epm_query_by_fulltitle(q)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return an easyPubMed object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @export
epm_query_by_fulltitle <- function(fulltitle, field = "[Title]", 
                                   api_key = NULL, verbose = TRUE) {
  
  # Initial arg checks
  stopifnot(is.character(fulltitle), length(fulltitle) == 1, 
            is.logical(verbose), length(verbose) == 1, 
            !is.na(verbose))
  
  if (!is.null(api_key)) {
    stopifnot(is.character(api_key), length(api_key) == 1, 
              !is.na(api_key), nchar(api_key) > 0)  
    waitime <- 0.15
  } else {
    waitime <- 0.35
  }
  
  # Get stopWords
  stopwords <- easyPubMed::epm_stopwords
  
  # Try with fulltitle as is 
  qry1 <- paste('"', fulltitle, '"', field,  sep = "")
  out1 <- epm_query(query_string = qry1, api_key = api_key, verbose = FALSE)
  out1_meta <- getEPMMeta(out1)
  
  if (as.numeric(out1_meta$exp_count) > 0) {
    
    if (verbose && out1_meta$exp_count > 1) {
      message(paste0('Note: multiple records matched this fulltitle. ', 
                     'We recommend searching by PMID when specific ', 
                     'PubMed records should be retrieved.'))
    }
    return(out1)
  }
  
  # Try with stopword-removed fulltitle
  keys <- strsplit(fulltitle, split = "[[:space:]]")[[1]]
  keys <- tolower(keys)
  keys <- keys[!keys %in% stopwords]
  Sys.sleep(0.35)
  qry2 <- paste(unique(keys), field, sep = "", collapse = " AND ")
  out2 <- epm_query(query_string = qry2, api_key = api_key, verbose = FALSE)
  out2_meta <- getEPMMeta(out2)
  
  if (out2_meta$exp_count > 0) {
    
    if (verbose && out2_meta$exp_count > 1) {
      message(paste0('Note: multiple records matched this fulltitle. ', 
                     'We recommend searching by PMID when specific ', 
                     'PubMed records should be retrieved.'))
    }
    out2 <- setEPMQuery(x = out2, qry1)
    return(out2)
  }
  
  # Else, return NULL
  if (verbose && out1_meta$exp_count < 1 && out2_meta$exp_count < 1) {
    message(paste0('Unfortunately, no PubMed record was found. ', 
                   'We recommend searching by PMID when specific ', 
                   'PubMed records should be retrieved.'))
  }
  return(NULL)
}





#' Extract Information from a Raw PubMed Record.
#' 
#' Read a raw PubMed record, identify XML tags, extract information 
#' and cast it into a structured `data.frame`. The expected input is
#' an XML-tag-decorated string corresponding to a single PubMed
#' record. Information about article title, authors, affiliations, 
#' journal name and abbreviation, publication date, references, and
#' keywords are returned.  
#' 
#' 
#' @param pubmedArticle String, this is an XML-tag-decorated raw PubMed record. 
#' @param max_authors Numeric, maximum number of authors to retrieve. If this
#' is set to -1, only the last author is extracted. If this is set to 1, 
#' only the first author is returned. If this is set to 2, the first and the 
#' last authors are extracted. If this is set to any other positive 
#' number (i), up to the leading (n-1) authors are retrieved together with the 
#' last author. If this is set to a number larger than the number of authors in
#' a record, all authors are returned. Note that at least 1 author has to be
#' retrieved, therefore a value of 0 is not accepted (coerced to -1).
#' @param autofill_address Logical, shall author affiliations be 
#' propagated within each record to fill missing values. 
#' @param compact_output Logical, shall record data be returned in a 
#' compact format where each row is a single record and author names are
#' collapsed together. If `FALSE`, each row corresponds to a single author of
#' the publication and the record-specific data are recycled for all included
#' authors.
#' @param include_abstract Logical, shall abstract text be included in the 
#' output data.frame. If `FALSE`, the abstract text column is populated 
#' with a missing value.
#' @param max_references Numeric, maximum number of references to return (for
#' each PubMed record).
#' @param ref_id_type String, must be one of the
#' following values: `c('pmid', 'doi')`. Type of identifier used to describe 
#' citation references. 
#' 
#' 
#' 
#' 
#' @examples 
#' data(epm_samples)
#' x <- epm_samples$bladder_cancer_2018$demo_data_03$raw[[1]]
#' epm_parse_record(x)
#' 
#' 
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return a data.frame including information extracted from a 
#' raw XML PubMed record.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @export
epm_parse_record <- function(pubmedArticle, 
                             max_authors = 15, 
                             autofill_address = TRUE, 
                             compact_output = TRUE,
                             include_abstract = TRUE,
                             max_references = 1000,
                             ref_id_type = "pmid") {

  stopifnot(is.character(pubmedArticle), length(pubmedArticle) == 1, 
            !is.na(pubmedArticle), nchar(pubmedArticle) > 1, 
            is.numeric(max_authors), length(max_authors) == 1, 
            !is.na(max_authors))
  
  # Record - Citation Data & sub-chunks
  x100 <- EPM_custom_grep(xml_data = pubmedArticle, 
                          tag = 'MedlineCitation', format = 'char')

  x105 <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'Journal', format = 'char')}, 
    error = function(e) { NA })
  
  x115 <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'AuthorList', format = 'char')}, 
    error = function(e) { NA })

  x125 <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'GrantList', format = 'char')}, 
    error = function(e) { NA })

  # x135 <- tryCatch({
  #   EPM_custom_grep(xml_data = x100, tag = 'PublicationTypeList', 
  #                   format = 'char')}, 
  #   error = function(e) { NA })

  x145 <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'MeshHeadingList', format = 'char')}, 
    error = function(e) { NA })
  
  
    
  # Record - References & sub-chunks
  x200 <- tryCatch({
    EPM_custom_grep(xml_data = pubmedArticle, 
                          tag = 'PubmedData', format = 'char')},
    error = function(e) { NULL })
  
  x210 <- tryCatch({
    tmp0 <- NULL
    tmp0 <- strsplit(x = x200, split = '<ReferenceList>', fixed = TRUE)
    if (is.list(tmp0) && length(tmp0) == 1 && length(tmp0[[1]]) == 2) {
      list(
        record = tryCatch({
          tmp1 <- EPM_custom_grep(xml_data = tmp0[[1]][1], 
                                  tag = 'ArticleIdList', format = 'char')
        }, error = function(e) { NA }), 
        references = paste0('<ReferenceList>', tmp0[[1]][2]))} 
    
    }, error = function(e) { list() })
  

  ## --- Section 1
  
  # PMID
  pmid <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'PMID', format = 'char')}, 
    error = function(e) { NA })
  
  # Dates
  date_cmpl <- tryCatch({
    tmp_d0 <- NULL
    tmp_d0 <- EPM_custom_grep(xml_data = x100, 
                              tag = 'DateCompleted', format = 'char')
    EPM_date_parse(tmp_d0)
  }, error = function(e) { NULL })
  
  date_rvsd <- tryCatch({
    tmp_d0 <- NULL
    tmp_d0 <- EPM_custom_grep(xml_data = x100, 
                              tag = 'DateRevised', format = 'char')
    EPM_date_parse(tmp_d0)
  }, error = function(e) { NULL })
  
  date_publ <- tryCatch({
    tmp_d0 <- NULL
    tmp_d0 <- EPM_custom_grep(xml_data = x105, 
                              tag = 'PubDate', format = 'char')
    EPM_date_parse(tmp_d0)
  }, error = function(e) { NULL })
  
  # Converge on a date
  if (!is.null(date_publ)) {
    date_out <- date_publ
  } else if (!is.null(date_cmpl)) {
    date_out <- date_cmpl
  } else if (!is.null(date_rvsd)) {
    date_out <- date_rvsd
  } else {
    date_out <- list()
  }

    
  # Journal Name and Abbrv
  jtitle <- tryCatch({
    EPM_custom_grep(xml_data = x105, tag = 'Title', format = 'char')}, 
    error = function(e) { NA })
  
  jabbrv <- tryCatch({
    EPM_custom_grep(xml_data = x105, tag = 'ISOAbbreviation', format = 'char')}, 
    error = function(e) { NA })
  
  artdoi <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'ELocationID', 
                    xclass = 'doi', format = 'char')}, 
    error = function(e) { NA })

  artpmc <- tryCatch({
    EPM_custom_grep(xml_data = x210$record, tag = 'ArticleId', 
                    xclass = 'pmc', format = 'char')}, 
    error = function(e) { NA })
  
  arttitle <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'ArticleTitle', format = 'char')}, 
    error = function(e) { NA })
  
  artabstr <- NA
  if (include_abstract) {
    artabstr <- tryCatch({
      EPM_custom_grep(xml_data = x100, tag = 'AbstractText', format = 'char')}, 
      error = function(e) { NA })
  } 
  
  artlang <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'Language', format = 'char')}, 
    error = function(e) { NA })
  
  artcoi <- tryCatch({
    EPM_custom_grep(xml_data = x100, tag = 'CoiStatement', format = 'char')}, 
    error = function(e) { NA })
  
  grantids <- tryCatch({
    EPM_custom_grep(xml_data = x125, tag = 'GrantID', format = 'char')}, 
    error = function(e) { NA })
  
  mesh_trms <- tryCatch({
    EPM_mesh_parse(x = x145)}, error = function(e) { list()})
  
  refes <- tryCatch({
    EPM_reference_parse(x = x210$references, 
                        max_references = max_references,
                        id_type = ref_id_type)
  }, error = function(e) { NA })
  
  artauths <- tryCatch({
    EPM_auth_parse(x = x115, max_authors = max_authors, 
                   autofill = autofill_address) 
    }, error = function(e) { 
      data.frame(lastname=NA, forename=NA,
                 address=NA, email=NA, 
                 stringsAsFactors = FALSE) 
    })
  
  #
  # Put things together, prep for output
  if (compact_output) {
    
    Y <- data.frame(
      pmid = ifelse(!is.null(pmid), pmid, NA),
      doi = ifelse(!is.null(artdoi), artdoi, NA), 
      pmc = ifelse(!is.null(artpmc), artpmc, NA), 
      journal = ifelse(!is.null(jtitle), jtitle, NA), 
      jabbrv = ifelse(!is.null(jabbrv), jabbrv, NA), 
      lang = ifelse(!is.null(artlang), artlang, NA),
      year = ifelse(!is.null(date_out$year), date_out$year, NA),
      month = ifelse(!is.null(date_out$month), date_out$month, NA), 
      day = ifelse(!is.null(date_out$day), date_out$day, NA), 
      title = ifelse(!is.null(arttitle), arttitle, NA), 
      
      abstract = ifelse(!is.null(artabstr) && 
                          length(artabstr)>0 && 
                          !is.na(artabstr[[1]]), 
                        paste(artabstr, collapse = ' '), NA),
      mesh_codes = ifelse(!is.null(mesh_trms$mesh_codes) &&
                            length(mesh_trms$mesh_codes) > 0 &&
                            !is.na(mesh_trms$mesh_codes[[1]]), 
                          paste(mesh_trms$mesh_codes, collapse = '; '), NA),
      mesh_terms = ifelse(!is.null(mesh_trms$mesh_terms) &&
                            length(mesh_trms$mesh_terms) > 0 &&
                            !is.na(mesh_trms$mesh_terms[[1]]), 
                          paste(mesh_trms$mesh_terms, collapse = '; '), NA),
      grant_ids = ifelse(!is.null(grantids) && 
                           length(grantids) > 0 && 
                           !is.na(grantids[[1]]), 
                         paste(grantids, collapse = '; '), NA), 
      references = ifelse(!is.null(refes) && 
                            length(refes) > 0 && 
                            !is.na(refes[[1]]), 
                          paste(refes, collapse = '; '), NA),
      coi = ifelse(!is.null(artcoi) && 
                     length(artcoi) > 0 && 
                     !is.na(artcoi[[1]]), 
                   paste(artcoi, collapse = '; '), NA),
      
      authors = ifelse(!is.null(artauths$collapsed$authors),
                       artauths$collapsed$authors, NA),
      affiliation = ifelse(!is.null(artauths$collapsed$address),
                           artauths$collapsed$address, NA),
      stringsAsFactors = FALSE)
    
  } else {
    
    if (is.null(artauths$authors)) {
      
      Y <- NULL
      
    } else {
      
      Y <- data.frame(
        pmid = ifelse(!is.null(pmid), pmid, NA),
        doi = ifelse(!is.null(artdoi), artdoi, NA), 
        pmc = ifelse(!is.null(artpmc), artpmc, NA), 
        journal = ifelse(!is.null(jtitle), jtitle, NA), 
        jabbrv = ifelse(!is.null(jabbrv), jabbrv, NA), 
        lang = ifelse(!is.null(artlang), artlang, NA),
        year = ifelse(!is.null(date_out$year), date_out$year, NA),
        month = ifelse(!is.null(date_out$month), date_out$month, NA), 
        day = ifelse(!is.null(date_out$day), date_out$day, NA), 
        title = ifelse(!is.null(arttitle), arttitle, NA), 
        
        abstract = ifelse(!is.null(artabstr) && 
                            length(artabstr)>0 && 
                            !is.na(artabstr[[1]]), 
                          paste(artabstr, collapse = ' '), NA),
        mesh_codes = ifelse(!is.null(mesh_trms$mesh_codes) &&
                              length(mesh_trms$mesh_codes) > 0 &&
                              !is.na(mesh_trms$mesh_codes[[1]]), 
                            paste(mesh_trms$mesh_codes, collapse = '; '), NA),
        mesh_terms = ifelse(!is.null(mesh_trms$mesh_terms) &&
                              length(mesh_trms$mesh_terms) > 0 &&
                              !is.na(mesh_trms$mesh_terms[[1]]), 
                            paste(mesh_trms$mesh_terms, collapse = '; '), NA),
        grant_ids = ifelse(!is.null(grantids) && 
                             length(grantids) > 0 && 
                             !is.na(grantids[[1]]), 
                           paste(grantids, collapse = '; '), NA), 
        references = ifelse(!is.null(refes) && 
                              length(refes) > 0 && 
                              !is.na(refes[[1]]), 
                            paste(refes, collapse = '; '), NA),
        coi = ifelse(!is.null(artcoi) && 
                       length(artcoi) > 0 && 
                       !is.na(artcoi[[1]]), 
                     paste(artcoi, collapse = '; '), NA),
        
        last_name = artauths$authors$lastname, 
        first_name = artauths$authors$forename,
        affiliation = artauths$authors$address, 
        email = artauths$authors$email,
        stringsAsFactors = FALSE)
    }
  }
  
  # Out
  tryCatch({rownames(Y) <- NULL}, error = function(e) { NULL })
  return(Y)
}







#' Fetch Raw Records from Pubmed.
#' 
#' Fetch raw PubMed records from PubMed. Records can be downloaded 
#' in text or xml format and stored into a local object or written to 
#' local files. 
#' 
#' 
#' @param x An `easyPubMed` object. 
#' @param format String, the desired format for the raw records. 
#' This argument must take one of the following
#' values: `c("uilist", "medline", "xml")` and defaults to `"xml"`.
#' @param api_key String, corresponding to the NCBI API token (if available). 
#' NCBI token strings can be requested from NCBI. Record download will be 
#' faster if a valid NCBI token is used. This argument can be `NULL`. 
#' @param write_to_file Logical of length 1. Shall raw records be written to 
#' a file on the local machine. It defaults to `FALSE`.
#' @param outfile_path Path to the folder on the local machine where files 
#' will be saved (if `write_to_file` is `TRUE`). It must point to an
#' already existing directory. If `NULL`, the working directory will be used. 
#' @param outfile_prefix String, prefix that will be added to the name
#' of each file written to the local machine. This argument is parsed only 
#' when `write_to_file` is `TRUE`. If `NULL`, an arbitrary prefix will be added 
#' (easypubmed_job_YYYYMMDDHHMM). 
#' @param store_contents Logical of length 1. Shall raw records be stored
#' in the `easyPubMed` object. It defaults to `TRUE`. It may convenient to 
#' switch this to `FALSE` when downloading large number of records. 
#' If `store_contents` is `FALSE`, `write_to_file` must be `TRUE`. 
#' @param encoding String, the encoding of the records retrieved from PubMed. 
#' Typically, this is 'UTF-8'. 
#' @param verbose Logical, shall details about the 
#' progress of the operation be printed to console.
#' 
#' 
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   x <- epm_fetch(x = x, format = 'uilist')
#'   x
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return an easyPubMed object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @export
epm_fetch <- function(x, format = 'xml',
                      api_key = NULL, 
                      write_to_file = FALSE, 
                      outfile_path = NULL, 
                      outfile_prefix = NULL,
                      store_contents = TRUE, 
                      encoding = 'UTF-8', 
                      verbose = TRUE) {
  
  # x must be an easyPubMed object
  stopifnot(inherits(x = x, what = "easyPubMed"))
  
  # validate all other params
  par_list <- list(format = format, api_key = api_key, 
                   write_to_file = write_to_file, outfile_path = outfile_path, 
                   outfile_prefix = outfile_prefix,
                   store_contents = store_contents, encoding = encoding, 
                   verbose = verbose)
  par_list <- EPM_validate_fetch_params(par_list)
  
  # execute job and return
  y <- fetchEPMData(x = x, params = par_list)

  return(y)
}











#' Extract Information from a Raw PubMed Record.
#' 
#' Read a raw PubMed record, identify XML tags, extract information 
#' and cast it into a structured data.frame. The expected input is
#' an XML-tag-decorated string corresponding to a single PubMed
#' record. Information about article title, authors, affiliations, 
#' journal name and abbreviation, publication date, references, and
#' keywords are returned.  
#' 
#' 
#' @param x An `easyPubMed` object. The object must include raw records (n>0) 
#' downloaded in the 'xml' format.
#' @param max_authors Numeric, maximum number of authors to retrieve. If this
#' is set to -1, only the last author is extracted. If this is set to 1, 
#' only the first author is returned. If this is set to 2, the first and the 
#' last authors are extracted. If this is set to any other positive 
#' number (i), up to the leading (n-1) authors are retrieved together with the 
#' last author. If this is set to a number larger than the number of authors in
#' a record, all authors are returned. Note that at least 1 author has to be
#' retrieved, therefore a value of 0 is not accepted (coerced to -1).
#' @param autofill_address Logical, shall author affiliations be 
#' propagated within each record to fill missing values. 
#' @param compact_output Logical, shall record data be returned in a 
#' compact format where each row is a single record and author names are
#' collapsed together. If `FALSE`, each row corresponds to a single author of
#' the publication and the record-specific data are recycled for all included
#' authors (legacy approach).
#' @param include_abstract Logical, shall abstract text be included in the 
#' output data.frame. If `FALSE`, the abstract text column is populated 
#' with a missing value.
#' @param max_references Numeric, maximum number of references to return (for
#' each PubMed record).
#' @param ref_id_type String, must be one of the
#' following values: `c('pmid', 'doi')`. Type of identifier used to describe 
#' citation references. 
#' @param verbose Logical, shall details about the 
#' progress of the operation be printed to console.
#' 
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   x <- epm_fetch(x = x, format = 'xml')
#'   x <- epm_parse(x, include_abstract = FALSE, max_authors = 1)
#'   get_epm_data(x)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#'  
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return an easyPubMed object including a data.frame (`data` slot) that 
#' stores information extracted from its raw XML PubMed records.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @export
epm_parse <- function(x, max_authors = 10, 
                      autofill_address = TRUE, 
                      compact_output = TRUE, 
                      include_abstract = TRUE, 
                      max_references = 150, 
                      ref_id_type = 'doi', 
                      verbose = TRUE) {
  
  # x must be an easyPubMed object
  stopifnot(inherits(x = x, what = "easyPubMed"))
  
  # read mets and check if everything looks OK
  my_meta <- getEPMMeta(x)
  stopifnot(!is.null(my_meta$raw_format), 
            !is.na(my_meta$raw_format), 
            my_meta$raw_format[1] == 'xml', 
            
            !is.null(my_meta$raw_data_embedded), 
            !is.na(my_meta$raw_data_embedded), 
            my_meta$raw_data_embedded == 1, 
            
            !is.null(my_meta$raw_record_num), 
            !is.na(my_meta$raw_record_num), 
            my_meta$raw_record_num > 0)
  
  # validate all other params
  par_list <- list(max_authors = max_authors, 
                   autofill_address = autofill_address, 
                   compact_output = compact_output, 
                   include_abstract = include_abstract, 
                   max_references = max_references, 
                   ref_id_type = ref_id_type, 
                   verbose = verbose)
  par_list <- EPM_validate_parse_params(par_list)
  
  
  # execute job and return
  y <- parseEPMData(x = x, params = par_list)
  return(y)
}




#' Import PubMed Records from Local Files.
#' 
#' Read one or more text files including XML-decorated raw PubMed records 
#' and rebuild an `easyPubMed` object. The function expects all files to be
#' generated from the same query using `easyPubMed>3.0` and the 
#' `epm_fetch()` function setting `write_to_file` to `TRUE`. This 
#' function can import a fraction or all of the files 
#' resulting from a single query. Files resulting
#' from non-compatible fetch jobs will be dropped. 
#' 
#' 
#' @param x Character vector, the paths to text files including 
#' XML-decorated raw PubMed records saved using `easyPubMed>3.0`.
#' 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   x <- epm_fetch(x = x, format = 'xml', write_to_file = TRUE, 
#'                  outfile_prefix = 'test', store_contents = FALSE)
#'   y <- epm_import_xml('test_batch_01.txt')
#'   tryCatch({unlink('test_batch_01.txt')}, error = function(e) { NULL }) 
#'   print(paste0('       Raw Record Num (fetched): ', 
#'                getEPMMeta(x)$raw_record_num))
#'   print(paste0('Raw Record Num (read & rebuilt): ', 
#'                getEPMMeta(y)$raw_record_num))
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#'  
#' 
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return an `easyPubMed` object including raw XML PubMed records.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'  
#'
#' 
#' @export
epm_import_xml <- function(x) {
  
  # check if file exist
  stopifnot(is.character(x), length(x) > 0, 
            sum(is.na(x)) == 0, 
            sum(vapply(x, file.exists, FUN.VALUE = TRUE)) == length(x))
  
  # initialize collector
  TMP <- list(guide = list(), meta = list(), 
              job_info = list(), contents = list())
  
  inp_cnt <- 0
  
  # loop and fetch data
  for (i in seq_len(length(x))) {
    
    xi <- x[i]
    tryCatch({
      tmp_i <- EPM_read_xml(x = xi)
      
      # rebuild guide
      tmp_gui <- data.frame(index = i, 
                            file = xi, 
                            JobUniqueId = tmp_i$guide$JobUniqueId,
                            JobQuery = tmp_i$guide$JobQuery, 
                            JobBatch = tmp_i$guide$JobBatch,
                            JobBatchNum = tmp_i$guide$JobBatchNum,
                            stringsAsFactors = FALSE)
      
      # Add to TMP
      TMP$guide[[length(TMP$guide) + 1]] <- tmp_gui
      TMP$meta[[length(TMP$meta) + 1]] <- tmp_i$meta
      TMP$job_info[[length(TMP$job_info) + 1]] <- tmp_i$job_info
      TMP$contents[[length(TMP$contents) + 1]] <- tmp_i$contents
      
      inp_cnt <- inp_cnt + 1
    }, error = function(e) { NULL })
  }
  
  # Check guide, are there problems?
  final_guide <- do.call(rbind, TMP$guide)
  final_guide <- EPM_check_guide(x = final_guide)
  
  # Filter based on final guide check
  if (sum(final_guide$pass == 1) < 1) {
    stop(paste0('None of the files could be imported. Were these files ', 
                'written in XML format using `easyPubMed>=3.0` ?'))
  }
  
  if (sum(final_guide$pass == 0) > 0) {
    wmsg <- paste0('Some of the XML files were discarded because of ', 
                   'query incompatibility or potential record duplication. ', 
                   'Please, import the following files separately: "', 
                   paste(final_guide$file[final_guide$pass == 0], 
                         collapse = '", "'), '".')
    rlang::warn(message = wmsg)
  }
  if (sum(final_guide$pass == 1) !=  
      final_guide$JobBatchNum[final_guide$pass == 1][1]) {
    wmsg <- paste0('Note that only part of the records returned by the ', 
                   'original query were imported.')
    rlang::warn(message = wmsg)
  }
  
  # Put things together & build an easyPubMed object.
  keep <- final_guide$pass == 1
  new_query_gui <- do.call(rbind, TMP$job_info)[keep, ]
  new_query_str <- new_query_gui$query_string[1]
  new_meta <- TMP$meta[[  which.min(keep) ]]
  
  new_job_list <- list(meta = new_meta, query_guide = new_query_gui)
  
  new_contents <- list()
  for (i in which(keep)) {
    if (!is.null(names(TMP$contents[[i]]))) {
      for (j in names(TMP$contents[[i]])) {
        new_contents[[j]] <- TMP$contents[[i]][[j]]     
      }
    } else {
      for (j in seq_along(TMP$contents[[i]])) {
        new_contents[[length(new_contents) + 1]] <- TMP$contents[[i]][[j]]     
      }
    }
  }
  
  if (!is.null(names(new_contents))) {
    new_pmids <- names(new_contents)
  } else {
    new_pmids <- tryCatch({EPM_detect_pmid(new_contents, as.list = FALSE)}, 
                          error = function(e) { NULL })
  }
  
  # Check records
  new_num_recs <- length(new_contents)
  new_job_list$meta$raw_record_num <- new_num_recs
  raw_n_diff <- new_job_list$meta$exp_count - new_num_recs
  if (new_job_list$meta$exp_count > new_num_recs) {
    new_job_list$meta$all_records_covered <- FALSE
    new_job_list$meta$exp_missed_records  <- raw_n_diff
  } else if (new_job_list$meta$exp_count < new_num_recs) {
    new_job_list$meta$all_records_covered <- NA
    new_job_list$meta$exp_missed_records  <- raw_n_diff
    wmsg <- paste0('The XML files you just imported had more PubMed ', 
                   'records than expected. This is an odd exception!\n',
                   'Please, be careful if you decide to proceed, since this ', 
                   'easyPubMed object may include inaccurate data...')
    rlang::warn(message = wmsg)
  }
  
  # New Misc
  new_misc <- list(job_list = new_job_list$query_guide, 
                   import_xml_epm_params = list(
                     encoding = new_job_list$meta$raw_encoding, 
                     format = new_job_list$meta$raw_format, 
                     input_files = x,
                     input_files_pass = keep, 
                     expected_record_num = new_job_list$meta$exp_count, 
                     actual_record_num = new_job_list$meta$raw_record_num, 
                     record_num_diff = raw_n_diff,
                     record_pmid_avail = length(new_pmids) > 0))
  
  # Build object
  y <- new(Class = 'easyPubMed', 
           query_string = new_query_str, job_info = new_job_list)
  
  # update contents
  y <- setEPMMeta(y, new_job_list$meta)
  y <- setEPMMisc(y, new_misc)
  tryCatch({
    if (is.list(new_pmids)) {
      y <- setEPMUilist(x = y, y = new_pmids)
    } else if (is.vector(new_pmids)) {
      y <- setEPMUilist(x = y, y = lapply(new_pmids, function(i) { i }))
    }
  }, error = function(e) { NULL })
  y <- setEPMRaw(y, new_contents)
  
  # return(y)
  return(y)  
}









# ----




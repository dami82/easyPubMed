#' @title Extract Data from a PubMed Record
#'
#' @description Extract publication-specific information from a PubMed record driven by XML tags. 
#' The input record is a string (character-class vector of length 1) and includes 
#' PubMed-specific XML tags. Data are returned as a data frame where each row corresponds 
#' to one of the authors of the PubMed article.
#'
#' @usage article_to_df(pubmedArticle, autofill = FALSE, 
#'                      max_chars = 500, getKeywords = FALSE, 
#'                      getAuthors = TRUE)
#'                     
#' @param pubmedArticle String including one PubMed record.
#' @param autofill Logical. If TRUE, missing affiliations are automatically imputed based on other non-NA 
#' addresses from the same record.
#' @param max_chars Numeric (integer). Maximum number of characters to be extracted from the Article 
#' Abstract field. Set max_chars to -1 for extracting the full-length abstract. Set max_chars to 0 to 
#' extract no abstract.
#' @param getKeywords Logical. If TRUE, an attempt to extract article Keywords will be made.
#' @param getAuthors Logical. If FALSE, author information won't be extracted. This will considerably 
#' speed up the operation.
#'
#' @details 
#' Given one Pubmed Article record, this function will automatically extract a set of features. 
#' Extracted information include: PMID, DOI, article title, article abstract, publication date (year, month, day), 
#' journal name (title, abbreviation), keywords, and a set of author-specific info (names, affiliation, email address). 
#' Each row of the output data frame corresponds to one of the authors of the PubMed record. Author-independent info 
#' (publication ID, title, journal, date) are identical across all rows. If information about authors are not required, 
#' set 'getAuthors' = TRUE.
#'
#' @return Data frame including the extracted features. Each row correspond a different author.
#'
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' try({
#'   ## Display some contents
#'   data("EPMsamples")
#'   #display Query String used for collecting the data
#'   print(EPMsamples$NUBL_1618$qry_st)
#'   #Get records
#'   BL_list <- EPMsamples$NUBL_1618$rec_lst
#'   cat(BL_list[[1]])
#'   # cast PM recort to data.frame
#'   BL_df <- article_to_df(BL_list[[1]], max_chars = 0)
#'   print(BL_df)
#' }, silent = TRUE)
#' 
#' \dontrun{
#' ## Query PubMed, retrieve a selected citation and format it as a data frame
#' dami_query <- "Damiano Fantini[AU] AND 2017[PDAT]"
#' dami_on_pubmed <- get_pubmed_ids(dami_query)
#' dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed)
#' dami_abstracts_list <- articles_to_list(dami_abstracts_xml)
#' article_to_df(pubmedArticle = dami_abstracts_list[[1]], autofill = FALSE)
#' article_to_df(pubmedArticle = dami_abstracts_list[[2]], autofill = TRUE, max_chars = 300)[1:2,]
#' }
#'
#' @export
article_to_df <- function (pubmedArticle, autofill = FALSE, 
                           max_chars = 500, getKeywords = FALSE, 
                           getAuthors = TRUE) 
{
  currWarn <- options()$warn
  options(warn = -1)
  fix_pdate <- function(date_vec) {
    
    month_dict <- data.frame(id = c("jan", "january", "01", "1", 
                                    "feb", "february", "02", "2", 
                                    "mar", "march", "03", "3", 
                                    "apr", "april", "04", "4", 
                                    "mar", "march", "05", "5", 
                                    "jun", "june", "06", "6", 
                                    "jul", "july", "07", "7", 
                                    "aug", "august", "08", "8", 
                                    "sep", "september", "09", "9", 
                                    "oct", "october", "10", 
                                    "nov", "november", "11", 
                                    "dec", "december", "12"), 
                             val = c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4),
                                     rep(5, 4), rep(6, 4), rep(7, 4), rep(8, 4),
                                     rep(9, 4), rep(10, 3), rep(11, 3), rep(12, 3)))
    
    date_vec <- tolower(date_vec)
    
    if (!is.na(date_vec[[1]])) {
      month <- tryCatch({
        match.arg(arg = tolower(date_vec[[2]]), choices = month_dict$id, 
                  several.ok = FALSE)
      }, error = function(e) {
        "january"
      })
      month <- month_dict$val[month_dict$id == month]
    } else {
      month <- NA
    }
    if (is.na(date_vec[[3]]) && !is.na(month) && !is.na(date_vec[[1]])) {
      date_vec[[3]] <- 1
    }
    
    # attempt getting MedlineDate Year
    if (is.na(date_vec[[1]])) {
      if (!is.na(date_vec[[4]])) {
        nuYR <- sub("^.*([[:digit:]]{4}).*$", "\\1", date_vec[[4]])
        if (nchar(nuYR) == 4 && grepl("(19)|(20)[[:digit:]]{2}", nuYR)) {
          date_vec[[1]] <- nuYR 
        }
      }
    }
    
    return(c(Year = as.numeric(date_vec[[1]]), Month = as.numeric(month), 
             Day = as.numeric(date_vec[[3]])))
  }
  
  
  if (class(pubmedArticle)[1] != "character" || length(pubmedArticle) != 
      1 || regexpr("(<PubmedArticle)(.+)(\\/PubmedArticle>)", 
                   pubmedArticle) < 0) {
    message("An error occurred")
    return(NULL)
  }
  if (!is.numeric(max_chars)) {
    max_chars <- 500
  } else if (max_chars < 0) {
    max_chars <- -1
  }
  tryCatch({
    pubmedArticle <- gsub("&amp;", "&", pubmedArticle, ignore.case = TRUE)
    tmp.article <- custom_grep(xml_data = pubmedArticle, 
                               tag = "PubmedArticle", format = "char")
    if (is.null(tmp.article)) {
      message("An error occurred")
      return(NULL)
    }
    tmp.title <- custom_grep(xml_data = tmp.article, tag = "ArticleTitle", 
                             format = "char")
    if (length(tmp.title) > 1) {
      tmp.title <- paste(tmp.title, collapse = " ", sep = " ")
    } else if (length(tmp.title) < 1) {
      tmp.title <- NA
    }
    tmp.abstract <- custom_grep(xml_data = tmp.article, tag = "AbstractText", 
                                format = "char")
    if (length(tmp.abstract) > 1) {
      tmp.abstract <- paste(tmp.abstract, collapse = " ", 
                            sep = " ")
      if (max_chars >= 0) {
        tmp.abstract <- gsub("</{0,1}i>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- gsub("</{0,1}b>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- gsub("</{0,1}sub>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- gsub("</{0,1}exp>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- substr(tmp.abstract, 0, max_chars)
      }
    } else if (length(tmp.abstract) < 1) {
      tmp.abstract <- NA
    } else {
      if (max_chars >= 0) {
        tmp.abstract <- substr(tmp.abstract, 0, max_chars)
        tmp.abstract <- gsub("</{0,1}i>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- gsub("</{0,1}b>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- gsub("</{0,1}sub>", "", tmp.abstract, 
                             ignore.case = T)
        tmp.abstract <- gsub("</{0,1}exp>", "", tmp.abstract, 
                             ignore.case = T)
      }
    }
    my.dateType <- c("PubDate", "DateCompleted", "DateCreated", 
                     "DateRevised")
    sel.dateType <- which(sapply(my.dateType, (function(xi) {
      regexpr(xi, tmp.article) > 0
    })))
    if (length(sel.dateType) < 1) {
      tmp.date <- c(Year = NA, Month = NA, Day = NA)
    } else {
      sel.dateType <- sel.dateType[1]
      tmp.date <- custom_grep(xml_data = tmp.article, tag = my.dateType[sel.dateType], 
                              format = "char")
      
      # Modify -- MEDLINEDATE??
      tmp.date <- sapply(c("Year", "Month", "Day", "MedlineDate"), (function(tt) {
        tdat.el <- custom_grep(xml_data = tmp.date, tag = tt, 
                               format = "char")
        ifelse(is.null(tdat.el), NA, tdat.el[1])
      }))
    }
    tmp.date <- fix_pdate(tmp.date)
    tmp.paperID <- custom_grep(xml_data = tmp.article, tag = "ArticleIdList", 
                               format = "char")
    if (is.null(tmp.paperID)) {
      message("An error occurred")
      return(NULL)
    } else {
      tmp.paperID <- gsub("[[:space:]]", "", tmp.paperID[1])
    }
    tmp.PMID <- gsub("^(.*ArticleIdIdType=\\\"pubmed\\\")([[:space:]]|[[:alnum:]]){0,20}>", 
                     "", tmp.paperID)
    tmp.PMID <- gsub("<.*$", "", tmp.PMID)
    tmp.DOI <- gsub("^(.*ArticleIdIdType=\\\"doi\\\")([[:space:]]|[[:alnum:]]){0,20}>", 
                    "", tmp.paperID)
    tmp.DOI <- gsub("<.*$", "", tmp.DOI)
    tmp.jabbrv <- custom_grep(xml_data = tmp.article, tag = "ISOAbbreviation", 
                              format = "char")
    tmp.jabbrv <- ifelse(is.null(tmp.jabbrv), NA, tmp.jabbrv)
    tmp.journal <- custom_grep(xml_data = tmp.article, tag = "Title", 
                               format = "char")
    tmp.journal <- ifelse(is.null(tmp.journal), NA, tmp.journal)
    tmp.keys <- tryCatch({
      if (getKeywords) {
        tmp.keys <- custom_grep(xml_data = tmp.article, 
                                tag = "Keyword", format = "char")
        tmp.mesh <- custom_grep(xml_data = tmp.article, 
                                tag = "MeshHeading", format = "char")
        if (length(tmp.mesh) > 0) {
          tmp.mesh <- sapply(tmp.mesh, function(xxm) {
            custom_grep(xml_data = xxm, tag = "DescriptorName", 
                        format = "char")
          })
        }
        tmp.keys <- c(tmp.keys, tmp.mesh)
        if (length(tmp.keys) > 1) {
          tmp.keys <- paste(tmp.keys, collapse = "; ")
        } else if (length(tmp.keys) < 1) {
          tmp.keys <- NA
        }
      } else {
        NA
      }
    }, error = function(e) {
      NA
    })
    tmp.resout <- c(pmid = tmp.PMID, doi = tmp.DOI, title = tmp.title, 
                    abstract = tmp.abstract, year = as.vector(tmp.date[1]), 
                    month = as.vector(tmp.date[2]), day = as.vector(tmp.date[3]), 
                    jabbrv = tmp.jabbrv, journal = tmp.journal, keywords = tmp.keys)
    tmp.authors <- custom_grep(xml_data = tmp.article, tag = "AuthorList", 
                               format = "char")
    if (length(tmp.authors) < 1 | !getAuthors) {
      final.mat <- data.frame(rbind(c(tmp.resout, lastname = NA, 
                                      firstname = NA, address = NA, email = NA)), stringsAsFactors = FALSE)
    } else {
      author.list <- custom_grep(xml_data = tmp.authors, 
                                 tag = "Author", format = "char")
      final.mat <- do.call(rbind, lapply(author.list, (function(al) {
        tmp.lastnm <- custom_grep(xml_data = al, tag = "LastName", 
                                  format = "char")
        tmp.firstnm <- custom_grep(xml_data = al, tag = "ForeName", 
                                   format = "char")
        email.PAT <- "([[:alnum:]]|\\.|\\-|\\_){3,200}@([[:alnum:]]|\\.|\\-|\\_){3,200}(\\.)([[:alnum:]]){2,6}"
        tmp.email <- regexpr(email.PAT, al)
        if (tmp.email > 0) {
          tmp.email <- substr(al, tmp.email, tmp.email + 
                                attributes(tmp.email)$match.length - 1)
          al <- gsub(email.PAT, "", al)
        } else {
          tmp.email <- NA
        }
        if (regexpr("Affiliation", al) > 0) {
          tmp.add <- custom_grep(al, "Affiliation", format = "char")[1]
          tmp.add <- trim_address(tmp.add)
        } else {
          tmp.add <- NA
        }
        c(tmp.resout, lastname = tmp.lastnm, firstname = tmp.firstnm, 
          address = tmp.add, email = tmp.email)
      })))
      rownames(final.mat) <- NULL
      final.mat <- data.frame(final.mat, stringsAsFactors = FALSE)
      DESELECT <- is.na(final.mat$lastname) | is.na(final.mat$firstname)
      if (length(DESELECT) > 0 & sum(DESELECT) > 0) 
        final.mat <- final.mat[!DESELECT, ]
      if (autofill) {
        tmp.address <- final.mat[, "address"]
        na.pos <- is.na(tmp.address)
        if (sum(na.pos) != length(tmp.address)) {
          tmp.list <- lapply(tmp.address, function(x) {
            x
          })
          cur.add <- tmp.list[[(which(!na.pos)[1])]]
          for (i in 1:length(na.pos)) {
            if (na.pos[i]) {
              tmp.list[[i]] <- cur.add
            } else {
              cur.add <- tmp.list[[i]]
            }
          }
          final.mat[, "address"] <- do.call(c, tmp.list)
        }
      }
    }
    if (ncol(final.mat) != 14) {
      final.mat <- NULL
    }
  }, error = function(e) {
    NULL
  }, finally = {
    options(warn = currWarn)
    return(final.mat)
  })
}




#' @title Cast PubMed Data into a List of Articles
#'
#' @description Convert an XML object of PubMed records into a list of strings 
#' (character vector of length 1) corresponding to individual PubMed articles. 
#' PubMed records are identified by a "/PubmedArticle" XML tag. This automatically casts 
#' all the content of each PubMed record to a character-class object without removing XML tags.
#' 
#' @usage articles_to_list(pubmed_data, encoding = "UTF8", simplify = TRUE)
#' 
#' @param pubmed_data String corresponding to the name of an XML file (typically, 
#' the result of a batch_pubmed_download() call). Alternatively, a string including 
#' PubMed records with XML tags, such as the object returned by a fetch_pubmed_data() call.
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the 
#' function base::iconv(). See iconv() help page for how to find out more about encodings 
#' that can be used on your platform. "UTF-8" is recommended.
#' @param simplify Logical; should the result be simplified to a character vector. 
#' If FALSE, results are returned as a list.
#' 
#' @details 
#' The input is an XML object or a string including PubMed records (with XML tags). These are the
#' output of easyPubMed functions: fetch_pubmed_data() or batch_pubmed_download(). 
#' The function returns a list or a character vector where each element is a different PubMed record.
#' 
#' @return 
#' List or character vector including all the records from the original XML object in text format. 
#' Elements in the list are not named and are only accessible via their numeric index.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' try({
#'   ## Retrieve PubMed data and return a list ot articles
#'   my_query <- "Damiano Fantini[AU]"
#'   my_query <- get_pubmed_ids(pubmed_query_string = my_query)
#'   my_data <- fetch_pubmed_data(my_query, encoding = "ASCII")
#'   listed_articles <- articles_to_list(my_data)
#'   custom_grep(listed_articles[[2]], "ArticleTitle", "char")
#' }, silent = TRUE)
#' 
#' \dontrun{
#' ## Download PubMed data and return a list ot articles
#' dami_query <- "Damiano Fantini[AU] AND 2018[PDAT]"
#' outfile <- batch_pubmed_download(dami_query, dest_file_prefix = "easyPM_ex001_")
#' listed_articles <- articles_to_list(pubmed_data = outfile)
#' custom_grep(listed_articles[[2]], "ArticleTitle", "char")
#' }
#' 
#' @export
articles_to_list <-
  function(pubmed_data, encoding = "UTF8", simplify = TRUE) 
  {
    
    # Define a nested (core) function handling standardized input
    # This assumes a string as a input
    easyPM_exec_art_to_lst <- function(pm_dataa, simply = TRUE, 
                                       max_pmrec_nchar = 500000) 
    {
      # it's a string to process
      #pm_datab <- strsplit(pm_dataa, "<PubmedArticle(>|[[:space:]]+?.*>)")[[1]][-1]
      
      # New chunk
      pm_datab <- tryCatch({
        strsplit(pm_dataa, "<PubmedArticle(>|[[:space:]]+?.*>)")[[1]][-1]}, 
        error = function(e) { NULL })
  
      if (is.null(pm_datab)) {

        salv_splts <- tryCatch({
          c(1, 
            as.numeric(gregexpr(pattern = "</PubmedArticle>", myXML)[[1]]))}, 
          error = function(e) {NULL})
    
        if (salv_splts[2] > salv_splts[1]) {
          CLCT <- list()
          for (ii in 2:length(salv_splts)) {
            pos_1 <- salv_splts[(ii - 1)]
            pos_2 <- salv_splts[(ii)]
            pos_diff <- (pos_2 - pos_1)
            #message(pos_diff)
        
            if (pos_diff <= max_pmrec_nchar){
              tmp <- substr(myXML, start = pos_1+16, stop = pos_2+15)
              CLCT[[length(CLCT) + 1]] <- tmp
            }
          }
          pm_datab <- do.call(c, CLCT)
        }
      }
      # end of new chunk
      
      pm_datab <- sapply(pm_datab, function(x) {
        #trim extra stuff at the end of the record
        if (!grepl("</PubmedArticle>$", x))
          x <- sub("(^.*</PubmedArticle>).*$", "\\1", x) 
        
        # Rebuid XML structure and proceed
        x <- paste("<PubmedArticle>", x)
        gsub("[[:space:]]{2,}", " ", x)}, 
        USE.NAMES = FALSE, simplify = simply)
      
      pm_datab
    }
    
    # Execute f(x)
    # Handle inputs of different type
    # check if it is a XMLAbstractDocument or a file
    TMP <- substr(pubmed_data[1], 1, 1000)
    if (grepl("<PubmedArticle", TMP)) {
      
      # it's a string to process
      out <- easyPM_exec_art_to_lst(pubmed_data[1], simply = simplify)
      
    } else if (file.exists(pubmed_data[1])) {
      
      # it's a file
      con1 <- file(pubmed_data[1], encoding = "UTF8")
      on.exit(close(con1))
      myXML <- readLines(con = con1, 
                         n = -1, ok = TRUE, encoding = "UTF8") 
      
      if (encoding != "UTF8")
        myXML <- base::iconv(myXML, from = "UTF8", to = encoding, sub = ".")
      
      myXML <- paste(myXML, collapse = "")
      
      # Run as above
      out <- easyPM_exec_art_to_lst(myXML, simply = simplify)
      
    } else {
      message("An error occurred")
      return(NULL)  
    }
    
    return(out)
  }  



#' @title Download PubMed Records in XML or TXT Format
#'
#' @description Performs a PubMed Query (via the get_pubmed_ids() function), downloads the 
#' resulting data (via multiple fetch_pubmed_data() calls) and then saves data in a series of 
#' xml or txt files on the local drive. The function is suitable for downloading 
#' a very large number of records.
#' 
#' @usage batch_pubmed_download(pubmed_query_string, dest_dir = NULL, 
#'                              dest_file_prefix = "easyPubMed_data_", 
#'                              format = "xml", api_key = NULL, 
#'                              batch_size = 400, res_cn = 1, 
#'                              encoding = "UTF8")
#' 
#' @param pubmed_query_string String (character-vector of length 1): this is the string 
#' used for querying PubMed (the standard PubMed Query synthax applies).
#' @param dest_dir String (character-vector of length 1): this string corresponds to the name 
#' of the existing folder where files will be saved. Existing files will be overwritten. 
#' If NULL, the current working directory will be used.
#' @param dest_file_prefix String (character-vector of length 1): this string is used as 
#' prefix for the files that are written locally. 
#' @param format String (character-vector of length 1): data will be requested from Entrez 
#' in this format. Acceptable values are: c("medline","uilist","abstract","asn.1", "xml"). 
#' When format != "xml", data will be saved as text files (txt).
#' @param api_key String (character vector of length 1): user-specific API key to increase 
#' the limit of queries per second. You can obtain your key from NCBI.
#' @param batch_size Integer (1 < batch_size < 5000): maximum number of records 
#' to be saved in a single xml or txt file.
#' @param res_cn Integer (> 0): numeric index of the data batch to start downloading from. 
#' This parameter is useful to resume an incomplete download job after a system crash.
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the 
#' function base::iconv(). See iconv() help page for how to find out more about encodings 
#' that can be used on your platform. Here, we recommend using "UTF-8".
#' 
#' @details 
#' Download large number of PubMed records as a set of xml or txt files that are saved in the 
#' folder specified by the user. This function enforces data integrity. If a batch of downloaded 
#' data is corrupted, it is discarded and downloaded again. Each download cycle is monitored until 
#' the download job is successfully completed. This function should allow to download a whole copy 
#' of PubMed, if desired. The function informs the user about the current progress by constantly 
#' printing to console the number of batches still in queue for download. pubmed_query_string 
#' accepts standard PubMed synthax. The function will query PubMed multiple times using the same 
#' query string. Therefore, it is recommended to use a [EDAT] or a [PDAT] filter in the query 
#' if you want to ensure reproducible results.
#' 
#' @return 
#' Character vector including the names of files downloaded to the local system  
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' \dontrun{
#' ## Example 01: retrieve data from PubMed and save as XML file
#' ml_query <- "Machine Learning[TI] AND 2016[PD]"
#' out1 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
#' readLines(out1[1])[1:30]
#' ##
#' ## Example 02: retrieve data from PubMed and save as TXT file
#' ml_query <- "Machine Learning[TI] AND 2016[PD]"
#' out2 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180, format = "medline")
#' readLines(out2[1])[1:30]
#' }
#' 
#' @export
batch_pubmed_download <-
  function (pubmed_query_string, 
            dest_dir = NULL,
            dest_file_prefix = "easyPubMed_data_",
            format = "xml",
            api_key = NULL,
            batch_size = 400, 
            res_cn = 1, 
            encoding = "UTF8") 
    
  {
    baseDir <- getwd()
    if (!is.null(dest_dir)) {
      setwd(as.character(dest_dir))
    }
    fileName.collector <- list()
    myQuery <- NULL
    my_rtime <- ifelse(is.null(api_key), 0.34, 0.11)
    
    cur_time <- Sys.time()
    while (is.null(myQuery)) {
      diff_time <- my_rtime - as.numeric(difftime(Sys.time(), cur_time, units = "secs"))
      if (diff_time > 0) {
        Sys.sleep(diff_time)
      }
      cur_time <- Sys.time()
      myQuery <- tryCatch(get_pubmed_ids(pubmed_query_string, api_key = api_key), 
                          error = function(e) NULL)
    }
    pubsNum <- as.numeric(myQuery$Count)
    tmpPapers <- NULL
    myRetstart <- 0
    myRetmax <- batch_size
    j = 1
    expTot <- pubsNum/batch_size
    if (expTot > as.integer(expTot)) {
      expTot <- as.integer(expTot) + 1
    } else {
      expTot <- as.integer(expTot)
    }
    while (myRetstart < pubsNum) {
      if (j < res_cn) {
        message(paste("cycle", j, "/", expTot, "skipped...", 
                      sep = " "))
      } else {
        cur_time <- Sys.time()
        while (is.null(myQuery) | is.null(tmpPapers)) {
          
          diff_time <- my_rtime - (as.numeric(Sys.time() - cur_time))
          if (diff_time > 0) {
            Sys.sleep(diff_time)
          }
          cur_time <- Sys.time()
          
          myQuery <- tryCatch(get_pubmed_ids(pubmed_query_string, api_key = api_key), 
                              error = function(e) NULL)
          
          diff_time <- my_rtime - as.numeric(difftime(Sys.time(), cur_time, units = "secs"))
          
          if (diff_time > 0) {
            Sys.sleep(diff_time)
          }
          cur_time <- Sys.time()
          
          # Force download as XML, but withoud collapsing strings
          if (format[1] == "xml") {
            format <- "batchxml"
          }
          
          tmpPapers <- tryCatch(fetch_pubmed_data(pubmed_id_list = myQuery, 
                                                  retstart = myRetstart, 
                                                  retmax = myRetmax,
                                                  format = format, 
                                                  encoding = encoding),
                                
                                error = function(e) NULL, 
                                finally = print(paste("PubMed data batch", 
                                                      j, "/", 
                                                      expTot, "downloaded...",
                                                      sep = " ")))
          if (is.null(tmpPapers)) {
            message("Data retrieval error. Retrying...")
          }
        }
        totDigits <- nchar(as.character(expTot)) + 1
        myExt <- paste(rep(0, totDigits - nchar(as.character(j))), 
                       collapse = "")
        
        tmp.dest.file <- paste(dest_file_prefix, myExt, j, ".txt", sep = "")
        con1 <- file(tmp.dest.file, encoding = encoding)
        doSaveData <- tryCatch(write(tmpPapers, tmp.dest.file), 
                               error = function(e) {"ERROR"}, 
                               finally = {close(con1)})
        if(is.null(doSaveData))
          doSaveData <- tmp.dest.file
        
        myQuery <- NULL
        tmpPapers <- NULL
        if (doSaveData == "ERROR") {
          myRetstart <- myRetstart - myRetmax
          j <- j - 1
          message("An error occurred... Trying to download data from PubMed again...")
        } else {
          fileName.collector[[1+length(fileName.collector)]] <- doSaveData
        }
      }
      myRetstart <- myRetstart + myRetmax
      j <- j + 1
    }
    setwd(baseDir)
    tryCatch(do.call(c, fileName.collector), error = function(e){NULL})
  }


#' @title Retrieve Text Between XML Tags
#'
#' @description Extract text form a string containing XML or HTML tags. Text 
#' included between tags of interest will be returned. If multiple tagged substrings are found, 
#' they will be returned as different elements of a list or character vector.
#' 
#' @usage custom_grep(xml_data, tag, format = "list")
#' 
#' @param xml_data String (of class character and length 1): corresponds to the PubMed 
#' record or any string including XML/HTML tags.
#' @param tag String (of class character and length 1): the tag of interest (does NOT include < > chars).
#' @param format c("list", "char"): specifies the format for the output.
#' 
#' @details 
#' The input string has to be a character string (length 1) containing tags (HTML or XML format). 
#' If an XML Document is provided as input, the function will rise an error.
#' 
#' @return 
#' List or vector where each element corresponds to an in-tag substring.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' try({
#'   ## extract substrings based on regular expressions
#'   string_01 <- "I can't wait to watch the <strong>Late Night Show with" 
#'   string_01 <- paste(string_01, "Seth Meyers</strong> tonight at <strong>11:30</strong>pm CT!")
#'   print(string_01)
#'   custom_grep(xml_data = string_01, tag = "strong", format = "char")
#'   custom_grep(xml_data = string_01, tag = "strong", format = "list")
#' }, silent = TRUE)
#' 
#' @export
custom_grep <-
  function(xml_data, 
           tag, 
           format = "list")
  {
    x <- xml_data[[1]]
    tag.op <- paste("\\<", tag, "((\\>)|([[:space:]]([^[<]]*)\\>))", sep = "")
    tag.cl <- paste("(<\\/)", tag, "(\\>)", sep = "")
    #
    out.result <- list()
    i = 1
    while (!is.null(x) &&
           !is.na(x) &&
           x != "" &&
           nchar(x) > 0 &&
           regexpr(tag.op, x) > 0 &&
           regexpr(tag.cl, x) > 0){
      tag.op.pos <- regexpr(tag.op, x)
      nu.x <- substr(x, (tag.op.pos - 1), nchar(x))
      inner.trim <- regexpr(">", nu.x, fixed = TRUE)
      nu.x <- substr(nu.x, (inner.trim + 1), nchar(nu.x))
      #
      tag.cl.pos <- regexpr(tag.cl, nu.x)
      tag.cl.full <- tag.cl.pos + attributes(tag.cl.pos)$match.length + 1
      x <- substr(nu.x, tag.cl.full, nchar(x))
      nu.x <- substr(nu.x, 1, (tag.cl.pos - 1))
      #
      out.result[[i]] <- nu.x
      i <- i + 1
    }
    if (format != "list") {
      out.result <- do.call(c, out.result)
    }
    return(out.result)
  }


#' @title Retrieve All PubMed Record Identifiers Returned by a Query
#'
#' @description Retrieve PubMed record identifiers from Entrez following a search performed 
#' via the get_pubmed_ids() function. Identifiers are returned as a character vector.
#' 
#' @usage fetch_all_pubmed_ids(pubmed_id_list)
#' 
#' @param pubmed_id_list List: the result of a get_pubmed_ids() call.
#' 
#' @details 
#' Retrieve PubMed identifiers, without any other information (such as article title, 
#' authors, publication date, and so on). The PubMed IDs can be stored or used with other software.
#' 
#' @return 
#' Character vector including all PMID (PubMed Identifiers) returned by the current query.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' \dontrun{
#' ## Fetch only PubMed Record IDs (PMIDs)
#' dami_query_string <- "Damiano Fantini[AU]"
#' dami_on_pubmed <- get_pubmed_ids(dami_query_string)
#' dami_pmids <- fetch_all_pubmed_ids(dami_on_pubmed)
#' print(dami_pmids)
#' 
#' }
#' 
#' @export
fetch_all_pubmed_ids <-
  function(pubmed_id_list)
  {
    # expected records, set retmax
    exp_num <- as.numeric(pubmed_id_list$Count)
    if (is.numeric(exp_num) && exp_num > 0) {
      my_retmax <- exp_num + 1
    } else {
      my_retmax <- 100000
    }
    
    # query, and then extract IDs
    myPubmedURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
                         "db=pubmed&retmax=", my_retmax, "&term=", pubmed_id_list$OriginalQuery, "&usehistory=n", sep = "")
    IDconnect <- url(myPubmedURL, open = "rb", encoding = "UTF8")
    on.exit(close(IDconnect))
    idXML <- readLines(IDconnect, warn = FALSE, encoding = "UTF8") 
    
    collect_ids <- list()
    for (i in 1:length(idXML)) {
      if (grepl("^(\\\\t){0,1}(\\t){0,1}<Id>", idXML[i])) {
        xx <- custom_grep(idXML[i], tag = "Id", format = "char")
        collect_ids[[length(collect_ids) + 1]] <- as.character(xx[1])
      }
    }
    myIDlist <- as.character(do.call(c, collect_ids))
    
    # final check and return
    if(length(myIDlist) != exp_num)
      message(paste("Note that only ", length(myIDlist), " PubMed IDs were retrieved (", 
                    exp_num, " were expected).", sep = ""))
    
    return(myIDlist)
  }


#' @title Retrieve PubMed Data in XML or TXT Format
#'
#' @description Retrieve PubMed records from Entrez following a search performed via the 
#' get_pubmed_ids() function. Data are downloaded in the XML or TXT format and are 
#' retrieved in batches of up to 5000 records.
#' 
#' @usage fetch_pubmed_data(pubmed_id_list, 
#'                          retstart = 0, 
#'                          retmax = 500, 
#'                          format = "xml", 
#'                          encoding = "UTF8")
#' 
#' @param pubmed_id_list List: the result of a get_pubmed_ids() call.
#' @param retstart Integer (>=0): index of the first UID in the retrieved PubMed Search Result set 
#' to be included in the output (default=0, corresponding to the first record of the entire set).
#' @param retmax Integer (>=1): size of the batch of PubMed records to be retrieved at one time.
#' @param format Character: element specifying the output format. The following values are allowed: 
#' c("asn.1", "xml", "medline", "uilist", "abstract").
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the function base::iconv(). 
#' See iconv() help page for how to find out more about encodings that can be used on your platform. 
#' Here, we recommend using "UTF-8".
#' 
#' @details 
#' Retrieve PubMed records based on the results of a get_pubmed_ids() query. 
#' Records are retrieved from Entrez via the PubMed API efetch function. The first entry to be retrieved 
#' may be adjusted via the retastart parameter (this allows the user to download large batches of PubMed 
#' data in multiple runs). The maximum number of entries to be retrieved can also be set adjusting the 
#' retmax parameter (1 < retmax < 5000). Data will be downloaded on the fly (no files are saved 
#' locally).
#' 
#' @return 
#' An object (vector) of class "character". If format is set to "xml" (default), a single String including all 
#' PubMed records (with XML tags embedded) is returned. If a different format is selected, a vector of strings 
#' is returned, where each row corresponds to a line of the output document.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/}
#'
#' @examples
#' try({ 
#'   ## Example 01: retrieve data in TXT format
#'   library("easyPubMed")
#'   dami_query_string <- "Damiano Fantini[AU] AND 2018[PDAT]"
#'   dami_on_pubmed <- get_pubmed_ids(dami_query_string)
#'   Sys.sleep(1) # avoid server timeout
#'   dami_papers <- fetch_pubmed_data(dami_on_pubmed, format = "abstract")
#'   dami_papers[dami_papers == ""] <- "\n"
#'   cat(paste(dami_papers[1:65], collapse = ""))
#' }, silent = TRUE)
#' 
#' \dontrun{
#' ## Example 02: retrieve data in XML format
#' library("easyPubMed")
#' dami_query_string <- "Damiano Fantini[AU]"
#' dami_on_pubmed <- get_pubmed_ids(dami_query_string)
#' dami_papers <- fetch_pubmed_data(dami_on_pubmed)
#' titles <- custom_grep(dami_papers, "ArticleTitle", "char")
#' print(titles)
#' }
#' 
#' @importFrom utils head
#' @export
fetch_pubmed_data <-
  function (pubmed_id_list,
            retstart = 0,
            retmax = 500,
            format = "xml", 
            encoding = "UTF8") 
  {
    myIDlist <- pubmed_id_list
    if ((!is.list(myIDlist)) | is.na(myIDlist$WebEnv) | is.na(myIDlist$QueryKey) | 
        is.na(myIDlist$Count) | !is.integer(as.integer(retstart)) | 
        !is.integer(as.integer(retmax))) {
      message("There is an issue with the PubMed ID list you supplied. Please, call the function again and supply the result of a <get_pubmed_ids()> call as argument. Thank you.")
      return(NULL)
    } else {
      myWebEnv <- myIDlist$WebEnv
      myKey <- myIDlist$QueryKey
      
      # ~~~ Add internal check ~~~
      if (is.na(as.numeric(myIDlist$Count))){
        return(NULL)
      } else if (as.numeric(myIDlist$Count) < 1) {
        message('There are 0 results to fetch!')
        return(NULL)
      }  
      # ~~~ end ~~~
      
      myCount <- as.numeric(as.character(myIDlist$Count))
      myRetstart = as.integer(retstart)
      if (myRetstart < 0) {
        myRetstart = 0
      }
      myRetmax <- as.integer(retmax)
      if (myRetmax > 5000) {
        myRetmax = 5000
      }
      if (myRetmax < 1) {
        myRetmax = 1
      }
      if (format[1] %in% c("medline","uilist","abstract","asn.1", "xml")) {
        myFormat <- format[1]
      } else {
        myFormat <- "xml"
      }
      typeMode <- switch(myFormat, 
                         "asn.1" = c("null", "asn.1"),
                         "xml" = c("null", "xml"),
                         "medline" = c("medline", "text"),
                         "uilist" = c("uilist", "text"),
                         "abstract" = c("abstract", "text"))
      efetch_url = paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", 
                         "db=pubmed&WebEnv=", myWebEnv, "&query_key=", myKey, 
                         "&retstart=", myRetstart, "&retmax=", myRetmax, 
                         "&rettype=", typeMode[1],"&retmode=", typeMode[2], 
                         sep = "")
      
      # api_key retrieval
      api_key <- pubmed_id_list$APIkey
      if (!is.null(api_key)) {
        efetch_url <- paste(efetch_url, "&api_key=", api_key, sep = "")
      }
      
      # initialize output
      out.data <- NULL
      
      # initialize extra params
      try_num <- 1
      t_0 <- Sys.time()
      
      # Try to fetch results
      while(is.null(out.data)) {
        
        # Timing check: kill at 2 min
        if (try_num > 1)
          Sys.sleep(time = 2)
        
        t_1 <- Sys.time()
        
        if(as.numeric(difftime(t_1, t_0, units = "mins")) > 5){
          message("Killing the request! Something is not working. Please, try again later")
          return(NULL)
        }
        
        # ENTREZ server connect
        out.data <- tryCatch({    
          tmpConnect <- suppressWarnings(url(efetch_url, open = "rb", encoding = "UTF8"))
          suppressWarnings(readLines(tmpConnect, warn = FALSE, encoding = "UTF8"))
        }, error = function(e) {
          # message(e)
          NULL
        }, finally = {
          try(suppressWarnings(close(tmpConnect)), silent = TRUE)
        })  
        
        # Check if error
        if (!is.null(out.data) && 
            class(out.data) == "character" &&
            grepl("<ERROR>", substr(paste(utils::head(out.data, n = 100), collapse = ""), 1, 250))) {
          # message(out.data)
          out.data <- NULL
        }
        try_num <- try_num + 1
      }
      
      if (is.null(out.data)) {
        message("Killing the request! Something is not working. Please, try again later")
        return(NULL)
      }
      
      if (encoding != "UTF8")
        out.data <- base::iconv(out.data, from = "UTF8", to = encoding, sub = ".")
      
      if (format[1] == "xml") {
        out.data <- paste(out.data, collapse = "")
      }
      
      return(out.data)
    }
  }


#' @title Simple PubMed Record Search by Full-length Title
#'
#' @description Query PubMed (Entrez) in a simple way via the PubMed API eSearch function. 
#' This function is designed to query PubMed using a full-length publication title as query string. 
#' It performs stopword removal from the query string before querying the PubMed server. 
#' Calling this function results in posting the results on the PubMed History Server. 
#' This allows later access to the resulting data via the fetch_pubmed_data() function, 
#' or other easyPubMed functions.
#' 
#' @usage get_pubmed_ids_by_fulltitle(fulltitle, field = "[Title]", api_key = NULL)
#' 
#' @param fulltitle String (character vector of length 1) that corresponds to the full-length 
#' publication title used for querying PubMed (titles should be used as is, without 
#' adding extra filters/tags).
#' @param field String (character vector of length 1) with a tag indicating the PubMed 
#' record field where the full-length string (fulltitle) should be searched in. By default, 
#' this points to the 'Title' field. This field can be changed (use fields supported by PubMed) 
#' as required by the user (for example, to attempt an exact-match query using a specific sentence 
#' included in the abstract of a record).
#' @param api_key String (character vector of length 1): user-specific API key to increase 
#' the limit of queries per second. You can obtain your key from NCBI.
#' 
#' @details 
#' This function will use the String provided as argument for querying PubMed via the eSearch 
#' function of the PubMed API. The Query Term should include a full-length publication title, 
#' without other PubMed operators (AND, OR, NOT) nor tags (i.e., [AU], [PDAT], 
#' [Affiliation], and so on). ESearch will post the UIDs resulting from the search operation 
#' onto the History server so that they can be used directly in a subsequent fetchPubmedData() call.
#' 
#' @return 
#' The function returns a list. The list includes the number of records found on PubMed and the first 
#' 20 PubMed IDs (UID) retrieved by the query. The list also includes QueryKey and WebEnv that are 
#' required for a subsequent fetch_pubmed_data() call.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' \dontrun{
#' ##  Search for a scientific article matching a full-length title
#' my_query <- "Body mass index and cancer risk among Chinese patients with type 2 diabetes mellitus"
#' my_field <- "[Title]"
#' # Full-length title query (designed to query titles)
#' res0 <- get_pubmed_ids(my_query)
#' print(as.numeric(res0$Count))
#' # Weird count!
#' res <- get_pubmed_ids_by_fulltitle(my_query, field = my_field)
#' # Num results = 1 as expected
#' print(as.numeric(res$Count))
#' 
#' }
#' 
#' @export
get_pubmed_ids_by_fulltitle <- 
  function(fulltitle, field = "[Title]", api_key = NULL) 
  {
    out <- get_pubmed_ids(paste("\"", fulltitle, "\"", field, sep = ""), api_key = api_key)
    if (as.numeric(out$Count) > 0) {
      return (out)
    }
    
    stopwords <- easyPubMed::PubMed_stopwords
    keys <- strsplit(fulltitle, split = "[[:space:]]")[[1]]
    keys <- tolower(keys)
    keys <- keys[!keys %in% stopwords]
    Sys.sleep(0.34)
    new_query <- paste(keys, field, sep = "", collapse = " AND ")
    return(get_pubmed_ids(new_query, api_key = api_key))
  }


#' @title Simple PubMed Record Search
#'
#' @description Query PubMed (Entrez) in a simple way via the PubMed API eSearch function. 
#' Calling this function results in posting the query results on the PubMed History Server. 
#' This allows later access to the resulting data via the fetch_pubmed_data() function, 
#' or other easyPubMed functions.
#' 
#' @usage get_pubmed_ids(pubmed_query_string, api_key = NULL)
#' 
#' @param pubmed_query_string is a string (character vector of length 1) that is used 
#' for querying PubMed (standard PubMed synthax, see reference for details).
#' @param api_key String (character vector of length 1): user-specific API key to 
#' increase the limit of queries per second. You can obtain your key from NCBI.
#' 
#' @details 
#' This function will use the String provided as argument for querying PubMed via the eSearch 
#' function of the PubMed API. The Query Term can include one or multiple words, as well as the standard 
#' PubMed operators (AND, OR, NOT) and tags (i.e., [AU], [PDAT], [Affiliation], and so on). ESearch will 
#' post the UIDs resulting from the search operation onto the History server so that they can be used directly 
#' in a subsequent fetchPubmedData() call.
#' 
#' @return 
#' The function returns a list. The list includes the number of records found on PubMed and 
#' the first 20 PubMed IDs (UID) retrieved by the query. The list also includes QueryKey and WebEnv 
#' that are required for a subsequent fetch_pubmed_data() call.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' \url{https://www.ncbi.nlm.nih.gov/books/NBK3827/#_pubmedhelp_Search_Field_Descriptions_and_}
#'
#' @examples 
#' try({
#'   ##  Search for scientific articles written by Damiano Fantini
#'   ##  and print the number of retrieved records to screen.
#'   ##  Also print the retrieved UIDs to screen.
#'   ##
#'   dami_on_pubmed <- get_pubmed_ids("Damiano Fantini[AU]")
#'   print(dami_on_pubmed$Count)
#'   print(unlist(dami_on_pubmed$IdList))
#' }, silent = TRUE)
#' 
#' @export
get_pubmed_ids <- function (pubmed_query_string, 
                            api_key = NULL) 
{
  # Silence warnings
  old_warn <- options()$warn
  options(warn = -1)
  
  # Timing
  t_0 <- Sys.time()
  
  myQuery <- as.character(pubmed_query_string)
  myQuery <- gsub(" ", "+", myQuery, fixed = TRUE)
  myPubmedURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
                       "db=pubmed&term=", myQuery, "&usehistory=y", sep = "")
  if (!is.null(api_key)) {
    myPubmedURL <- paste(myPubmedURL, "&api_key=", api_key, sep = "")
  }
  
  idXML <- NULL
  try_num <- 1
  while(is.null(idXML)) {
    
    # Timing check: kill at 2 min
    if (try_num > 1)
      Sys.sleep(time = 2)
    
    t_1 <- Sys.time()
    
    if(as.numeric(difftime(t_1, t_0, units = "mins")) > 2){
      message("Killing the request! Something is not working. Please, try again later")
      return()
    }
    
    # ENTREZ server connect
    idXML <- tryCatch({    
      IDconnect <- suppressWarnings(url(myPubmedURL, open = "rb", encoding = "UTF8"))
      idXML <- suppressWarnings(readLines(IDconnect, warn = FALSE, encoding = "UTF8"))
      idXML <- paste(idXML, collapse = "")
      if (grepl("<ERROR>", substr(idXML, 1, 250))) {
        # message(idXML)
        NULL
      } else {
        idXML
      }
    }, error = function(e) {
      # message(e)
      NULL
    }, finally = {
      try(suppressWarnings(close(IDconnect)), silent = TRUE)
    })  
    
    # Data processing (if result not null)
    myIDlist <- NULL
    
    if (!is.null(idXML)) {
      tryCatch({
        
        # Initialize collector
        myIDlist <- list()
        
        my_tags <- c("Count", "RetMax", "RetStart", 
                     "QueryKey", "WebEnv", "IdList", 
                     "TranslationSet", "QueryTranslation")
        
        # First pass
        for (j in 1:length(my_tags)) {
          ttag <- my_tags[j]
          xx <- custom_grep(idXML, tag = ttag, "char")
          myIDlist[[ttag]] <- xx[1]
        }
        
        # Try to expand IdList
        nutag <- "Id"
        xx <- myIDlist[["IdList"]]
        xx <- custom_grep(xx, "Id", format = "list")
        names(xx) <- rep("Id", length(xx))
        myIDlist[["IdList"]] <- xx
        
        # Try to expand TranslationSet
        xx <- myIDlist[["TranslationSet"]]
        myIDlist[["TranslationSet"]] <- list()
        nutag <- c("From", "To")
        for (z in nutag) {
          yy <- custom_grep(xx, z, format = "char")
          myIDlist[["TranslationSet"]][[z]] <- yy[1]
        }
        
      }, error = function(e) {
        idXML <- NULL
      })
    }
    
    # Final check!
    if(!is.list(myIDlist)) {
      idXML <- NULL  
    }
    
    try_num <- try_num + 1
  }
  
  # Wrap up and return
  myIDlist[['OriginalQuery']] <- myQuery
  myIDlist[['APIkey']] <- api_key
  
  # Restore warnings
  options(warn = old_warn)
  
  return(myIDlist)
}


#' @title Extract Publication and Affiliation Data from PubMed Records
#'
#' @description Extract Publication Info from PubMed records and cast data into a 
#' data.frame where each row corresponds to a different author. It is possible to limit
#' data extraction to first authors or last authors only, or get information about 
#' all authors of each PubMed record.
#' 
#' @usage table_articles_byAuth(pubmed_data, 
#'                              included_authors = "all", 
#'                              max_chars = 500, 
#'                              autofill = TRUE, 
#'                              dest_file = NULL, 
#'                              getKeywords = TRUE, 
#'                              encoding = "UTF8")
#' 
#' @param pubmed_data PubMed Data in XML format: typically, an XML file resulting from a 
#' batch_pubmed_download() call or an XML object, result of a fetch_pubmed_data() call.
#' @param included_authors Character: c("first", "last", "all"). Only includes information 
#' from the first, the last or all authors of a PubMed record.
#' @param max_chars Numeric: maximum number of chars to extract from the AbstractText field.
#' @param autofill Logical. If TRUE, missing affiliations are imputed according to the available 
#' values (from the same article).
#' @param dest_file String (character of length 1). Name of the file that will be written for 
#' storing the output. If NULL, no file will be saved.
#' @param getKeywords Logical. If TRUE, the operation will attempt to extract PubMed record 
#' keywords (MESH topics, keywords).
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the function 
#' base::iconv(). See iconv() help page for how to find out more about encodings that can be 
#' used on your platform. Here, we recommend using "UTF-8".
#' 
#' @details 
#' Retrieve publication and author information from PubMed data, and cast them as a data.frame.
#' 
#' @return 
#' Data frame including the following fields: c("article.title","article.abstract", "date.year", 
#' "date.month", "date.day", "journal.abbrv", "journal.title", "keywords", "auth.last", 
#' "auth.fore", "auth.address", "auth.email").
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' \dontrun{
#' ## Cast PubMed record info into a data.frame
#'
## Auto-fill enabled
#' dami_query <- "Damiano Fantini[AU]"
#' dami_on_pubmed <- get_pubmed_ids(dami_query)
#' dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed, encoding = "ASCII")
#' xx <- table_articles_byAuth(pubmed_data = dami_abstracts_xml, 
#'                             included_authors = "first", 
#'                             max_chars = 100, 
#'                             autofill = TRUE)
#' 
#' print(xx[1:5, c("pmid", "lastname", "jabbrv")])
#' #
#' ## Download records first
#' ## Also, auto-fill disabled
#' dami_query <- "Damiano Fantini[AU]"
#' curr.file <- batch_pubmed_download(dami_query, dest_file_prefix = "test_bpd_", encoding = "ASCII")
#' xx <- table_articles_byAuth(pubmed_data = curr.file[1], 
#'                             included_authors = "all", 
#'                             max_chars = 20, 
#'                             autofill = FALSE)
#' print(xx[1:5, c("pmid", "lastname", "jabbrv")])
#' 
#' }
#' @importFrom utils write.table
#' 
#' @export
table_articles_byAuth <-
  function (pubmed_data, 
            included_authors = "all", 
            max_chars = 500, 
            autofill = TRUE, 
            dest_file = NULL, 
            getKeywords = TRUE, 
            encoding = "UTF8") 
  {
    if (!included_authors %in% c("all", "first", "last")) 
      stop("Method is not supported!")
    message("Processing PubMed data ", appendLF = FALSE)
    
    #
    if(is.character(pubmed_data) && length(pubmed_data) == 1) {
      
      paper.data <- articles_to_list(pubmed_data = pubmed_data, encoding = encoding, simplify = TRUE)
      
    } else if (is.list(pubmed_data) && 
               sum(sapply(pubmed_data, function(x) {
                 grepl("PubmedArticle", substr(x, 1, 20))
               })) == length(pubmed_data)) {
      
      paper.data <- do.call(c, pubmed_data)
      
    } else {
      stop("Bad input")
    }
    
    expFields <- c("pmid", "doi", "title", "abstract", "year", "month", "day", "jabbrv", 
                   "journal", "keywords", "mesh", "lastname", "firstname", "address", "email")
    
    curWarn <- options()$warn
    options(warn = -1)
    
    papers.authors.list <- lapply(1:length(paper.data), (function(i) {
      if (length(paper.data) > 50) {
        rep.dot <- as.integer(seq(1, length(paper.data), 
                                  length.out = 50))
        if (i %in% rep.dot) 
          message(".", appendLF = FALSE)
      } else {
        message(".", appendLF = FALSE)
      }
      art <- paper.data[[i]]
      out <- tryCatch({suppressMessages(
        article_to_df(pubmedArticle = art,
                      autofill = autofill, 
                      max_chars = max_chars, 
                      getKeywords = getKeywords, 
                      getAuthors = TRUE))},
        error = function(e) { NULL })
      
      if (regexpr("MeshHeadingList", text = art) > 0) {
        out.mesh <- tryCatch({easyPubMed:::custom_grep(xml_data = art, tag = "MeshHeadingList", format = "character")}, 
                             error = function(e) NA)  
        out.mesh <- tryCatch({easyPubMed:::custom_grep(xml_data = out.mesh[1], tag = "MeshHeading", format = "list")}, 
                             error = function(e) NA)
        out.mesh <- tryCatch({lapply(out.mesh, FUN = function(XX) {
          sub('^.{10,20}UI=\\"([[:alnum:]]+)\\".+$', "\\1", x = XX)})}, 
          error = function(e) NA)
        out.mesh <- tryCatch({paste(do.call(what = c, args = out.mesh), collapse = "; ")}, 
                             error = function(e) NA)
        
      } else {
        out.mesh <- NA
      }
      
      if (is.null(out)) {
        out <- data.frame(pmid = NA, doi = NA, title = NA, 
                          abstract = NA, year = NA, month = NA, day = NA, 
                          jabbrv = NA, journal = NA, keywords = NA, mesh = NA,
                          lastname = NA, firstname = NA, 
                          address = NA, email = NA)
      } else {
        out <- data.frame(pmid = out$pmid, 
                          doi = out$doi, 
                          title = out$title, 
                          abstract = out$abstract, 
                          year = out$year, 
                          month = out$month, 
                          day = out$day, 
                          jabbrv = out$jabbrv, 
                          journal = out$journal, 
                          keywords = out$keywords, 
                          mesh = out.mesh[1],
                          lastname = out$lastname, 
                          firstname = out$firstname, 
                          address = out$address, 
                          email = out$email, 
                          stringsAsFactors = FALSE)
      }
      
      if (included_authors == "first") {
        out <- out[1, ]
      } else if (included_authors == "last") {
        out <- out[nrow(out), ]
      } 
      
      # Handle missing fields exception
      out2 <- data.frame(rebuild = (1:nrow(out))) 
      for (jj in 1:length(expFields)) {
        if (expFields[jj] %in% names(out)) {
          out2[,expFields[jj]] <- out[,expFields[jj]]
        } else {
          out2[,expFields[jj]] <- NA
        }
      }
      out2[,-1]
    }))
    message(" done!")
    
    #y <- names(papers.authors.list[[1]])
    #kp <- sapply(papers.authors.list, function(x) {sum(! y %in% names(x)) == 0 })
    #class(papers.authors.list[!kp])
    #names(papers.authors.list[!kp][[1]])
    #sum(!kp)
    
    papers.authors.df <- do.call(rbind, papers.authors.list)
    keep.rw <- apply(papers.authors.df, 1, (function(rw) {
      sum(is.na(rw)) < length(rw)
    }))
    papers.authors.df <- papers.authors.df[keep.rw, ]
    if (!is.null(dest_file)) {
      if (class(dest_file) == "character" & length(dest_file) == 1) {
        tryCatch(utils::write.table(papers.authors.df, dest_file, fileEncoding = encoding), 
                 error = function(e) {
                   NULL
                 })
      }
    }
    
    options(warn = curWarn)
    return(papers.authors.df)
  }


#' @title Trim and Format Address Information
#'
#' @description Set of rules for trimming and standardizing the format of address information 
#' retrieved from PubMed records. Affiliations including more than one address will be trimmend 
#' and only the first address will be returned.
#' 
#' @usage trim_address(addr)
#' 
#' @param addr Character string including an address as extracted from PubMed records.
#' 
#' @return 
#' Character string including a formatted and trimmed address (if available).
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples 
#' addr_string <- " 2 Dept of Urology, Feinberg School of Medicine," 
#' addr_string <- paste(addr_string, "Chicago, US; Dept of Mol Bio as well...")
#' print(addr_string)
#' print(trim_address(addr = addr_string))
#' 
#' @export
trim_address <-
  function(addr) 
  {
    out.addr <- gsub("(^([[:space:]]{0,5}[[:digit:]]{1,2}[[:space:]]{0,2}))", "", addr)
    out.addr <- gsub("(\\.).{0,5}$", "", out.addr)
    out.addr <- gsub("(\\;).*$", "", out.addr)
    return (out.addr)
  }




#' @title Fetch PubMed data based on PMIDs
#'
#' @description Retrieve PubMed Data using a list of PubMed Ids (PMIDs) as input. 
#' Results are returned as a list of PubMed articles.
#' 
#' @usage fetch_PMID_data(pmids, format = "xml", encoding = "UTF-8", verbose = TRUE)
#' 
#' @param pmids list of PubMed Ids (PMIDs) 
#' @param format Character: element specifying the output format. The following values are allowed: 
#' c("asn.1", "xml", "medline", "uilist", "abstract").
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the function base::iconv(). 
#' See iconv() help page for how to find out more about encodings that can be used on your platform. 
#' Here, we recommend using "UTF-8".
#' @param verbose logical, shall info about the processing status be printed to console 
#' 
#' @return 
#' List of PubMed articles .
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' 
#' @export
fetch_PMID_data <- function(pmids, format = "xml", encoding = "UTF-8", verbose = TRUE) 
{
  
  # Hardcoded
  pmids.per.batch = 100
  
  # custom f(x)
  query_and_fetch <- function(x, format = format, encoding = encoding) {
    q0 <- paste(paste0(x, "[PMID]"), collapse = " OR ")
    
    q1 <- easyPubMed::get_pubmed_ids(pubmed_query_string = q0)
    r1 <- easyPubMed::fetch_pubmed_data(pubmed_id_list = q1, format = format, encoding = encoding)
    r2 <- easyPubMed::articles_to_list(pubmed_data = r1, encoding = encoding, simplify = FALSE)
    
    nm2 <- lapply(r2, function(xx) { 
      y <- xx[1]; y <- gsub("</PMID>.+$", "</PMID>", y); custom_grep(y, "PMID", format = "list")[[1]]
    })
    
    
    nm2 <- as.character(do.call(c, nm2))
    
    # Order
    nuORD <- as.numeric(sapply(x, function(z) {
      which(nm2 == z)
    }))      
    r3 <- r2[nuORD]
    return(r3)
  }
  
  # standardize input
  if(is.character(pmids)) {
    pmids <- pmids[!is.na(pmids)]
  } else if (is.list(pmids)) {
    pmids <- do.call(c, pmids)
    pmids <- as.character(pmids)
    pmids <- pmids[!is.na(pmids)]
  }
  
  # split ID in small batches
  out <- list()
  ti <- Sys.time() - 2
  
  if (length(pmids) < pmids.per.batch) {
    
    out <- list(query_and_fetch(pmids, format = format, encoding = encoding))
    
  } else {
    all_i <- seq(1, (length(pmids) - 1), by = pmids.per.batch)
    all_i <- c(all_i, (length(pmids) + 1))
    
    for(j in 1:(length(all_i) - 1)) {
      
      if (verbose)
        message(".", appendLF = FALSE)
      
      i0 <- all_i[j]
      i1 <- (all_i[(j+1)] - 1)
      tdf <- as.numeric(difftime(time1 = Sys.time(), time2 =  ti, units = "sec"))
      
      # check time
      if (tdf < 1) {
        Sys.sleep(time = (1 - tdf))
      } 
      
      # generate query string
      tmp <- pmids[i0:i1]
      tmp2 <- query_and_fetch(x = tmp, format = format, encoding = encoding)
      names(tmp2) <- NULL
      
      ti <- Sys.time()
      out[[length(out) + 1]] <- tmp2
    }
  }
  
  # loop over and recompose
  OUT <- list()
  for(l1 in 1:length(out)) {
    TMP <- out[[l1]]
    for(l2 in 1:length(TMP)) {
      OUT[[length(OUT) + 1]] <- TMP[[l2]]
    }
  }
  
  if(verbose) {
    message("", appendLF = TRUE)
    message("Done!", appendLF = TRUE)
  }
  
  return(OUT)
}               


#' @title Extract Article Ids
#'
#' @description Extract PubMed Article Ids (PMID, DOI, PMCID) from a list of PubMed Article Data
#' 
#' @usage extract_article_ids(pubmed_data_list)
#' 
#' @param pubmed_data_list List of PubMed records.
#' 
#' @return 
#' Data Frame including PMID, DOI and PMCID identifiers for each PubMed data included in the input PubMed data list.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' 
#' @export
extract_article_ids <- function(pubmed_data_list) {
  
  x <- pubmed_data_list
  
  # Tmp f(x)
  my_grep <- function(X, idtype = "pubmed") {
    
    myPAT <- paste0("^.*", idtype)
    out <- tryCatch({
      if (!grepl(myPAT, X)) {
        Y <- NA
      } else {
        Y <- sub(myPAT, "", X)
        Y <- sub("</ArticleId>.*$", "", Y)
        Y <- sub("^.*>", "", Y)
      }
      Y
    }, error = function(e) NA)
    
    return(out)  
  }
  
  # dddd
  out <- lapply(x, function(xx) {
    y <- sub("</ArticleIdList>.+", "</ArticleIdList>", xx)
    y <- sub("^.+(<ArticleIdList.+$)", "\\1", y)
    
    data.frame(
      PMID = my_grep(X = y, idtype = "pubmed"),
      DOI = my_grep(X = y, idtype = "doi"), 
      PMCID = my_grep(X = y, idtype = "pmc"), 
      stringsAsFactors = FALSE)
  })
  
  # Return
  out <- do.call(rbind, out)
  return(out)
}


#' @title Retrieve PubMed Data in XML or TXT Format based on a list of PMIDs
#'
#' @description Retrieve PubMed records from Entrez starting from a list of PumBed IDs (PMIDs). 
#' Data are downloaded in the XML or TXT format and are retrieved in batches of up to 55 PMIDs
#' at a time.
#' 
#' @usage fetch_pubmed_data_by_PMID(pmids, 
#'                                  batch = 50, 
#'                                  format = "xml", 
#'                                  encoding = "UTF-8", 
#'                                  verbose = TRUE)
#' 
#' @param pmids file name or PMID character vector; this can point to a file including a list of PMIDs 
#' (for example, a PMID file downloaded from the PubMed website), or a character vector including
#' a list of PMIDs (example: c("31867428", "31564976", "31564903")).
#' @param batch Integer (>10 and <55): number of PMIDs to be retrieved at each iteration.
#' @param format Character: element specifying the output format. The following values are allowed: 
#' c("asn.1", "xml", "medline", "uilist", "abstract"). Note that "uilist" will only retrieve PMIDs.
#' @param encoding The encoding of an input/output connection can be specified by name 
#' (for example, "ASCII", or "UTF-8", in the same way as it would be given to the function base::iconv(). 
#' See iconv() help page for how to find out more about encodings that can be used on your platform. 
#' Here, we recommend using "UTF-8".
#' @param verbose logical, shall info about the progress of the data download be printed/messaged to console.
#' 
#' @details 
#' Retrieve PubMed records based on a list of PMIDs. 
#' Records are retrieved from Entrez via the PubMed API efetch function. 
#' Data will be downloaded on the fly (no files are saved locally).
#' 
#' @return 
#' An object (vector) of class "character". If format is set to "xml" (default), a single String including all 
#' PubMed records (with XML tags embedded) is returned. If a different format is selected, a vector of strings 
#' is returned, where each row corresponds to a line of the output document.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/}
#'
#' @examples
#' try({ 
#'   ## Example 01: retrieve data in XML format
#'   library("easyPubMed")
#'   PMIDs <- c("31506076", "31410434", "31299087", "31552170", 
#'              "31546103", "31518875", "31518741", "31510010",
#'              "31509964", "31509747", "31506606", "31503007")
#'   pmxml <- fetch_pubmed_data_by_PMID(PMIDs)
#'   # length is 1
#'   length(pmxml)
#'   # show an excerpt
#'   substr(pmxml, 1, 299)
#' }, silent = TRUE)
#' 
#' @importFrom utils head
#' @export
fetch_pubmed_data_by_PMID <- function(pmids, 
                                      batch = 50, 
                                      format = "xml", 
                                      encoding = "UTF-8", 
                                      verbose = TRUE) {
  
  if (length(pmids) == 1 && file.exists(pmids)) {
    pmids <- suppressMessages(
      suppressWarnings(
        readLines(con = pmids)))
  }
  
  stopifnot(is.vector(pmids) && length(pmids) > 0)
  stopifnot(is.numeric(batch) && length(batch) == 1)
  stopifnot(batch > 10 && batch < 55)
  
  # total n
  tot.n <- length(pmids)
  all.i <- seq(1, tot.n, by = batch)
  
  Y <- list()
  for (ij in 1:length(all.i)) {
    
    i <- all.i[ij]
    
    if (verbose){
      if (ij %% 20 == 0) {
        message(".", appendLF = TRUE)
      } else {
        message(".", appendLF = FALSE)
      }
    }
    
    i1 <- i
    i2 <- (i + batch - 1)
    if (i2 > tot.n) { i2 <- tot.n }
    tmp.ids <- pmids[(i1:i2)] 
    tmp.ids <- paste0(tmp.ids, "[PMID]")
    tmp.ids <- paste(tmp.ids, collapse = " OR ")
    
    tmp.q <- get_pubmed_ids(pubmed_query_string = tmp.ids)
    
    tmp.d <- fetch_pubmed_data(pubmed_id_list = tmp.q, 
                               retstart = 0, retmax = 1000, 
                               format = format, encoding = encoding)
    
    Y[[length(Y) + 1]] <- tmp.d
    Sys.sleep(time = 0.35)
  }
  
  if (length(Y) > 1) {
    
    if (format == "xml") {
      
      YY <- list()
      for (j in 1:length(Y)) {
        if (length(YY) == 0) {
          
          tmp.ei <- regexpr(pattern = "<PubmedArticleSet>", text = Y[[j]])
          tmp.ii <- regexpr(pattern = "<\\?xml", text = Y[[j]])
          
          tmp.top <- substr(x = Y[[j]], start = as.numeric(tmp.ii), 
                            stop = (as.numeric(tmp.ei) + as.numeric(attributes(tmp.ei)$match.length) - 1))
          
          YY[[1]] <- tmp.top
        }
        
        str.fi <- regexpr(pattern = "<PubmedArticleSet>", text = Y[[j]]) 
        str.fe <- regexpr(pattern = "<\\/PubmedArticleSet>", text = Y[[j]]) 
        
        if (as.numeric(str.fi) > 0 && as.numeric(str.fe) > 0) {
          yz <- substr(x = Y[[j]], 
                       start = as.numeric(str.fi) + as.numeric(attributes(str.fi)$match.length), 
                       stop = (as.numeric(str.fe) - 1))
          
          YY[[length(YY) + 1]] <- yz
        }
      }
      YY[[length(YY) + 1]] <- "</PubmedArticleSet>"
      
      # merge together
      yy <- do.call(paste0, YY)
      
    } else if (format == "uilist") {
      
      yy <- do.call(c, Y)
      yy <- unique(yy)
      
    } else {
      
      yy <- do.call(c, lapply(Y, function(xx) {c(xx, "")}))
    }
    
  } else {
    
    yy <- Y[[1]]
  }
  
  # return
  return(yy)
}


#' @title PubMed Records downloaded and analyzed via easyPubMed
#'
#' @description This dataset includes a collection of 4 examples showing how to download and analyze records 
#' from PubMed by using easyPubMed. Each element in the EPMsamples list corresponds to a different query 
#' and/or analysis. Also, each element of EPMsamples is a list including intermediates and notes about the analysis.
#'
#' @usage data("EPMsamples")
#' 
#' @format The dataset is formatted as a list including 4 elements: 
#' 
#' * `DF_papers_abs`: List of 4
#' 
#' * `DF_papers_std`: List of 4
#' 
#' * `NUBL_dw18`: List of 3
#' 
#' * `NUBL_1618`: List of 5
#' 
#'   
#' @details The dataset was built as described in this vignette: \url{https://www.data-pulse.com/projects/Rlibs/vignettes/building_the_easyPubMed_EPMsamples_dataset.html}
#'
#' @examples 
#' ## Display some contents
#' data("EPMsamples")
#' # The following examples are focused on example query #4 (i.e., NUBL_1618)
#' # Display Query String used for collecting the data
#' print(EPMsamples$NUBL_1618$qry_st)
#' # show one PubMed record element from the IL vector
#' NU_records <- EPMsamples$NUBL_1618$rec_lst
#' cat(NU_records[[1]])
#' # cast PM recort to data.frame
#' BL_df <- article_to_df(NU_records[[6]], max_chars = 0)
#' print(BL_df)
"EPMsamples"


#' @title PubMed Query Stopwords 
#'
#' @description Collection of 133 Stopwords that can be removed from query strings to improve the 
#'              accuracy of exact-match PubMed queries.
#'
#' @usage data("PubMed_stopwords")
#' 
#' @format A character vector including all PubMed stopwords tat are typically filtered out from queries. 
#'   
#' @details Number of stopwords included, n=133.
#'
#' @examples 
#' ## Display some contents
#' data("PubMed_stopwords")
#' head(PubMed_stopwords)
"PubMed_stopwords"

#' @title Retrieve and Process Scientific Publication Records from Pubmed
#' 
#' @description Query NCBI Entrez and retrieve PubMed records in XML or TXT format. PubMed records 
#' can be downloaded and saved as XML or text files. Data integrity is enforced during data download, 
#' allowing to retrieve and save very large number of records effortlessly. PubMed records can be processed 
#' to extract publication- and author-specific information.
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @references \url{https://www.data-pulse.com/dev_site/easypubmed/}
#'
#' @examples
#' try({
#'   ## Example 01: retrieve data in TXT format
#'   dami_query_string <- "Damiano Fantini[AU]"
#'   dami_on_pubmed <- get_pubmed_ids(dami_query_string)
#'   dami_papers <- fetch_pubmed_data(dami_on_pubmed, format = "abstract")
#'   dami_papers[dami_papers == ""] <- "\n"
#'   cat(paste(dami_papers[1:65], collapse = ""))
#'   #
#' }, silent = TRUE)
#' 
#' \dontrun{
#' ## Example 02: retrieve data in XML format
#' library("easyPubMed")
#' dami_query_string <- "Damiano Fantini[AU] AND 2018[PDAT]"
#' dami_on_pubmed <- get_pubmed_ids(dami_query_string)
#' dami_papers <- fetch_pubmed_data(dami_on_pubmed)
#' titles <- sapply(dami_papers, custom_grep, tag = "ArticleTitle", format = "char", USE.NAMES = FALSE)
#' print(titles)
#' #
#' ## Example 03: retrieve data from PubMed and save as XML file
#' ml_query <- "Machine Learning[TI] AND 2016[PD]"
#' out1 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
#' x <- paste(readLines(out1[1], n = 10), collapse = "\n")
#' cat(x)
#' #
#' ## Example 04: retrieve data from PubMed and save as TXT file
#' ml_query <- "Machine Learning[TI] AND 2016[PD]"
#' out2 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180, format = "medline")
#' x <- paste(readLines(out1[1], n = 30), collapse = "\n")
#' cat(x)
#' #
#' ## Example 05: extract information from a single PubMed record 
#' ml_query <- "Machine Learning[TI] AND 2016[PD]"
#' out3 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
#' PM_data <- articles_to_list(out3[1])
#' PM_record_df <- article_to_df(PM_data[[80]])
#' print(PM_record_df[1,])
#' print(PM_record_df[,"address"])
#' #
#' ## Example 06: query PubMed and extract information from multiple records in one step 
#' ml_query <- "Machine Learning[TI] AND 2016[PD]"
#' out4 <- batch_pubmed_download(pubmed_query_string = ml_query, batch_size = 180)
#' PM_tab <- table_articles_byAuth(out4[1], autofill = TRUE, included_authors = "last")
#' PM_tab$address <- substr(PM_tab$address, 1, 12)
#' PM_tab[50:70,c("pmid", "jabbrv", "year", "lastname", "address")]
#' }
#' 
#' @keywords internal
"_PACKAGE"


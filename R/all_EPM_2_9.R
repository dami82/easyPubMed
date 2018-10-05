articles_to_list <- function(pubmed_data) 
{
  # check if it is a XMLAbstractDocument or a file
  options(warn = -1)
  if (sum(class(pubmed_data) == "XMLAbstractDocument") > 0) {
    out <- tryCatch(XML::xpathSApply(pubmed_data, "//PubmedArticle", XML::saveXML), error = function(e) { NULL })
  } else if (class(pubmed_data) == "character") {
    out <- tryCatch(XML::xmlParse(pubmed_data), error = function(e) {NULL})    
    if (!is.null(out))
      out <- tryCatch(XML::xpathApply(out, "//PubmedArticle", XML::saveXML), error = function(e) { NULL })    
  } else {
    out <- NULL
  }
  if (is.null(out)) {
    message("An error occurred")
  }
  options(warn = 0)
  return(out)
}

article_to_df <- function(pubmedArticle, 
                          autofill = FALSE, 
                          max_chars = 500, 
                          getKeywords = FALSE,
                          getAuthors = TRUE) 
{
  #
  options(warn = -1)
  
  # initial check
  # expected cols = 14
  # "pmid", "doi", "title", "abstract", "year", "month", "day", "jabbrv", "journal", 
  # "keywords", "lastname", "firstname", "address", "email" 
  
  # Global Check!
  if (class(pubmedArticle) != "character" |
      regexpr("(<PubmedArticle)(.+)(\\/PubmedArticle>)", pubmedArticle) < 0 )
  {
    message("An error occurred")
    return(NULL)
  }
  
  # max_chars Check
  if (!is.numeric(max_chars)) {
    max_chars <- 500
  } else if (max_chars < 0) {
    max_chars <- -1  
  }
  
  # Get started
  tryCatch({
    
    tmp.article <- custom_grep(xml_data = pubmedArticle, tag = "PubmedArticle", format = "char")
    if (is.null(tmp.article)) 
    {
      message("An error occurred")
      return(NULL)
    }
    
    # Title
    tmp.title <- custom_grep(xml_data = tmp.article, tag = "ArticleTitle", format = "char")
    if (length(tmp.title) > 1){
      tmp.title <- paste(tmp.title, collapse = " ", sep = " ")
    } else if (length(tmp.title) < 1) {
      tmp.title <- NA
    }
    
    # Abstract
    tmp.abstract <- custom_grep(xml_data = tmp.article, tag = "AbstractText", format = "char")
    if (length(tmp.abstract) > 1){
      tmp.abstract <- paste(tmp.abstract, collapse = " ", sep = " ")
      if(max_chars >= 0) {
        tmp.abstract <- gsub("</{0,1}i>", "", tmp.abstract, ignore.case = T)
        tmp.abstract <- gsub("</{0,1}b>", "", tmp.abstract, ignore.case = T)
        tmp.abstract <- gsub("</{0,1}sub>", "", tmp.abstract, ignore.case = T)
        tmp.abstract <- gsub("</{0,1}exp>", "", tmp.abstract, ignore.case = T)
        
        tmp.abstract <- substr(tmp.abstract, 0, max_chars)
      }
    } else if (length(tmp.abstract) < 1) {
      tmp.abstract <- NA
    } else {
      if(max_chars >= 0) {
        tmp.abstract <- substr(tmp.abstract, 0, max_chars)
        tmp.abstract <- gsub("</{0,1}i>", "", tmp.abstract, ignore.case = T)
        tmp.abstract <- gsub("</{0,1}b>", "", tmp.abstract, ignore.case = T)
        tmp.abstract <- gsub("</{0,1}sub>", "", tmp.abstract, ignore.case = T)
        tmp.abstract <- gsub("</{0,1}exp>", "", tmp.abstract, ignore.case = T)
        
      }
    }
    
    # Dates, if any
    my.dateType <- c("DateCompleted",  "DateCreated",  "DateRevised",  "PubDate")
    sel.dateType <-which(sapply(my.dateType, (function(xi) {
      regexpr(xi, tmp.article) > 0
    })))
    if (length(sel.dateType) < 1) {
      tmp.date <- c(Year=NA, Month=NA, Day=NA)
    } else {
      sel.dateType <- sel.dateType[1]
      tmp.date <- custom_grep(xml_data = tmp.article, tag = my.dateType[sel.dateType], format = "char")
      tmp.date <- sapply(c("Year", "Month", "Day"), (function(tt){
        tdat.el <- custom_grep(xml_data = tmp.date, tag = tt, format = "char")
        ifelse(is.null(tdat.el), NA, tdat.el[1])
      }))
    }
    
    # Fetch ID string
    tmp.paperID  <- custom_grep(xml_data = tmp.article, tag = "ArticleIdList", format = "char")
    if (is.null(tmp.paperID)) 
    {
      message("An error occurred")
      return(NULL)
    } else {
      tmp.paperID <- gsub("[[:space:]]", "", tmp.paperID[1])
    }
    
    # Get PMID
    tmp.PMID <- gsub("^(.*ArticleIdIdType=\\\"pubmed\\\")([[:space:]]|[[:alnum:]]){0,20}>", "", tmp.paperID)
    tmp.PMID <- gsub("<.*$", "", tmp.PMID)
    
    # Get DOI
    tmp.DOI <- gsub("^(.*ArticleIdIdType=\\\"doi\\\")([[:space:]]|[[:alnum:]]){0,20}>", "", tmp.paperID)
    tmp.DOI <- gsub("<.*$", "", tmp.DOI)
    
    # Get Journal Abbrv
    tmp.jabbrv  <- custom_grep(xml_data = tmp.article, tag = "ISOAbbreviation", format = "char")
    tmp.jabbrv <- ifelse(is.null(tmp.jabbrv), NA, tmp.jabbrv)
    
    # Get Title
    tmp.journal <- custom_grep(xml_data = tmp.article, tag = "Title", format = "char")
    tmp.journal <- ifelse(is.null(tmp.journal), NA, tmp.journal)
    
    # Fetch Keywords ----MeshHeading
    tmp.keys <- tryCatch({
      if (getKeywords) {
        tmp.keys <- custom_grep(xml_data = tmp.article, 
                                tag = "Keyword", 
                                format = "char")
        
        tmp.mesh <- custom_grep(xml_data = tmp.article, 
                                tag = "MeshHeading", 
                                format = "char")
        
        if (length(tmp.mesh) > 0) {
          tmp.mesh <- sapply(tmp.mesh, function(xxm) {
            custom_grep(xml_data = xxm, 
                        tag = "DescriptorName", 
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
    }, error = function(e) {NA})
    
    # vector with all unique fields extracted o far
    tmp.resout <- c(pmid=tmp.PMID, 
                    doi=tmp.DOI, 
                    title=tmp.title,
                    abstract=tmp.abstract,
                    year = as.vector(tmp.date[1]),
                    month = as.vector(tmp.date[2]),
                    day = as.vector(tmp.date[3]),
                    jabbrv=tmp.jabbrv,
                    journal=tmp.journal,
                    keywords=tmp.keys)
    
    # Slow part - authors
    tmp.authors <- custom_grep(xml_data = tmp.article, tag = "AuthorList", format = "char")
    
    if (length(tmp.authors) < 1 | !getAuthors) {
      # Set every placeholder with NA
      final.mat <- data.frame(rbind(c(tmp.resout, 
                                      lastname=NA, 
                                      firstname=NA, 
                                      address=NA, 
                                      email=NA)), stringsAsFactors = FALSE)
    } else {
      author.list <- custom_grep(xml_data = tmp.authors, tag = "Author", format = "char")
      final.mat <- do.call(rbind, lapply(author.list, (function(al) {
        tmp.lastnm <- custom_grep(xml_data = al, tag = "LastName", format = "char")
        tmp.firstnm <- custom_grep(xml_data = al, tag = "ForeName", format = "char")
        tmp.email <- regexpr("([[:alnum:]]|\\.|\\-\\_){3,200}@([[:alnum:]]|\\.|\\-\\_){3,200}(\\.)([[:alnum:]]){2,6}", al)
        if (tmp.email > 0) {
          tmp.email <- substr(al, tmp.email, tmp.email + attributes(tmp.email)$match.length -1 )
        } else {
          tmp.email <- NA
        }
        #
        if (regexpr("Affiliation", al) > 0) {
          tmp.add <- custom_grep(al, "Affiliation", format = "char")[1]
          tmp.add <- trim_address(tmp.add)
        } else {
          tmp.add <- NA
        }
        c(tmp.resout, 
          lastname=tmp.lastnm, 
          firstname=tmp.firstnm, 
          address=tmp.add, 
          email=tmp.email)
        #
      })))
      rownames(final.mat) <- NULL
      
      final.mat <- data.frame(final.mat, stringsAsFactors = FALSE)
      DESELECT <- is.na(final.mat$lastname) | is.na(final.mat$firstname)
      if (length(DESELECT) > 0 & sum(DESELECT) > 0)
        final.mat <- final.mat[!DESELECT, ]
      #
      if (autofill){
        tmp.address <- final.mat[,"address"]
        na.pos <- is.na(tmp.address)
        if (sum(na.pos) != length(tmp.address)) {
          tmp.list <- lapply(tmp.address, function(x) {x} ) 
          cur.add <-  tmp.list[[(which(!na.pos)[1])]]
          for (i in 1:length(na.pos)){
            if(na.pos[i]){
              tmp.list[[i]] <- cur.add
            } else {
              cur.add <- tmp.list[[i]]
            }
          }
          final.mat[,"address"] <- do.call(c, tmp.list)
        }
      }
    }
    
    # Final check and return
    if (ncol(final.mat) != 14) {
      final.mat <- NULL
    }
  }, error = function(e) {NULL}, 
  finally = {
    options(warn = 0)
    return(final.mat)
  })
}

get_pubmed_ids <- function (pubmed_query_string, api_key = NULL) 
{
  myQuery <- as.character(pubmed_query_string)
  myQuery <- gsub(" ", "+", myQuery, fixed = TRUE)
  myPubmedURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
                       "db=pubmed&term=", myQuery, "&usehistory=y", sep = "")
  if (!is.null(api_key)) {
    myPubmedURL <- paste(myPubmedURL, "&api_key=", api_key, sep = "")
  }
  IDconnect <- url(myPubmedURL, open = "rb")
  idXML <- readLines(IDconnect, warn = FALSE, encoding = "UTF-8")
  idXML <- XML::xmlTreeParse(idXML)
  close.connection(IDconnect)
  myIDlist <- XML::xmlToList(idXML)
  return(myIDlist)
}

batch_pubmed_download <- function (pubmed_query_string, 
                                   dest_dir = NULL,
                                   dest_file_prefix = "easyPubMed_data_",
                                   format = "xml",
                                   api_key = NULL,
                                   batch_size = 400, 
                                   res_cn = 1) 
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
    diff_time <- my_rtime - (as.numeric(Sys.time() - cur_time))
    if (diff_time > 0) {
      Sys.sleep(diff_time)
    }
    cur_time <- Sys.time()
    myQuery <- tryCatch(easyPubMed::get_pubmed_ids(pubmed_query_string, api_key = api_key), 
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
          
        myQuery <- tryCatch(easyPubMed::get_pubmed_ids(pubmed_query_string), 
                            error = function(e) NULL)
        
        diff_time <- my_rtime - (as.numeric(Sys.time() - cur_time))
        if (diff_time > 0) {
          Sys.sleep(diff_time)
        }
        cur_time <- Sys.time()
        
        tmpPapers <- tryCatch(easyPubMed::fetch_pubmed_data(pubmed_id_list = myQuery, 
                                                            retstart = myRetstart, 
                                                            retmax = myRetmax,
                                                            api_key = api_key,
                                                            format = format), 
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
      if (format == "xml") {
        doSaveData <- tryCatch(XML::saveXML(tmpPapers,
                                            paste(dest_file_prefix, myExt, j, ".xml", sep = "")), 
                               error = function(e) "ERROR")
      } else {
        tmp.dest.file <- paste(dest_file_prefix, myExt, j, ".txt", sep = "")
        doSaveData <- tryCatch(write(tmpPapers, tmp.dest.file), 
                               error = function(e) "ERROR")
        if(is.null(doSaveData))
          doSaveData <- tmp.dest.file
      }
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

custom_grep <- function(xml_data, 
                        tag, 
                        format = "list")
{
  x <- xml_data
  tag.op <- paste("\\<", tag, "((\\>)|([[:space:]]([^[<]]*)\\>))", sep = "")
  tag.cl <- paste("(<\\/)", tag, "(\\>)", sep = "")
  #
  out.result <- list()
  i = 1
  while (nchar(x) > 0 &
         regexpr(tag.op, x) > 0 &
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

fetch_pubmed_data <- function (pubmed_id_list, 
                               retstart = 0, 
                               retmax = 500,
                               api_key = NULL,
                               format = "xml") 
{
  myIDlist <- pubmed_id_list
  if ((!is.list(myIDlist)) | is.na(myIDlist$WebEnv) | is.na(myIDlist$QueryKey) | 
      is.na(myIDlist$Count) | !is.integer(as.integer(retstart)) | 
      !is.integer(as.integer(retmax))) {
    stop("There is an issue with the PubMed ID list you supplied. Please, call the function again and supply the result of a <get_pubmed_ids()> call as argument. Thank you.")
  } else {
    myWebEnv <- myIDlist$WebEnv
    myKey <- myIDlist$QueryKey
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
    if (!is.null(api_key)) {
      efetch_url <- paste(efetch_url, "&api_key=", api_key, sep = "")
    }
    tmpConnect <- url(efetch_url, open = "rb")
    out.data <- readLines(tmpConnect, warn = FALSE, encoding = "UTF-8")
    if (myFormat == "xml") {
      out.data <- XML::xmlParse(out.data, encoding = "UTF-8")  
    }
    close.connection(tmpConnect)
    return(out.data)
  }
}

table_articles_byAuth <- function (pubmed_data, included_authors = "all", max_chars = 500, 
                                   autofill = TRUE, dest_file = NULL, getKeywords = TRUE) 
{
  if (!included_authors %in% c("all", "first", "last")) 
    stop("Method is not supported!")
  message("Processing PubMed data ", appendLF = FALSE)
  paper.data <- articles_to_list(pubmed_data)
  expFields <- c("pmid", "doi", "title", "abstract", "year", "month", "day", "jabbrv", 
                 "journal", "keywords", "lastname", "firstname", "address", "email")
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
    out <- tryCatch(article_to_df(pubmedArticle = art, 
                                  autofill = autofill, 
                                  max_chars = max_chars, 
                                  getKeywords = getKeywords, 
                                  getAuthors = TRUE), error = function(e) {
                                    NULL
                                  })
    if (is.null(out)) {
      out <- data.frame(pmid = NA, doi = NA, title = NA, 
                        abstract = NA, year = NA, month = NA, day = NA, 
                        jabbrv = NA, journal = NA, keywords = NA, lastname = NA, firstname = NA, 
                        address = NA, email = NA)
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
    if (class(dest_file) == "character" & length(dest_file) == 
        1) {
      tryCatch(utils::write.table(papers.authors.df, dest_file), 
               error = function(e) {
                 NULL
               })
    }
  }
  return(papers.authors.df)
}

trim_address <- function(addr) 
{
  out.addr <- gsub("(^([[:space:]]{0,5}[[:digit:]]{1,2}[[:space:]]{0,2}))", "", addr)
  out.addr <- gsub("(\\.).{0,5}$", "", out.addr)
  out.addr <- gsub("(\\;).*$", "", out.addr)
  return (out.addr)
}


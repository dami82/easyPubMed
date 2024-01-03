
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



## ---------------
## --- Classes ---
## ---------------

#' Class easyPubMed.
#'
#' Class easyPubMed defines objects that represent PubMed Query jobs and
#' the corresponding results. Briefly, these objects are initialized using
#' information that will guide the communication with the NCBI Entrez server. 
#' Also, easyPubMed objects are used to store raw and processed data retrieved 
#' from Pubmed. 
#' 
#' 
#' @slot query String (character vector of length 1) corresponding to 
#' the PubMed request submitted by the user.
#'
#' @slot meta List including meta information about the PubMed Query job.
#' 
#' @slot uilist List including all unique identifiers corresponding to 
#' the Pubmed records returned by the query. Can be empty.
#'
#' @slot raw List including the raw data (in `xml` or `medline` format) 
#' retrieved from the NCBI eFetch server. Can be empty.
#' 
#' @slot data Data.frame including processed data based on 
#' the xml raw data retrieved from PubMed. 
#' 
#' @slot misc List including additional information.
#' 
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#' 
#' @name easyPubMed-class
#' @rdname easyPubMed-class
#' @exportClass easyPubMed
#' @export
setClass("easyPubMed",
         slots = list(
           query="character",
           meta="list",
           uilist = "list",
           raw="list",
           data="data.frame",
           misc="list"))



## ---------------------
## --- Basic Methods ---
## ---------------------



#' Constructor method of the easyPubMed Class.
#'
#' @rdname easyPubMed-class
#' @importFrom methods callNextMethod
#' @importFrom utils installed.packages
#' 
#' @param .Object The easyPubMed object being built.
#' @param query_string String (character vector of length 1) corresponding 
#' to the user-provided text of the query to be submitted to PubMed. 
#' @param job_info List, this should be the output of `EPM_job_split()`.
#' 
#' 
#' @aliases initialize,easyPubMed-method
setMethod("initialize", "easyPubMed",
          function(.Object, query_string, job_info) {
            .Object <- callNextMethod(.Object)
            
            # Check args
            stopifnot(is.character(query_string), length(query_string) == 1, 
                      !is.na(query_string),  nchar(query_string) > 0, 
                      is.list(job_info), 
                      sum(c("meta", "query_guide") %in% names(job_info)) == 2, 
                      is.list(job_info$meta), 
                      is.data.frame(job_info$query_guide), 
                      nrow(job_info$query_guide) > 0)
            
            # Collect params
            epmver <- tryCatch({installed.packages()['easyPubMed','Version']}, 
                               error = function(e) { NA })
            
            mymeta <- job_info$meta[-1]
            mymeta$query_date <- format(Sys.time(), tz = 'GMT')
            
            # Raw data params
            mymeta$raw_format <- NA
            mymeta$raw_encoding <- NA
            mymeta$raw_data_embedded <- NA
            mymeta$raw_record_num <- NA
            mymeta$raw_date <- NA
            
            # Processed data params
            mymeta$processed_compact_output <- NA
            mymeta$processed_record_num <- NA
            mymeta$processed_ref_type <- NA
            
            # Unique Key
            mymeta$UID <- EPM_init_unique_key()
            
            # EPM version
            mymeta$EPM_version <- epmver
            
            # Initialize the object
            .Object@query <- query_string
            .Object@meta <- mymeta
            .Object@uilist <- list()
            .Object@raw <- list()
            .Object@data <- data.frame()
            .Object@misc <- list(job_list = job_info$query_guide)
            .Object
          })



## ---------------------
## --- Basic Methods ---
## ---------------------


#' Show method of the easyPubMed Class.
#'
#' @rdname easyPubMed-show
#' 
#' @param object the `easyPubMed` object being shown.
#' 
#' @aliases show,easyPubMed-method
setMethod("show", signature(object = "easyPubMed"),
          function(object) {
            
            # get query string
            str0 <- tryCatch(object@query, 
                             error = function(e) { NA })
            if(!is.na(str0) && nchar(str0) < 1) { str0 <- NA }
            
            # get expected counts
            cnt0 <- tryCatch(object@meta$exp_count, 
                             error = function(e) { NA })
            
            # get number of uilist elems
            uin0 <- tryCatch(length(object@uilist), 
                             error = function(e) { NA })
            
            # expected number of sub-jobs
            sjn0 <- tryCatch(nrow(object@misc$job_list), 
                             error = function(e) { NA })
            
            # fetched data
            fdn0 <- tryCatch(length(object@raw), 
                             error = function(e) { NA })
            
            # fetched data format
            ffr0 <- tryCatch(object@meta$raw_format, 
                             error = function(e) { NA })
            if (is.na(ffr0)) { ffr0 <- '' }
            
            # proc data
            prc0 <- tryCatch(nrow(object@data), 
                             error = function(e) { NA })
            prc1 <- tryCatch(as.integer(object@meta$processed_record_num), 
                             error = function(e) { NA })
            if (is.na(prc1)) { prc1 <- 0 }
            
            # EPM ver
            evr0 <- tryCatch(object@meta$EPM_version, 
                             error = function(e) { NA })
            
            #
            # print info to console
            cat(       '   --- easyPubMed object --- \n \n')
            cat(paste0('    query_string: <', str0, '> \n'))
            cat(paste0('  expected_count: n=', cnt0, ' \n'))
            cat(paste0('      PMID_count: n=', uin0, ' \n'))
            cat(paste0('   fetch_subjobs: n=', sjn0, ' \n'))
            cat(paste0(' fetched_records: n=', fdn0, ' \n'))
            cat(paste0('   record_format: ', ffr0, ' \n'))
            
            cat(paste0('  processed_data: n=', prc1, ' (records)\n'))
            cat(paste0('  processed_data: n=', prc0, ' (rows)\n'))
            
            cat(paste0('  easyPubMed_ver: ', evr0, ' \n\n'))
            
          })



#' Print method of the easyPubMed Class.
#'
#' @rdname easyPubMed-print
#' 
#' @param x the `easyPubMed` object being shown.
#' @importFrom methods show
#' 
#' @aliases print,easyPubMed-method
setMethod("print", signature(x = "easyPubMed"),
          function(x) {
            show(x)
          })



## ---------------
## --- Getters ---
## ---------------




#' Method getEPMData.
#' 
#' Retrieve processed data from an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`.
#' 
#' @rdname getEPMData
#' @exportMethod getEPMData  
setGeneric("getEPMData", function(x) {
  standardGeneric("getEPMData")
})

#' @rdname getEPMData
#' @aliases getEPMData,easyPubMed-method
setMethod("getEPMData", signature(x="easyPubMed"),
          function(x) {
            
            out <- x@data
            out
          })



#' Method getEPMQuery.
#' 
#' Retrieve the user-provided query string from an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' 
#' @rdname getEPMQuery
#' @exportMethod getEPMQuery  
setGeneric("getEPMQuery", function(x) {
  standardGeneric("getEPMQuery")
})

#' @rdname getEPMQuery
#' @aliases getEPMQuery,easyPubMed-method
setMethod("getEPMQuery", signature(x="easyPubMed"),
          function(x) {
            
            out <- x@query
            out
          })



#' Method getEPMMeta.
#' 
#' Retrieve meta data from an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' 
#' @rdname getEPMMeta
#' @exportMethod getEPMMeta  
setGeneric("getEPMMeta", function(x) {
  standardGeneric("getEPMMeta")
})

#' @rdname getEPMMeta
#' @aliases getEPMMeta,easyPubMed-method
setMethod("getEPMMeta", signature(x="easyPubMed"),
          function(x) {
            
            out <- x@meta
            out
          })



#' Method getEPMUilist.
#' 
#' Retrieve the list of unique record identifiers (PMIDs) 
#' from an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' 
#' @rdname getEPMUilist
#' @exportMethod getEPMUilist  
setGeneric("getEPMUilist", function(x) {
  standardGeneric("getEPMUilist")
})

#' @rdname getEPMUilist
#' @aliases getEPMUilist,easyPubMed-method
setMethod("getEPMUilist", signature(x="easyPubMed"),
          function(x) {
            
            out <- x@uilist
            tryCatch({
              out <- do.call(c, out)
              names(out) <- NULL
            }, error = function(e) { NULL })
            out
          })



#' Method getEPMRaw.
#' 
#' Retrieve the raw PubMed record data stored in an `easyPubMed` object. 
#' 
#' @param x an object of class `easyPubMed`. 
#' 
#' @rdname getEPMRaw
#' @exportMethod getEPMRaw  
setGeneric("getEPMRaw", function(x) {
  standardGeneric("getEPMRaw")
})

#' @rdname getEPMRaw
#' @aliases getEPMRaw,easyPubMed-method
setMethod("getEPMRaw", signature(x="easyPubMed"),
          function(x) {
            
            out <- x@raw
            out
          })



#' Method getEPMMisc.
#' 
#' Retrieve miscellaneous information stored in an `easyPubMed` object. 
#' 
#' @param x an object of class `easyPubMed`. 
#' 
#' @rdname getEPMMisc
#' @exportMethod getEPMMisc  
setGeneric("getEPMMisc", function(x) {
  standardGeneric("getEPMMisc")
})

#' @rdname getEPMMisc
#' @aliases getEPMMisc,easyPubMed-method
setMethod("getEPMMisc", signature(x="easyPubMed"),
          function(x) {
            
            out <- x@misc
            out
          })



#' Method getEPMJobList.
#' 
#' Retrieve the list of record retrieval sub-jobs 
#' from an `easyPubMed` object. 
#' Record retrieval sub-jobs are stored in a `data.frame` and each
#' row corresponds to an independent non-overlapping 
#' PubMed query. This `data.frame` guides the record retrieval process. 
#' The `data.frame` is obtained from the 
#' 'misc' slot of an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' 
#' @rdname getEPMJobList
#' @exportMethod getEPMJobList  
setGeneric("getEPMJobList", function(x) {
  standardGeneric("getEPMJobList")
})

#' @rdname getEPMJobList
#' @aliases getEPMJobList,easyPubMed-method
setMethod("getEPMJobList", signature(x="easyPubMed"),
          function(x) {
            
            out <- tryCatch({x@misc$job_list}, 
                            error = function(e) { NA })
            out
          })



## ---------------
## --- Setters ---
## ---------------


#' Method setEPMData.
#' 
#' Attach (or replace) processed data to an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`.
#' @param y `data.frame` including processed data.
#' 
#' @rdname setEPMData
#' @exportMethod setEPMData  
setGeneric("setEPMData", function(x, y) {
  standardGeneric("setEPMData")
})

#' @rdname setEPMData
#' @aliases setEPMData,easyPubMed,data.frame-method
setMethod("setEPMData", signature(x="easyPubMed", y="data.frame"),
          function(x, y) {
            
            x@data <- y
            x
          })



#' Method setEPMQuery.
#' 
#' Attach (or replace) a user-provided query string 
#' to an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' @param y string (character vector of length 1) corresponding to
#' a PubMed query string.
#' 
#' @rdname setEPMQuery
#' @exportMethod setEPMQuery  
setGeneric("setEPMQuery", function(x, y) {
  standardGeneric("setEPMQuery")
})

#' @rdname setEPMQuery
#' @aliases setEPMQuery,easyPubMed,character-method
setMethod("setEPMQuery", signature(x="easyPubMed", y="character"),
          function(x, y) {
            
            stopifnot(length(y) == 1, 
                      sum(is.na(y)) == 0)
            
            x@query <- y
            x
          })



#' Method setEPMMeta.
#' 
#' Attach (or replace) meta data to an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' @param y list including meta data information.
#' 
#' @rdname setEPMMeta
#' @exportMethod setEPMMeta  
setGeneric("setEPMMeta", function(x, y) {
  standardGeneric("setEPMMeta")
})

#' @rdname setEPMMeta
#' @aliases setEPMMeta,easyPubMed,list-method
setMethod("setEPMMeta", signature(x="easyPubMed", y="list"),
          function(x, y) {
            
            x@meta <- y
            x
          })



#' Method setEPMUilist.
#' 
#' Attach (or replace)  the list of unique record identifiers (PMIDs) 
#' to an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' @param y list of unique PubMed record identifiers (PMIDs).
#' 
#' @rdname setEPMUilist
#' @exportMethod setEPMUilist  
setGeneric("setEPMUilist", function(x, y) {
  standardGeneric("setEPMUilist")
})

#' @rdname setEPMUilist
#' @aliases setEPMUilist,easyPubMed,list-method
setMethod("setEPMUilist", signature(x="easyPubMed", y="list"),
          function(x, y) {
            
            x@uilist <- y
            x
          })



#' Method setEPMRaw.
#' 
#' Attach (or replace) raw PubMed record data 
#' to an `easyPubMed` object. 
#' 
#' @param x an object of class `easyPubMed`. 
#' @param y list of PubMed records (raw data).
#' 
#' @rdname setEPMRaw
#' @exportMethod setEPMRaw  
setGeneric("setEPMRaw", function(x, y) {
  standardGeneric("setEPMRaw")
})

#' @rdname setEPMRaw
#' @aliases setEPMRaw,easyPubMed,list-method
setMethod("setEPMRaw", signature(x="easyPubMed", y="list"),
          function(x, y) {
            
            x@raw <- y
            x
          })



#' Method setEPMMisc.
#' 
#' Attach (or replace)  miscellaneous information 
#' to an `easyPubMed` object. 
#' 
#' @param x an object of class `easyPubMed`.
#' @param y list including miscellaneous data and information. 
#' 
#' @rdname setEPMMisc
#' @exportMethod setEPMMisc  
setGeneric("setEPMMisc", function(x, y) {
  standardGeneric("setEPMMisc")
})

#' @rdname setEPMMisc
#' @aliases setEPMMisc,easyPubMed,list-method
setMethod("setEPMMisc", signature(x="easyPubMed", y="list"),
          function(x, y) {
            
            x@misc <- y
            x
          })



#' Method setEPMJobList.
#' 
#' Attach (or replace) the list of record retrieval sub-jobs 
#' to an `easyPubMed` object. 
#' Record retrieval sub-jobs are stored in a data.frame and each
#' row corresponds to an independent non-overlapping 
#' PubMed query. This `data.frame` guides the record retrieval process. 
#' The `data.frame` is written into the  
#' `misc` slot of an `easyPubMed` object.
#' 
#' @param x an object of class `easyPubMed`. 
#' @param y `data.frame` including the list of 
#' PubMed record retrieaval sub-jobs.
#' 
#' @rdname setEPMJobList
#' @exportMethod setEPMJobList  
setGeneric("setEPMJobList", function(x, y) {
  standardGeneric("setEPMJobList")
})

#' @rdname setEPMJobList
#' @aliases setEPMJobList,easyPubMed,data.frame-method
setMethod("setEPMJobList", signature(x="easyPubMed", y="data.frame"),
          function(x, y) {
            
            tryCatch({
              x@misc$job_list <- y
            },error = function(e) { 
              message(
                paste0(
                  "The job_list data.frame couldn't be ", 
                  "written to this easyPubMed object...")); 
              NULL})
            x
          })



## ------------------------
## --- Advanced Methods ---
## ------------------------


#' Method fetchEPMData.
#' 
#' Retrieve PubMed records for an `easyPubMed` object.
#' 
#' 
#' @param x an easyPubMed-class object.
#' @param params list including parameters to tune the 
#' record retrieval job. For more info, see 
#' `?easyPunMed:::EPM_validate_fetch_params`.
#' 
#' 
#' @rdname fetchEPMData
#' @exportMethod fetchEPMData
setGeneric("fetchEPMData", function(x,params) {
  standardGeneric("fetchEPMData")
})


#' @rdname fetchEPMData
#' @aliases fetchEPMData,easyPubMed,list-method
setMethod("fetchEPMData", 
          signature(x="easyPubMed", params="list"),
          function(x, params) {
            
            # Should keal the job if params are
            # missing or non-compatible 
            my_params <- EPM_validate_fetch_params(params)
            
            # retrieve job list
            job_list <- getEPMJobList(x)
            
            # retrieve object meta
            my_meta <- getEPMMeta(x)
            
            # retrieve misc
            my_misc <- getEPMMisc(x)
            
            
            # prepare outfile information (if requested)
            if (my_params$write_to_file) {
              outfiles <- EPM_prep_outfile(job_list = job_list, 
                                           path = my_params$outfile_path, 
                                           prefix = my_params$outfile_prefix)
              job_list$outfile <- outfiles
            } else {
              job_list$outfile <- NA
            }
            
            # Message if verbose
            if (my_params$verbose)
              message('Retrieving records ', appendLF = FALSE)
            
            # init collectors
            out_recos <- list()
            out_pmids <- list()
            book_pmids <- list()
            
            for (i in seq(1, nrow(job_list), by = 1)) {
              
              tmp_qi <- NULL
              tmp_i <- NULL
              labs_i <- NULL
              cnt_i <- 0
              
              tmp_qi <- paste0('(', job_list$query_string[i], ') AND ("', 
                               job_list$init_date[i], '"[CRDT]:"', 
                               job_list$end_date[i], '"[CRDT])')
              
              while(is.null(tmp_i) && cnt_i < 10) {
                tmp_i <- EPM_retrieve_data(query_string = tmp_qi, 
                                           api_key = my_params$api_key, 
                                           format = my_params$format, 
                                           encoding = my_params$encoding) 
                if (is.null(tmp_i)) {
                  cnt_i <- cnt_i + 1
                  Sys.sleep(1)
                }
              }
              if (is.null(tmp_i)) {
                warning(paste0("Batch ", i, ": no records could be retrieved!"))
                
              } else {
                
                cnt_i <- 0
                
                if (is.list(tmp_i) && length(tmp_i) > 0) {
                  
                  # Parse PMIDs
                  if (my_params$format != 'uilist') {
                    labs_i <- EPM_detect_pmid(tmp_i, format = my_params$format, 
                                              as.list = TRUE)
                  } else {
                    labs_i <- tmp_i
                  }
                  
                  # try detecting Book Chapters
                  if (my_params$format == 'xml') {
                    tryCatch({
                      book_lab_i <- grepl('<BookDocument', 
                                          lapply(tmp_i, substr, 
                                                 start = 1, stop = 30))
                      if (sum(book_lab_i) > 0) {
                        book_pmids[[length(book_pmids) + 1]] <-
                          do.call(c, labs_i[book_lab_i])
                      }
                      
                    }, error = function(e) { NULL })
                  } else {
                    book_pmids[[1]] <- NA
                  }
                  
                  
                  # Assign PMIDs as element names
                  tryCatch({
                    names(tmp_i) <- do.call(c, labs_i)
                  }, error = function(e) { NULL })
                  
                  # Dump into collectors (if requested)
                  out_pmids <- c(out_pmids, labs_i)
                  
                  if (params$store_contents) {
                    out_recos <- c(out_recos, tmp_i)
                  }
                  
                  # Dump into outfile (if requested)
                  if (params$write_to_file && !is.na(job_list$outfile[i])) {
                    
                    # if format is XML, add META section
                    addon_meta <- NULL
                    if (my_params$format == 'xml') {
                      addon_meta <- EPM_encode_meta_to_xml(
                        meta = my_meta, job_list = job_list, 
                        i = i, encoding = my_params$encoding)
                    }
                    
                    # write to file!!!
                    zz <- EPM_write_to_file(x = tmp_i, 
                                            to = job_list$outfile[i], 
                                            format = my_params$format, 
                                            addon = addon_meta,
                                            verbose = my_params$verbose)
                    
                    # Inform if there was an error
                    if (zz != 1 && my_params$verbose){
                      warning(paste0(
                        "An error has occurred while writing PubMed records ", 
                        "to the local disc. Please check the permissions, ", 
                        "make sure there is enough space on the volume, ", 
                        "and then try again."))
                    }
                  }
                }
              }
              if (my_params$verbose)
                message('.', appendLF = FALSE)
            }
            
            # Update meta
            my_meta$raw_format <- my_params$format
            my_meta$raw_encoding <- my_params$encoding
            my_meta$raw_record_num <- tryCatch({length(out_pmids)}, 
                                               error = function(e) { NA })
            my_meta$raw_data_embedded <- as.numeric(my_params$store_contents)
            my_meta$raw_date <- format(Sys.time(), tz = 'GMT')
            
            # Update misc
            book_pmids <- tryCatch({do.call(c, book_pmids)}, 
                                   error = function(e) { NA })
            my_misc$fetch_epm_params <- my_params
            my_misc$book_chapter_pmids <- book_pmids
            
            # update contents
            x <- setEPMMeta(x, my_meta)
            x <- setEPMMisc(x, my_misc)
            x <- setEPMUilist(x, out_pmids)
            x <- setEPMRaw(x, out_recos)
            
            # make sure to reset 
            x <- setEPMData(x, data.frame())
            
            # Close the messaging
            if (my_params$verbose)
              message('', appendLF = TRUE)
            
            # return
            x
          })



#' Method parseEPMData.
#' 
#' Extract, parse and format information from raw PubMed records 
#' stored in an `easyPubMed` object.
#' 
#' @param x an easyPubMed-class object
#' @param params list including parameters to tune the 
#' record data parsing job. For more info, see 
#' `?easyPunMed:::EPM_validate_parse_params`.
#' 
#' @rdname parseEPMData
#' @exportMethod parseEPMData
setGeneric("parseEPMData", function(x, params) {
  standardGeneric("parseEPMData")
})

#' @rdname parseEPMData
#' @aliases parseEPMData,easyPubMed,list-method
setMethod("parseEPMData", signature(x = "easyPubMed", 
                                    params = "list"),
          function(x, params) {
            
            # Should keal the job if params are
            # missing or non-compatible 
            my_params <- EPM_validate_parse_params(params)
            
            # Check if easyPubMed object is compatible (XML, w/ raw data)
            my_meta <- getEPMMeta(x)
            
            # get misc
            my_misc <- getEPMMisc(x)
            
            # Check if the EPM object has data in it
            if (is.na(my_meta$raw_data_embedded) ||
                my_meta$raw_data_embedded != 1 ||
                is.na(my_meta$raw_format) || 
                is.na(my_meta$raw_record_num) ||
                my_meta$raw_record_num < 1) {
              
              stop('This easyPubMed object does NOT contain any raw data.')
              
            } else if (my_meta$raw_format != 'xml') {
              
              stop('Raw data are NOT in xml format.')
            }
            
            # if checks are passed, proceed
            x_raw <- getEPMRaw(x)
            
            # Message if verbose
            if (my_params$verbose)
              message('Parsing record information', appendLF = TRUE)
            
            # prep for messaging (if requested)            
            if (length(x_raw) > 20) {
              prc_seq <- unique(as.integer(
                seq(1, length(x_raw), length.out = 19)))
            } else {
              prc_seq <- length(x_raw)
            }
            
            if (my_params$verbose) {
              vrb_li1 <- '|----+----+----+----| 100%'
              vrb_li2 <- '...................|'
              message(vrb_li1)
              message('|', appendLF = FALSE)
            }
            
            # collector (init) and loop
            y <- list()

            for (i in seq_len(length(x_raw))) {
              
              # message
              if (my_params$verbose && length(prc_seq) > 1 && i %in% prc_seq){
                message('.', appendLF = FALSE)  
              }
              
              # parse
              y[[length(y) + 1 ]] <- tryCatch({
                
                epm_parse_record(x_raw[[i]], 
                                 max_authors = my_params$max_authors, 
                                 autofill_address = my_params$autofill_address, 
                                 compact_output = my_params$compact_output,
                                 include_abstract = my_params$include_abstract,
                                 max_references = my_params$max_references,
                                 ref_id_type = my_params$ref_id_type)
              }, error = function(e) { NULL })
            }
            
            # message
            if (my_params$verbose) {
              if (length(prc_seq) > 1) {
                message('|', appendLF = TRUE)  
              } else {
                message(vrb_li2, appendLF = TRUE)
              }
            }
            
            # collapse
            tot_processed_n <- length(y)
            y <- tryCatch({
              ytmp <- do.call(rbind, y)
              ytmp <- ytmp[!is.na(ytmp[, 1]), ]
              rownames(ytmp) <- NULL
              ytmp
            }, error = function(e) { data.frame() })
            
            if (nrow(y) < 1 && my_params$verbose) {
              message('None of these PubMed records could be parsed.')
            }

            # update misc & meta
            my_misc$parse_epm_params <- my_params
            
            my_meta$processed_compact_output <- my_params$compact_output
            my_meta$processed_record_num <- tot_processed_n
            my_meta$processed_ref_type <- my_params$ref_id_type
            
            # update contents
            x <- setEPMMisc(x, my_misc)
            x <- setEPMMeta(x, my_meta)
            x <- setEPMData(x, y)
            
            # return
            x
          })





## --------------------------------------
## --- Exported Method-wrapping f(x)s ---
## --------------------------------------


#' Get Meta Data from an easyPubMed Object.
#' 
#' Request Meta Data from an `easyPubMed` object. This is a wrapper
#' function that calls the `getEPMMeta()` method. This function returns 
#' contents from the `meta` slot.
#' 
#' @param x An `easyPubMed` object. 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive. 
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   get_epm_meta(x)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#'  
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return a list including meta data from an `easyPubMed` object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' 
#' @export
get_epm_meta <- function(x) {

  stopifnot(inherits(x, "easyPubMed"))
  y <- tryCatch(getEPMMeta(x), error = function(e) { NULL })
  return(y)
  
}




#' Get Raw Data from an easyPubMed Object.
#' 
#' Request Raw Data from an `easyPubMed` object. This is a wrapper
#' function that calls the `getEPMRaw()` method. This function returns contents 
#' from the `raw` slot.
#' 
#' @param x An `easyPubMed` object. 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   x <- epm_fetch(x)
#'   get_epm_raw(x)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#'  
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return a list including raw data from an `easyPubMed` object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' 
#' @export
get_epm_raw <- function(x) {
  
  stopifnot(inherits(x, "easyPubMed"))
  y <- tryCatch(getEPMRaw(x), error = function(e) { NULL })
  return(y)
  
}




#' Get Processed Data from an easyPubMed Object.
#' 
#' Obtain Processed Data that were extracted from a list of PubMed records. 
#' This is a wrapper function that calls 
#' the `getEPMData()` method. This function returns contents 
#' from the `data` slot.
#' 
#' @param x An `easyPubMed` object. 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   x <- epm_fetch(x)
#'   x <- epm_parse(x, max_references = 5, max_authors = 5)
#'   get_epm_data(x)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#'  
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return a `data.frame` including processed data from an `easyPubMed` object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' 
#' @export
get_epm_data <- function(x) {
  
  stopifnot(inherits(x, "easyPubMed"))
  y <- tryCatch(getEPMData(x), error = function(e) { NULL })
  return(y)
  
}




#' Get PubMed Record Identifiers from an easyPubMed Object.
#' 
#' Request the list of unique PubMed Record Identifiers that are contained 
#' in an `easyPubMed` object. This function is a wrapper
#' function calling the `getEPMUilist()` method. This function returns contents 
#' from the `uilist` slot.
#' 
#' @param x An `easyPubMed` object. 
#' 
#' @examples 
#' # Note: a time limit can be set in order to kill the operation when/if 
#' # the NCBI/Entrez server becomes unresponsive.
#' setTimeLimit(elapsed = 4.9)
#' try({
#'   x <- epm_query(query_string = 'Damiano Fantini[AU] AND "2018"[PDAT]')
#'   x <- epm_fetch(x)
#'   get_epm_uilist(x)
#' }, silent = TRUE)
#' setTimeLimit(elapsed = Inf)
#' 
#'  
#' @author 
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @return a character vector including a list of unique record identifiers 
#' from an `easyPubMed` object.
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/easypubmed/}
#' 
#' @export
get_epm_uilist <- function(x) {
  
  stopifnot(inherits(x, "easyPubMed"))
  y <- tryCatch(getEPMUilist(x), error = function(e) { NULL })
  return(y)
  
}







% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epm_all_fx.R
\name{EPM_decode_xml_meta}
\alias{EPM_decode_xml_meta}
\title{Decode an XML String into the Corresponding Metadata.}
\usage{
EPM_decode_xml_meta(x)
}
\arguments{
\item{x}{String corresponding to the XML-decorated text including
metadata from an easyPubMed object/query job.}
}
\value{
String, chunck of XML-decorated text including meta information.
}
\description{
Decode an XML String including a list of meta information 
associated to an easyPubMed object whose contents were written to 
a text file on a local disk. These meta-information are used to keep track 
of easyPubMed query jobs and/or to re-build objects
starting from XML files saved on a local disk.
}
\examples{
xml <- paste0('<EPMxJobData><EPMxJobUniqueId>EPMJ_20231017151112_mi7xvol743', 
              'rvz5ry5z3n8qm0ww</EPMxJobUniqueId><EPMxJobBatchNum>4</EPMxJo', 
              'bBatchNum><EPMxJobBatch>1</EPMxJobBatch><EPMxQuery>Test_Quer', 
              'y</EPMxQuery><EPMxQBatchInitDate>1937/01/22</EPMxQBatchInitD', 
              'ate><EPMxQBatchEndDate>1980/08/01</EPMxQBatchEndDate><EPMxQB', 
              'atchDiffDays>15897</EPMxQBatchDiffDays><EPMxQBatchExpCount>2', 
              '13</EPMxQBatchExpCount><EPMxMaxRecordsPerBatch>1000</EPMxMax', 
              'RecordsPerBatch><EPMxExpCount>2083</EPMxExpCount><EPMxExpNum', 
              'OfBatches>4</EPMxExpNumOfBatches><EPMxAllRecordsCovered>TRUE', 
              '</EPMxAllRecordsCovered><EPMxExpMissedRecords>0</EPMxExpMiss', 
              'edRecords><EPMxQueryDate>2023-10-17 15:11:12</EPMxQueryDate>', 
              '<EPMxRawFormat>xml</EPMxRawFormat><EPMxRawEncoding>UTF-8</EP', 
              'MxRawEncoding><EPMxRawDate>2023-10-17 15:14:12</EPMxRawDate>', 
              '<EPMxLibVersion>3.01</EPMxLibVersion></EPMxJobData>')
easyPubMed:::EPM_decode_xml_meta(xml)



}
\references{
\url{https://www.data-pulse.com/dev_site/easypubmed/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
\keyword{internal}

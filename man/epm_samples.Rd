% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emp_misc.R
\docType{data}
\name{epm_samples}
\alias{epm_samples}
\title{Preprocessed PubMed Records and Data}
\format{
The dataset is formatted as a list including 4 elements: 

* `bladder_cancer_2018`: List of 4

* `bladder_cancer_40y`: List of 1

* `fx`: List of 5
}
\usage{
data("epm_samples")
}
\description{
This dataset includes a collection of sample 
data obtained from PubMed records and saved in different formats.  
This dataset is used to demonstrate specific functionalities of 
the `easyPubMed` R library. Each element in the `epm_samples` list 
corresponds to a different input or intermediate object.
}
\examples{
## Display some contents
data("epm_samples")
# Display Query String used for collecting the data
print(epm_samples$bladder_cancer_2018$demo_data_01)
}
\keyword{datasets}

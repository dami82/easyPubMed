\name{custom_grep}
\alias{custom_grep}
\title{
Retrieve Text Between XML Tags
}
\description{
Extract text form a string containing XML or HTML tags. Text included between tags of interest will be returned. If multiple tagged substrings are found, they will be returned as different elements of a list or character vector. 
}
\usage{
custom_grep(xml_data, tag, format = "list")
}
\arguments{
  \item{xml_data}{
String (of class character and length 1): corresponds to the PubMed record or any string including XML/HTML tags.
}
  \item{tag}{
String (of class character and length 1): the tag of interest (does NOT include < > chars). 
}
  \item{format}{
c("list", "char"): specifies the format for the output.
}
}
\details{
The input string has to be a character string (length 1) containing tags (HTML or XML format). If an XML Document is provided as input, the function will rise an error.
}
\value{
List or vector where each element corresponds to an in-tag substring. 
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
#
# extract substrings based on regular expressions
string_01 <- "I can't wait to watch the <strong>Late Night Show with" 
string_01 <- paste(string_01, "Seth Meyers</strong> tonight at <strong>11:30</strong>pm CT!")
custom_grep(xml_data = string_01, tag = "strong", format = "char")
custom_grep(xml_data = string_01, tag = "strong", format = "list")
}

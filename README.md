# easyPubMed, latest dev version (ver. 3.01)

**easyPubMed** is an **open-source R interface to the Entrez Programming Utilities** aimed at allowing programmatic access to [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/) in the R environment. 
The package is suitable for downloading large number of records, and includes a collection of functions to perform basic 
processing of the *Entrez/PubMed* query responses. The library supports either *XML* or *TXT* ("medline") format.

---



### New features of easyPubMed version 3.01

- **Simplified Pipeline**. The process of retrieving and analyzing Pubmed records 
has been updated and simplified. The revised pipeline includes three steps: 
*1)* submit a query;  *2)* fetch records;  *3)* extract information. The 
corresponding functions are discussed below in this vignette.


- **Automatic Job splitting into Sub-Queries**. The Entrez server imposes a strict 
n=10,000 limit to the number of records that can be programmatically 
retrieved from a single query. Whenever possible, the `easyPubMed` 
library automatically attempts to split queries returning large number of
records into lists of smaller, manageable queries.  


- **The easyPubMed S4 Class**. Here we introduce a new S4 class 
(`easyPubMed`) that was designed to better store and manage query information, 
retrieved records and associated meta-data. 
This S4 class comes with a series of methods and is 
aimed at improving data handling and analysis reproducibility. Raw or processed
data can be obtained from an `easyPubMed` object via the appropriate
*getter* functions (as discussed below).


- **Additional Parsed Fields**. The new `epm_parse()` function supports 
extraction of additional information from Pubmed records (compared to 
previous versions of this R library). Extracted fields now include
*mesh_codes*, *grant_ids*, *references* and *conflict of interest statements* 
(cois) among others. For more information, see the examples below. 


- **Compact Output**. Unlike previous versions of `easyPubMed`, it is now 
possible to collapse author information (*i.e.*, author names) into a 
single string. This way, the output (parsed `data.frame`) only 
includes one row per record.

---

### Installation

The latest version (3.01) of the library is hosted on *GitHub*, and you can install it using the `devtools` R library as follows.

```
devtools::install_github("dami82/easyPubMed")
```

---

### Notes & Info

- `easyPubMed` is an open-source software, under the GPL-3 license and comes with ABSOLUTELY NO WARRANTY. For more questions about the GPL-3 license terms, see [www.gnu.org/licenses](https://www.gnu.org/licenses/gpl-3.0.html).

- This R library was written based on the information included in the *Entrez Programming Utilities Help* manual authored by Eric Sayers, PhD and available on the NCBI Bookshelf ([NBK25500](https://www.ncbi.nlm.nih.gov/books/NBK25500/)).

- This R library is NOT endorsed, supported, maintained NOR affiliated with NCBI.

- There is only one person maintaining this R package: I work on code updates in my spare time and for free. Take-home message: please, be patient.

 

  

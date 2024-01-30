# RBPInper
RNA-binding protein interaction mapper

# Overview
RBPInper is a simple two-step method that extends the analysis of high-throughput RBP’s RNA binding activities. It runs on two input files, including a table of P-values with genes in rows and datasets in columns. The columns can include P-values of differential gene expression from RBP perturbations (knockdown, overexpression), signal enrichment from RBP immunoprecipitation (RIP-Seq, eCLIP, PAR-CLIP). The second input is a dataset annotation file that indicates the cell type grouping and the data experimental strategy. Depending on the project, the cell type grouping may also be used to specify many other groups such as cell or tissue lineages. Experimental strategies are automatically filtered to remove non-RNA strategies. The R implementation of the framework enables flexible inclusion of new integration methods and a set of optional parameters at run time. When available, the final output has two tables representing global and group-specific interactomes (see Materials and Methods).

![workflow](https://github.com/caanene1/RBPInper/blob/5e39afdf71e8fa5faa4515ddc0d245a44bd8f4f6/Workflow.jpg)



This tool is built in R and uses FMStable, GenomicRanges, rtracklayer, dplyr
These applications are automatical installed with RBPInper.

# Installation and running the tool
The best way to get RBPInper along with all the dependencies is to install the release from github at [github](https://github.com/AneneLab/RBPInper) with:

``` r
devtools::install_github('AneneLab/RBPInper')
```
or 

``` r
devtools::install_github('caanene1/RBPInper')
```

This will add two core functions as below:

| Function | Context | Usage |
| ---    | --- | --- |
| rbpinper.run | Run integration analysis | ```?rbpinper.run``` |
| prebed | Prepare bed files (optional) | ```?prebed``` |

# Input rbpinper.run()
Evidence Matrix - The p-value tables, specified by ```evi```, where genes in rows and datasets in columns.
Make sure the p-values are numeric and all less than 1.

| gene  | data1 | data2 | data3 |
| --- | --- | --- | --- |
| TGF | 0.0002 | 0.04 | 1 | 
| PDL1  |  0.1  | 0.056  |  0.6  |
| CD21  |  0.8  | 0.4  |  0.2 |

This input can contain meta information but they must be non numeric.

Information matrix - The dataset annotation table, specified by ```info```, where datasets in rows and annotations in column. The first two columns of this table should be ID and cell group. 

| ID  | cell types |
| --- | --- |
| Foxp1 | HepG2 |
| PD1  |  PM  |
| CD8  |  HepG2  |

You can also supply other type of grouping instead of cell types.


# Output rbpinper.run()
RBPInper object

| ID | Information |
| --- | --- |
| evidence | P-values from evi without the meta information |
|  meta | Meta information from evi argument (None numeric columns) |
|  cells  | Unique samples from info argument |
|  gene_id  | The column name of the gene ID in evi |
|  id | The column name of the dataset ID in info argument |
|  information |  List of integration groups |
| L1.result | Interactome at cell group level |
| L2.result | Global interactome |


# Input prebed()
Bed files - List of bed file file paths, specified by ```bed``` , the function can take a list or vector of bed file paths to increase speed.

GFF, GTF, GTF3 or GFF3 - Genome annotation file path, specified by ```gtf``` , this file will be used to process all the bed files specified above.


# Output prebed()
Data frame of P-value evidence for each bed file.

# Extras
R codes to reproduce the downstream analyses reported in the paper are inside the folder "inst/Manuscript". The file manuscript.R has everything code used in the manuscript.

# Test Datasets
Datasets to test the package and replicate the analysis in the manuscript are inside the folder "inst/data"


# Citation
Meta-Analysis of RNA Binding Protein Interaction profiles using the RBPInper tool
Joseph A. Cogan, Kuklinkova Rene, Benova Natalia, James R. Boyne, Chinedu A. Anene

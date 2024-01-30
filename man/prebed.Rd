\name{prebed}
\alias{prebed}

\title{
Prepare bed file for RBPInper
}
\description{
Annotate narrow peak bed files using gtf,gff,gff3 annotation file
}
\usage{
prebed(bed, gtf, fTarget = "gene", islog10 = TRUE, mMerge = "Bonferroni")
}

\arguments{
  \item{bed}{
Path to bed file to annotate and summaries. It can take a vector of paths and process them together for speed optimisation. This must have the narrow peaks columns from MACS2 or similar.
}
  \item{gtf}{
Path to gtf, gff, gff3 annotation file. This should the same genome version as bed file.
}
  \item{fTarget}{
The gtf feature to use for summarising the peaks, default "gene". It can also use exon, transcript, etc but the gtf used must have the "type"" column.
}
  \item{islog10}{
Is bed file p-value -log10 transformed, default TRUE. This should ideally be the case for narrow peak file.
}
  \item{mMerge}{
Method to merge peak p-values to gene p-values. Avilable methods include "Bonferroni" default, "Fisher", "Binomial" or "Harmonic mean". The defualt is the most robust for combining peak p-values due to dependency. The implimention of the fusion function allows for lower or upper case, short version and mispelling provided the first four letters are correct.
}
}
\details{
This function can operate on any source of peak file provided the matching gtf is provided and the standard bed columns present.
}
\value{
It returns a data frame with correct column names, ensemble id and gene symbols
\item{output }{Data frame of genes in rows and p-value in column}
}
\references{
Depends on:
 base R functions.
 Uses the R packages GenomicRanges, rtracklayer and dplyr
}
\author{
Chinedu A. Anene, PhD
}
\note{
This function expects that bed files and gtf files are matched.
}

\seealso{
The output can be used for a varity of downstream analysis and directly suitable for other RBPInper functions.
}
\examples{
 See the mansucript
  }


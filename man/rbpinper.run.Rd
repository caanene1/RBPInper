\name{rbpinper.run}
\alias{rbpinper.run}

\title{
Integrate RBP interaction profiles.
}
\description{
Meta-analysis of RNA binding protein interaction profiles.
}
\usage{
rbpinper.run(evi, info, gene.col = "gene_id", info.id = "ID", cell.col = "cell", L1 = "bonfe", L2 = "fishe", p = 0.05)
}

\arguments{
  \item{evi}{
Data frame of p-value evidence, genes in row and data in columns.
}
  \item{info}{
Data frame of information file, at least data ID and cell type included.
}
  \item{gene.col}{
Column name in "evi" p-value matrix with gene IDs.
}
  \item{info.id}{
Column name in "info" information matrix with matching sample IDs.
}
  \item{cell.col}{
Column name in "info" information matrix with cell type name.
}
  \item{L1}{
Method for p-value merging at cell type level, default "bonf".
}
  \item{L2}{
Method for p-value merging at global level, default "fisher".
}
  \item{p}{
P-value threshold for calls.
}
}
\details{
RBPInper object with cell type level results and global results.
}
\value{
RBPInper output.
}
\references{
Meta-Analysis of RNA Binding Protein Interaction profiles using the RBPInper tool.
}
\author{
Chinedu A. Anene
}
\note{
%%  ~~further notes~~
}

\seealso{
prebed, fusion
}
\examples{
 rbpinper.run(evi=predata, info = infoRNA)
}

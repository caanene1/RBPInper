\name{fs_matrix}
\alias{fs_matrix}

\title{
P-value fusion on data frame.
}
\description{
Apply fusion on data frame.
}
\usage{
fs_matrix(x, alpha = 0.05, type = "fisher")
}

\arguments{
  \item{x}{
Data frame of p-values to combine.
}
  \item{alpha}{
Alpha method to use for the Binomial method.
}
  \item{type}{
Method to use for merging, fisher, binomial, bonferroni, harmonic.
}
}
\details{
Derives directly from the fusion function.
}
\value{
Data frame of fused p-values with threshold labels.
}
\references{

}
\author{
Chinedu A. Anene, PhD
}
\note{

}

\seealso{
fusion
}
\examples{
  }

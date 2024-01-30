\name{fusion}
\alias{fusion}

\title{
P-value fusion
}
\description{
Combine multiple p-values to a single p-value.
}
\usage{
fusion(x = c(0.02,0.08,0.001), type = "all", alpha = 0.05).
}

\arguments{
  \item{x}{
Vector of p-values to combine.
}
  \item{alpha}{
Alpha method to use for the Binomial method.
}
  \item{type}{
Character to indicate if p-values in x are dependent/independent.
Defaults to "all", but can be set to Fisher, Binomial, Bonferroni,
Harmonic (all lower case letters will work. Option "both" does not
combine the methods, returns them, so it is possible to do filtering after.
}
}

\details{
This function provides a range of methods to combine multiple p-values.
}

\value{
A named list of combined p-values corresponding to the avilable methods.
 \item{Fisher}{The clasical Fisher method Ref:dio}
 \item{Brown}{The extension of fisher method by Brown Ref:dio}
}

\references{
Depends on:
 base R functions.
 Calcualtes the methods:
 Fisher
 Bionomial
 Bonferroni
 Harmonic mean

}

\author{
Chinedu A. Anene, PhD
}
\note{
This function expects your p-values to come from clasical
omics analysis tools, such Limma, edger macs2 etc.

}

\seealso{
The output of this function is intended for RBPInper, but can be used for other things.
}

\keyword{ ~Combine p-values, multi-omics integration }

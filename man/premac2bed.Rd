\name{premac2bed}
\alias{premac2bed}

\title{
Prepare macs2 bed file and save it to file.
}
\description{
Prepare macs2 peak call file with different column arrangement
}
\usage{
premac2bed(fie, nam="fie.bed", soure="base")
}

\arguments{
  \item{fie}{
Data frame of the bed columns.
}
  \item{nam}{
The name of the saved bed file
}
  \item{soure}{
Indicate the source of the bed file to allow column changes.
Either "base" for standard macs2 output or "encode" for encode qvalue call output.
}

}
\details{
This function is to process bed file from macs2 with a different column
arrnagment compared to the narrow peak bed file from Encode.
}
\value{
Bed file saved to the current working directory.
}
\references{
}
\author{
Chinedu A. Anene, PhD
}
\note{
Check that input has the expected other.
}


\seealso{
bedload
}
\examples{
}


\name{migration.cmr}
\alias{migration.cmr}
\title{Crude Migration Rate}
\usage{
  migration.cmr(m, PAR, k = 100)
}
\arguments{
  \item{m}{migration matrix}

  \item{PAR}{population at risk (estimated average
  population size)}

  \item{k}{scaling constant (set to \code{100} by default
  to result in percentage)}
}
\value{
  percentage (when \code{k=100})
}
\description{
  Crude Migration Rate
}
\examples{
data(migration.world)
migration.cmr(migration.world, 6e+9)
}
\references{
  \itemize{ \item Philip Rees, Martin Bell, Oliver
  Duke-Williams and Marcus Blake (2000) Problems and
  Solutions in the Measurement of Migration Intensities:
  Australia and Britain Compared. \emph{Population Studies}
  \bold{54}, 207--222 }
}


\name{check.migration.matrix}
\alias{check.migration.matrix}
\title{Check Migration Matrix}
\usage{
  check.migration.matrix(m)
}
\arguments{
  \item{m}{R object to check}
}
\value{
  (invisibly) TRUE
}
\description{
  Checks if provided R object looks like a migration
  matrix.
}
\details{
  A migration matrix is a symmetric matrix with the exact
  same row and column names. The diagonal equals to zero.
  The upper triangle shows the in- and the lower triangle
  shows the out-migration.
}
\keyword{internal}


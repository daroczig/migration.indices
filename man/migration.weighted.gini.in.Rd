\name{migration.weighted.gini.in}
\alias{migration.weighted.gini.in}
\title{Migration-weighted In-migration Gini Index}
\usage{
  migration.weighted.gini.in(m, mgi)
}
\arguments{
  \item{m}{migration matrix}

  \item{mgi}{optionally passed (precomputed) Migration
  In-migration Gini Index}
}
\value{
  number
}
\description{
  Migration-weighted In-migration Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.weighted.gini.in(migration.hyp)   # 0.1222222
migration.weighted.gini.in(migration.hyp2)  # 0.05238095
}
}
\references{
  \itemize{ \item M. Bell, M. Blake, P. Boyle, O.
  Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002)
  Cross-National Comparison of Internal Migration. Issues
  and Measures. \emph{Journal of the Royal Statistical
  Society. Series A (Statistics in Society)} \bold{165},
  435--464 }
}
\seealso{
  \code{\link{migration.gini}}
  \code{\link{migration.gini.in}}
  \code{\link{migration.weighted.gini.out}}
  \code{\link{migration.weighted.gini.mean}}
}


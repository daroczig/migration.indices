\name{migration.weighted.gini.mean}
\alias{migration.weighted.gini.mean}
\title{Migration-weighted Mean Gini Index}
\usage{
  migration.weighted.gini.mean(m, mwgi, mwgo)
}
\arguments{
  \item{m}{migration matrix}

  \item{mwgi}{optionally passed (precomputed)
  Migration-weighted In-migration Gini Index}

  \item{mwgo}{optionally passed (precomputed)
  Migration-weighted Out-migration Gini Index}
}
\value{
  number
}
\description{
  Migration-weighted Mean Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.weighted.gini.mean(migration.hyp)  # 0.06111111
migration.weighted.gini.mean(migration.hyp2) # 0.03660714
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
  \code{\link{migration.weighted.gini.in}}
  \code{\link{migration.weighted.gini.out}}
}


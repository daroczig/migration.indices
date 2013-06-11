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
  This combined index results in a number between 0 and 1
  where 0 means no spatial focusing and 1 shows maximum
  focusing.
}
\description{
  The Migration-weighted Mean Gini Index is simply the
  average of the Migration-weighted In-migration
  (\code{\link{migration.weighted.gini.in}}) and the
  Migration-weighted Out-migration
  (\code{\link{migration.weighted.gini.out}}) Gini Indices:
  \deqn{MWG^A = \frac{MWG^O + MWG^I}{2}}
}
\examples{
data(migration.hyp)
migration.weighted.gini.mean(migration.hyp)  # 0.06111111
migration.weighted.gini.mean(migration.hyp2) # 0.03660714
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


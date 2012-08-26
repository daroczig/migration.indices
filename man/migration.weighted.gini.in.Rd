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
\author{
  Gergely Daróczi
}
\references{
  M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees,
  J. Stillwell and G. Hugo (2002): Cross-National
  Comparison of Internal Migration. Issues and Measures.
  In. Journal of the Royal Statistical Society. Series A
  (Statistics in Society), Vol. 165, No. 3 (2002), pp.
  435-464
}
\seealso{
  \code{\link{migration.weighted.gini.out}}
  \code{\link{migration.weighted.gini.mean}}
}


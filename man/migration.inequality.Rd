\name{migration.inequality}
\alias{migration.inequality}
\title{Migration Inequality Index}
\usage{
  migration.inequality(m,
    expected = c("equal", "weighted"))
}
\arguments{
  \item{m}{migration matrix}

  \item{expected}{type of expected distribution}
}
\value{
  number between 0 and 1 where 1 shows greater inequality
}
\description{
  Measures the distance from an expected distribution.
}
\examples{
\dontrun{
data(migration.hyp)
migration.inequality(migration.hyp)
migration.inequality(migration.hyp, expected = 'weighted')
data(migration.world)
migration.inequality(migration.world)
)}
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


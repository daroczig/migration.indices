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
migration.weighted.gini.mean(migration.hyp)   #
migration.weighted.gini.mean(migration.hyp2)  #
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
  \code{\link{migration.weighted.gini.in}}
  \code{\link{migration.weighted.gini.out}}
}


\name{migration.weighted.gini.out}
\alias{migration.weighted.gini.out}
\title{Migration-weighted Out-migration Gini Index}
\usage{
  migration.weighted.gini.out(m, mgo)
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
  Migration-weighted Out-migration Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.weighted.gini.out(migration.hyp)   #
migration.weighted.gini.out(migration.hyp2)  #
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
  \code{\link{migration.weighted.gini.mean}}
}


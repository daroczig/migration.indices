\name{migration.connectivity}
\alias{migration.connectivity}
\title{Migration Conncetivity Index}
\usage{
  migration.connectivity(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number between 0 and 1 where zero shows no connections
  between regions
}
\description{
  Measures the proportion of the total number of potential
  interregional flows which are not zero.
}
\examples{
\dontrun{
data(migration.hyp)
migration.connectivity(migration.hyp)
data(migration.world)
migration.connectivity(migration.world)
}
}
\references{
  M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees,
  J. Stillwell and G. Hugo (2002) Cross-National Comparison
  of Internal Migration. Issues and Measures. \emph{Journal
  of the Royal Statistical Society. Series A (Statistics in
  Society)} \bold{165}, 435--464
}


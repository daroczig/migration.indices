\name{migration.effectiveness}
\alias{migration.effectiveness}
\title{Migration Effectiveness Index}
\usage{
  migration.effectiveness(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number between 0 and 100 where the higher number shows an
  efficient mechanism of population redistribution
}
\description{
  Measures the degree of (a)symmetry or (dis)equilibrium in
  the network of interregional migration flows.
}
\examples{
\dontrun{
data(migration.hyp)
migration.effectiveness(migration.hyp)
data(migration.world)
migration.effectiveness(migration.world)
}
}
\references{
  \itemize{ \item Martin Bell and Salut Muhidin (2009)
  {Cross-National Comparisons of Internal Migration}.
  Research Paper. UNDP.
  \url{http://hdr.undp.org/en/reports/global/hdr2009/papers/HDRP_2009_30.pdf}
  }
}


#' Migration Effectiveness Index
#'
#' Measures the degree of (a)symmetry or (dis)equilibrium in the network of interregional migration flows.
#' @param m migration matrix
#' @return number between 0 and 100 where the higher number shows an efficient mechanism of population redistribution
#' @references Martin Bell and Salut Muhidin (2009) {Cross-National Comparisons of Internal Migration}. Research Paper. UNDP. \url{http://hdr.undp.org/en/reports/global/hdr2009/papers/HDRP_2009_30.pdf}
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.effectiveness(migration.hyp)
#' data(migration.world)
#' migration.effectiveness(migration.world)
#' }
#' @export
migration.effectiveness <- function(m) {

    check.migration.matrix(m)

    D <- colSums(m, na.rm = TRUE)
    O <- rowSums(m, na.rm = TRUE)

    100 * sum(abs(D - O)) / sum(D + O)

}


#' Migration Conncetivity Index
#'
#' Measures the proportion of the total number of potential interregional flows which are not zero.
#' @param m migration matrix
#' @return number between 0 and 1 where zero shows no connections between regions
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.connectivity(migration.hyp)
#' data(migration.world)
#' migration.connectivity(migration.world)
#' }
#' @export
migration.connectivity <- function(m) {

    check.migration.matrix(m)

    diag(m) <- NA
    n       <- nrow(m)
    sum(m != 0, na.rm = TRUE) / (n * (n - 1))

}


#' Migration Inequality Index
#'
#' Measures the distance from an expected distribution.
#' @param m migration matrix
#' @param expected type of expected distribution
#' @return number between 0 and 1 where 1 shows greater inequality
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.inequality(migration.hyp)
#' migration.inequality(migration.hyp, expected = 'weighted')
#' data(migration.world)
#' migration.inequality(migration.world)
#' )}
#' @export
migration.inequality <- function(m, expected = c('equal', 'weighted')) {

    expected <- match.arg(expected)
    check.migration.matrix(m)

    diag(m) <- NA
    n       <- nrow(m)

    if (expected == 'equal') {
        m.expected       <- matrix(sum(m, na.rm = TRUE) / (n^2 - n), n, n)
        diag(m.expected) <- NA
    } else {
        m.expected <- m
        for (i in 1:nrow(m))
            m.expected[i, -i] <- colSums(m.expected[, -i], na.rm = TRUE) * rowSums(m, na.rm = TRUE)[i] / sum(colSums(m.expected[, -i], na.rm = TRUE))
    }

    sum(abs(m - m.expected), na.rm = TRUE) * 0.5

}

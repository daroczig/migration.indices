#' Migration Effectiveness Index
#'
#' Measures the degree of (a)symmetry or (dis)equilibrium in the network of interregional migration flows.
#' @param m migration matrix
#' @return number between 0 and 100 where the higher number shows an efficient mechanism of population redistribution
#' @references Martin Bell and Salut Muhidin (2009): Cross-National Comparisons of Internal Migration. Research Paper 2009/30. UNDP.
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.effectiveness(migration.hyp)
#' data(migration.world)
#' migration.effectiveness(migration.world)
#' }
#' @author Gergely Daróczi
#' @export
migration.effectiveness <- function(m) {

    check.migration.matrix(m)

    D <- colSums(m, na.rm = TRUE)
    O <- rowSums(m, na.rm = TRUE)

    100 * sum(abs(D - O)) / sum(D + O)

}


#' Migration Conncetivity Index
#'
#' Meauseres the proportion of the total number of potential interregional flows which are not zero.
#' @param m migration matrix
#' @return number between 0 and 1 where zero shows no connections between regions
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002): Cross-National Comparison of Internal Migration. Issues and Measures. In. Journal of the Royal Statistical Society. Series A (Statistics in Society), Vol. 165, No. 3 (2002), pp. 435-464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.connectivity(migration.hyp)
#' data(migration.world)
#' migration.connectivity(migration.world)
#' }
#' @author Gergely Daróczi
#' @export
migration.connectivity <- function(m) {

    check.migration.matrix(m)

    diag(m) <- NA
    n       <- nrow(m)
    sum(m != 0, na.rm = TRUE) / (n * (n - 1))

}


#' Migration Inequality Index
#'
#' Meauseres the proportion of the total number of potential interregional flows which are not zero.
#' @param m migration matrix
#' @return number between 0 and 1 where 1 shows greater inequality
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002): Cross-National Comparison of Internal Migration. Issues and Measures. In. Journal of the Royal Statistical Society. Series A (Statistics in Society), Vol. 165, No. 3 (2002), pp. 435-464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.inequality(migration.hyp)
#' data(migration.world)
#' migration.inequality(migration.world)
#' )}
#' @author Gergely Daróczi
#' @export
## migration.inequality <- function(m) {

##     check.migration.matrix(m)

##     diag(m)     <- NA
##     m.expected  <- ???  ## TODO: what is that expected matrix? p. 454
##     sum(m - m.expected na.rm = TRUE) * 0.5

## }

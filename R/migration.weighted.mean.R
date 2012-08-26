################################################################################
## Based on:
##  * M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002): Cross-National Comparison of Internal Migration: Issues and Measures. In. Journal of the Royal Statistical Society. Series A (Statistics in Society), Vol. 165, No. 3 (2002), pp. 435-464
################################################################################

#' .. content for
#'
#' .. content for
#' @param m
#' @references Rogers, A. and Raymer, J. (1998) The spatial focus of US interstate migration flows. Int. J PoplnGeogr,4, 63-80.
migration.gini.in.standardized <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    apply(m, 2, function(m.row) (sum(dist(m.row), na.rm = TRUE)*2)/(2*(n-2)*colSums(m, na.rm = TRUE)))


}

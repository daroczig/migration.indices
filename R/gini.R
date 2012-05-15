#' Check Migration Matrix
#'
#' Checks if provided R object looks like a migration matrix.
#'
#' A migration matrix is a symmetric matrix with the exact same row and column names. The diagonal equals to zero. The upper triangle shows the in- and the lower triangle shows the out-migration.
#' @param m R object to check
#' @return (invisibly) TRUE
#' @author Gergely Daróczi
#' @keywords internal
check.migration.matrix <- function(m) {

    ## dummy checks on provided matrix
    if (missing(m))
        stop('No data provided!')
    if (!is.matrix(m))
        stop('Wrong data type (!matrix) provided!')
    if (nrow(m) != ncol(m))
        stop('Wrong data tpye (!symmetrical matrix) provided!')
    if (!is.numeric(m))
        stop('Wrong data tpye (!numeric) provided!')

    return(invisible(TRUE))

}


#' Total Flows Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0), 3, 3)
#' migration.gini.total(m)        # 0.2222222
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.total(m2)       # 0.1875
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.total <- function(m) {

    check.migration.matrix(m)

    n <- nrow(m)
    ## get values from matrix except for the diagonal
    m.val       <- m[xor(upper.tri(m), lower.tri(m))]
    m.val.l     <- length(m.val)
    ## vectorized solution (cannot vectorize all because of memory limits)
    return(sum(apply(as.data.frame(m.val), 1, function(x) sum(abs(m.val-x))))/(2*n*(n-1)*sum(m)))

    ## big memory method
    return(sum(abs(rep(m.val, m.val.l) -unlist(lapply(m.val, rep, m.val.l))))/(2*n*(n-1)*sum(m)))

    ## original (loop) method
    for (i in l)
        for (j in setdiff(l, i))
            for (g in l)
                for (h in setdiff(l, g))
                    res <- res + abs(m[i, j] - m[g, h])
    res/(2*n*(n-1)*sum(m))

}


#' Rows Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0), 3, 3)
#' migration.gini.row(m)  # 0
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.row(m2) # 0.02083333
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.row <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    res <- sum(apply(m, 1, function(m.row) sum(dist(m.row), na.rm = TRUE)*2))
    return(res/(2*n*(n-1)*sum(m, na.rm = TRUE)))

    ## vectorized old solution
    res <- sum(apply(m, 1, function(m.row) {
        id <- !apply(m, 1, function(x) identical(x, m.row))
        sum(apply(as.data.frame(m.row[id]), 1, function(m.cell) sum(abs(m.row - m.cell), na.rm = TRUE)))
    }))

    ## original (loop) method
    n   <- nrow(m)
    l   <- 1:n
    res <- 0
    for (i in l)
        for (j in setdiff(l, i))
                for (h in setdiff(l, c(i, j)))
                    res <- res + abs(m[i, j] - m[i, h])
    res/(2*n*(n-1)*sum(m))

}


#' Standardized Rows Gini Index
#'
#' @param m migration matrix
#' @param migration.gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.row.standardized(m)     # 0
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.row.standardized(m2)    # 11.11111
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.row.standardized <- function(m, migration.gini.total = migration.gini.total(m)) {

    100 * migration.gini.row(m) / migration.gini.total
}


#' Columns Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.col(m)  # 0.05555556
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.col(m2) # 0.04166667
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.col <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    res <- sum(apply(m, 2, function(m.row) sum(dist(m.row), na.rm = TRUE)*2))
    return(res/(2*n*(n-1)*sum(m, na.rm = TRUE)))

    ## original (loop) method
    n   <- ncol(m)
    l   <- 1:n
    res <- 0
    for (j in l)
        for (i in setdiff(l, j))
                for (g in setdiff(l, c(i, j)))
                    res <- res + abs(m[i, j] - m[g, j])
    res/(2*n*(n-1)*sum(m))

}


#' Standardized Columns Gini Index
#'
#' @param m migration matrix
#' @param migration.gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.col.standardized(m)     # 25
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.col.standardized(m2)    # 22.22222
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.col.standardized <- function(m, migration.gini.total = migration.gini.total(m)) {

    100 * migration.gini.col(m) / migration.gini.total
}


#' Exchange Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.exchange(m)     # 0.05555556
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.exchange(m2)    # 0.04166667
#' }
#' @author Gergely Daróczi
migration.gini.exchange <- function(m) {

    check.migration.matrix(m)

    n           <- nrow(m)
    m.t         <- t(m)
    m.val       <- m[xor(upper.tri(m), lower.tri(m))]
    m.t.val     <- m.t[xor(upper.tri(m.t), lower.tri(m.t))]

    return(sum(abs(m.val - m.t.val))/(2*n*(n-1)*sum(m)))

    ## original (loop) method
    n   <- nrow(m)
    l   <- 1:n
    res <- 0
    for (i in l)
        for (j in setdiff(l, i))
            res <- res + abs(m[i, j] - m[j, i])
    res/(2*n*(n-1)*sum(m))

}


#' Standardized Exchange Gini Index
#'
#' @param m migration matrix
#' @param migration.gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.exchange.standardized(m) # 25
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.exchange.standardized(m2) # 22.22222
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.exchange.standardized <- function(m, migration.gini.total = migration.gini.total(m)) {

    100 * migration.gini.exchange(m) / migration.gini.total

}


#' Out-migration Field Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.out(m)  # 0 0 0
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.out(m2) # 0.000 0.125 0.000
#' }
#' @author Gergely Daróczi
migration.gini.out <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    return(apply(m, 1, function(m.row) sum(dist(m.row), na.rm = TRUE)*2)/(2*(n-1)*rowSums(m, na.rm = TRUE)))

    ## original (loop) solution
    n   <- nrow(m)
    l   <- 1:n
    res <- NULL
    for (k in l) {
        res.t <- 0
        for (j in setdiff(l, k))
            for (h in setdiff(l, k))
                res.t <- res.t + abs(m[k, j] - m[k, h])
        res <- c(res, res.t/(2*(n-1)*rowSums(m)[k]))
    }

    return(res)

}


#' In-migration Field Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997): Measuring Spatial Focusing in a Migration System. In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini.in(m)   # 0.1000000 0.2500000 0.1666667
#' m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
#' migration.gini.in(m2)  # 0.1000000 0.0000000 0.2142857
#' }
#' @author Gergely Daróczi
#' @export
migration.gini.in <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    return(apply(m, 2, function(m.row) sum(dist(m.row), na.rm = TRUE)*2)/(2*(n-1)*colSums(m, na.rm = TRUE)))

    ## original (loop) solution
    n   <- nrow(m)
    l   <- 1:n
    res <- NULL
    for (k in l) {
        res.t <- 0
        for (i in setdiff(l, k))
            for (g in setdiff(l, k))
                res.t <- res.t + abs(m[i, k] - m[g, k])
        res <- c(res, res.t/(2*(n-1)*colSums(m)[k]))
    }

    return(res)

}


#' Spatial Gini Indexes
#'
#' @param m migration matrix
#' @author Gergely Daróczi
#' @examples \dontrun{
#' m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
#' migration.gini(m)
#' }
#' @export
migration.gini <- function(m) {

    res <- list(
             migration.gini.total         = migration.gini.total(m),
             migration.gini.exchange      = migration.gini.exchange(m),
             migration.gini.row           = migration.gini.row(m),
             migration.gini.col           = migration.gini.col(m),
             migration.gini.in            = migration.gini.in(m),
             migration.gini.out           = migration.gini.out(m)
      )
    res$migration.gini.row.standardized           <- migration.gini.row.standardized(m, res$migration.gini.total)
    res$migration.gini.col.standardized           <- migration.gini.col.standardized(m, res$migration.gini.total)
    res$migration.gini.exchange.standardized      <- migration.gini.exchange.standardized(m, res$migration.gini.total)

    class(res) <- 'migration.gini'
    return(res)
}


#' @method print migration.gini
#' @S3method print migration.gini
print.migration.gini <- function(x) {

    cat('\n')
    cat('Total Flows Gini Index:\t\t\t', x$migration.gini.total, '\n')
    cat('Rows Gini Index:\t\t\t', x$migration.gini.row, '\n')
    cat('Standardized Rows Gini Index:\t\t', x$migration.gini.row.standardized, '\n')
    cat('Columns Gini Index:\t\t\t', x$migration.gini.col, '\n')
    cat('Standardized Columns Gini Index:\t', x$migration.gini.col.standardized, '\n')
    cat('Exchange Gini Index:\t\t\t', x$migration.gini.exchange, '\n')
    cat('Standardized Exchange Gini Index:\t', x$migration.gini.exchange.standardized, '\n')
    cat('In-migration Field Gini Index:\t\t', 'vector', '\n')
    cat('Out-migration Field Gini Index:\t\t', 'vector', '\n')
    cat('\n')

}

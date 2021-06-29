#' Ensure the working directory contains desired subdirectories
#' 
#' For each desired subdirectory, the function checks if the subdirectory
#' already exists, and if it does not, creates it.
#' 
#' @param needed_dirs A character vector giving the desired subdirectories
#' 
#' @return Invisibly returns a character vector of the directories the function
#'     created because they were not present
prepare_directories <- function(needed_dirs = c("output", "plots")) {
    created_dirs <- character()
    for ( directory in needed_dirs ) {
        if ( dir.exists(directory) ) {
            next()
        } else {
            dir.create(directory)
            created_dirs <- c(created_dirs, directory)
        }
    }
    return(invisible(created_dirs))
}


#' Determine if a roll call vote is "lopsided"
#' 
#' @param x An integer vector of responses
#' @param cutoff A numeric vector of length one; what is the smallest allowable
#'     proportion of votes that can be minority votes before we call the vote
#'     lopsided?
#' @return A logical vector of length one; if the proportion of minority votes
#'     in x is less than the value supplied in cutoff, it is TRUE, and
#'     otherwise it is FALSE.
is_lopsided <- function(x, cutoff = 0.01) { ## helper func: are votes lopsided?
    tab  <- table(x)            ## Tabulate the votes
    prop <- min(tab) / sum(tab) ## Get proportion of minority votes
    return(prop < cutoff)       ## Return TRUE if that proporion is < cutoff
}

#' Determine if a roll call vote is unanimous
#' 
#' @param x An integer vector of responses
#' @return A logical vector of length one; if there is only one unique non-NA
#'     value in x, it is TRUE, and otherwise it is FALSE. 
is_unanimous <- function(x) {
    return(length(unique(na.omit(x))) == 1)
}

#' Determine if a legislator has too few votes to analyze
#' 
#' @param x An integer vector of responses
#' @param allowed_missingness A numeric vector of length one; what is the
#'     largest allowable proportion of missing votes before we say the
#'     legislator has too few votes?
#' @return A logical vector of length one; if the proportion of missing votes
#'     in x is greater than the value supplied in allowed_missingness, it is
#'     TRUE, and otherwise it is FALSE.
has_few_votes <- function(x, allowed_missingness = 0.9) {
    return(sum(is.na(x)) > (allowed_missingness * length(x)))
}


## This function creates a plot with no tickmarks, where the axis ticks and
## labels are closer to the axes than the default settings
trimplot <- function(x, xlab = "", ylab = "", ...) {
    plot(x, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
    axis(side = 1, tick = FALSE, line = -0.75)
    axis(side = 2, tick = FALSE, line = -0.75)
    mtext(side = 1, text = xlab, line =  1.50)
    mtext(side = 2, text = ylab, line =  1.50)
    invisible(NULL)
}

## The code for this function was more or less taken from Stack Overflow:
## https://stackoverflow.com/a/39802705
## (the only innovations being putting it in a function and making the
## number of digits customizable)
## It rounds numbers toward +/- infinity (i.e. away from zero, or "out")
round_out <- function(x, digits = 2) {
    z <- 10^digits
    return(sign(x) * ceiling(abs(x) * z) / z)
}


#' Determine whether items' response functions are non-monotonic on a support
#' 
#' @param x An object of class \code{\link[bggum]{summary.ggum"}}
#' @param from A numeric vector of length one giving the lower bound of the
#'   support; if \code{NULL} (the default), the lowest theta value is used
#' @param to A numeric vector of length one giving the upper bound of the
#'   support; if \code{NULL} (the default), the highest theta value is used
#' @param n_points An integer vector of length one giving the number of points
#'   in the grid searched for non-monotonicity; the default is 100
#'
#' @return A logical vector of the same length as \code{x$estimates$alpha};
#'   each element is \code{TRUE} if the corresponding item is non-monotonic
#'   on the support and \code{FALSE} if the corresponding item is monotonic
is_nonmonotonic <- function(x, from = NULL, to = NULL, n_points = 100) {
    ## Sanity checks
    if ( !("summary.ggum" %in% class(x)) ) {
        stop("x must be a summary.ggum object")
    }
    if ( is.null(from) ) {
        from <- min(x$estimates$theta)
    }
    if ( is.null(to) ) {
        to <- max(x$estimates$theta)
    }
    if ( from > to ) {
        stop("to must be greater than from")
    }
    if ( n_points < 1 ) {
        stop("n_points must be a positive integer")
    }
    ## Extract item parameters for convenience
    alpha <- x$estimates$alpha
    delta <- x$estimates$delta
    tau   <- x$estimates$tau
    ## Create vector of theta values; a grid on the support
    theta <- seq(from = from, to = to, length.out = n_points)
    ## Allocate storage for result vector
    m   <- length(alpha)
    res <- logical(length = m)
    ## Get the number of options for each item
    if ( is.list(tau) ) {
        K <- sapply(tau, length)
    }
    else {
        K   <- length(tau)
        tau <- list(tau)
    }
    ## Check each item for non-monotonicity
    for ( j in 1:m ) {
        ## Get the expected response at each value of theta
        kk <- K[j] - 1
        p <- sapply(theta, function(x) {
            sum(0:kk * ggumProbability(0:kk, x, alpha[j], delta[j], tau[[j]]))
        })
        ## Check that all adjacent expected responses are either weakly larger
        ## or weakly smaller; we round to ensure ~very small~ fluctuations
        ## don't give us a false positive
        y <- round(p[-1], 3)
        z <- round(p[-length(p)], 3)
        res[j] <- !(all(y <= z) | all(y >= z))
    }
    ## Return result
    return(res)
}

##
prop_nonmonotonic <- function(x, from = NULL, to = NULL, n_points = 100) {
    nonmonotonic <- is_nonmonotonic(x, from, to, n_points)
    return(sum(nonmonotonic) / length(nonmonotonic))
}


#' Get probability of responses given parameter estimates for the CJR model
#' 
#' @return A numeric vector or matrix of response probabilities
cjrProbability <- function(response, theta, alpha, beta) {
    ## Sanity checks
    if ( !is.matrix(theta) ) {
        theta <- matrix(theta, ncol = 1)
    }
    if ( !is.matrix(beta) ) {
        beta <- matrix(beta, ncol = 1)
        # alpha <- -alpha
    }
    if ( length(alpha) != nrow(beta) ) {
        stop("The lengths of alpha and beta should be the same")
    }
    n <- nrow(theta)
    m <- length(alpha)
    if ( is.matrix(response) ) {
        if ( nrow(response) != n ) {
            stop("The length or # of rows of theta should equal nrow(response)")
        }
        if ( ncol(response) != m ) {
            stop("The length of alpha should equal ncol(response)")
        }
    } else {
        if ( n > 1 & m > 1 ) {
            stop("With multiple items & respondents, response must be a matrix")
        } else if ( n == 1 ) {
            if ( length(response) != m ) {
                stop("length(response) should equal length(alpha)")
            }
        } else if ( m == 1 ) {
            if ( length(response) != n ) {
                stop("length(response) should equal length or nrow(theta)")
            }
        }
        response <- matrix(response, nrow = n, ncol = m)
    }
    ## Allocate result storage
    res <- matrix(NA_real_, nrow = n, ncol = m)
    ## Calculate probabilities
    z <- t(apply(beta %*% t(theta), 2, "-", alpha))
    for ( j in 1:m ) {
        for ( i in 1:n ) {
            res[i, j] <- pnorm(z[i, j], lower.tail = response[i, j])
        }
    }
    ## Return result
    return(res)
}


## This is a slightly modified version of wnominate::plot.coords()
## that allows us to color the parties how they should be
plot_coords <- function(
    x,
    main.title = "",
    d1.title = "First Dimension", d2.title = "Second Dimension",
    dims = c(1, 2), plotBy = "party", 
    color = TRUE, shape = TRUE, cutline = NULL, Legend = TRUE, 
    legend.x = 0.75, legend.y = 1, ...
) {
    if (!class(x) == "nomObject") 
        stop("Input is not of class 'nomObject'.")
    if ( length(party) > 1 ) {
        types <- party
    }
    else if (!any(colnames(x$legislators) == plotBy)) {
        warning("Variable '", plotBy, "' does not exist in your W-NOMINATE object.")
        types <- rep("Leg", dim(x$legislators)[1])
    }
    else {
        types <- x$legislators[, plotBy]
    }
    if (length(dims) != 2 & x$dimensions != 1) 
        stop("'dims' must be an integer vector of length 2.")
    nparties <- length(unique(types))
    colorlist <- c("firebrick3", "dodgerblue4", "darkorchid")
    shapes <- c(16, 15, 17)
    if (color == FALSE) 
        colorlist <- sample(colors()[160:220], 50)
    if (shape == FALSE) 
        shapes <- rep(16, 50)
    if (x$dimensions == 1) {
        coord1D <- x$legislators[, "coord1D"]
        ranking <- rank(x$legislators[, "coord1D"])
        plot(seq(-1, 1, length = length(coord1D)), 1:length(coord1D), 
             type = "n", cex.main = 1.2, cex.lab = 1.2, font.main = 2, 
             xlab = "First Dimension Nominate", ylab = "Rank", 
             main = "1D W-NOMINATE Plot")
        if (Legend) 
            legend(0.67, 0.7 * length(coord1D), unique(types), 
                   pch = shapes[1:nparties], col = colorlist[1:nparties], 
                   cex = 0.7)
        for (i in 1:nparties) {
            suppressWarnings(
                points(
                    coord1D[types == unique(types)[i]],
                    ranking[types == unique(types)[i]], 
                    pch = shapes[i], col = colorlist[i], cex = 1.1, lwd = 2
                )
            )
        }
    }
    else {
        coord1D <- x$legislators[, paste("coord", dims[1], "D", sep = "")]
        coord2D <- x$legislators[, paste("coord", dims[2], "D", sep = "")]
        suppressWarnings(
            symbols(
                x = 0, y = 0, circles = 1, inches = FALSE, 
                asp = 1, main = main.title, xlab = d1.title, ylab = d2.title, 
                xlim = c(-1, 1), ylim = c(-1, 1), lwd = 2, fg = "grey", ...
            )
        )
        if (Legend) 
            legend(legend.x, legend.y, unique(types), pch = shapes[1:nparties], 
                   bty = "n", col = colorlist[1:nparties])
        for (i in 1:nparties) {
            suppressWarnings(
                points(
                    coord1D[types == unique(types)[i]],
                    coord2D[types == unique(types)[i]], 
                    pch = shapes[i], col = colorlist[i], lwd = 2
                )
            )
        }
    }
}

##### Setup -----
## Load packages
library(bggum)
## Ensure required directories exist
source("code/00-util.R")
prepare_directories()
## While I have a response probability function for the GGUM in our package,
## I need to set up functions for a monotonic IRT and the GPCM for these plots.
## Here's the IRT model from Clinton Jackman Rivers:
cjr <- function(theta, discrimination, location, response) {
    return(pnorm(location + theta*discrimination, lower.tail = response))
}
## Here's the General Partial Credit Model:
gpcm <- function(theta, discrimination, location, thresholds, response) {
    M = length(thresholds) - 1
    numerator <- 0
    denominator <- 0
    tausum <- cumsum(thresholds)
    for ( i in 0:M ) {
        tmp <- exp(discrimination * (i * (theta - location) - tausum[i+1]))
        if ( i == response ) {
            numerator <- tmp
        }
        denominator <- denominator + tmp
    }
    return(numerator / denominator)
}


##### Reproduce Figure 1 -----
## Figure 1 has two panels:
## (a) Example IRF for the standard 2PL IRT model
## (b) Example IRF for GPCM

##### Figure 1(a)
## Set example parameter values:
discrimination <- 1
location <- 0
theta <- seq(-3, 3, 0.01)
## Calculate probability of 1 and 0 responses:
p1 <- sapply(theta, cjr, discrimination, location, 1)
p0 <- sapply(theta, cjr, discrimination, location, 0)
## Setup the plotting device:
tiff(
    filename = "plots/fig1a.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1, xpd = TRUE)
## Plot the zero response probabilities (and set up the plotting area):
plot(
    x = theta, y = p0,
    type = "l", lty = 2, lwd = 2,
    ylim = c(0, 1.01), ylab = "", yaxt = "n", xaxt = "n", xlab = "", main = ""
)
## Draw a custom x axis:
axis(side = 1, at = -3:3, tick = FALSE, line = -0.75)
mtext(expression(theta), side = 1, line = 1.5)
## Draw a custom y axis:
axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), tick = FALSE, line = -0.75)
mtext("Probability", side = 2, line = 1.5)
## Plot the one response probabilities:
lines(theta, p1, lwd = 2)
## Add a legend
legend(
    "top", inset = c(0, -0.15), bty = "n", horiz = TRUE,
    lty = 2:1, lwd = 2, legend = c("Disagree", "Agree")
)
## Reset the margins
par(opar)
## And close the plotting device
dev.off()

##### Figure 1(b)
## Set example parameter values:
discrimination <- 1
location <- 0
thresholds <- c(0, -1, 0, 1)
theta <- seq(-3, 3, 0.01)
## Calculate probability of responses:
p3 <- sapply(theta, gpcm, discrimination, location, thresholds, 3)
p2 <- sapply(theta, gpcm, discrimination, location, thresholds, 2)
p1 <- sapply(theta, gpcm, discrimination, location, thresholds, 1)
p0 <- sapply(theta, gpcm, discrimination, location, thresholds, 0)
## Setup the plotting device:
tiff(
    filename = "plots/fig1b.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1, xpd = TRUE)
## Plot the zero response probabilities (and set up the plotting area):
plot(
    x = theta, y = p0,
    type = "l", lty = 4, lwd = 2,
    ylim = c(0, 1.01), ylab = "", yaxt = "n", xaxt = "n", xlab = "", main = ""
)
## Draw a custom x axis:
axis(side = 1, at = -3:3, tick = FALSE, line = -0.75)
mtext(expression(theta), side = 1, line = 1.5)
## Draw a custom y axis:
axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), tick = FALSE, line = -0.75)
mtext("Probability", side = 2, line = 1.5)
## Plot the other response probabilities:
lines(theta, p1, lwd = 2, lty = 3)
lines(theta, p2, lwd = 2, lty = 2)
lines(theta, p3, lwd = 2, lty = 1)
## Add a legend
legend(
    "top", inset = c(0, -0.2), bty = "n", lty = c(4, 1, 3:2), lwd = 2, ncol = 2,
    legend = c(
        "Disagree from below", "Disagree from above",
        "Agree from below", "Agree from above"
    )
)
## Reset the margins
par(opar)
## And close the plotting device
dev.off()


##### Reproduce Figure 2 -----
## Figure 2 has two panels:
## (a) Example IRF for GGUM, plotted near bliss point
## (b) Example IRF for GGUM, plotted far from bliss point

##### Figure 2(a)
## Set example parameter values:
alpha <- 1
delta <- 0
tau   <- c(0, -1)
theta <- seq(-3, 3, 0.01)
## Calculate probability of responses:
n  <- length(theta)
p1 <- ggumProbability(rep(1, n), theta, alpha, delta, tau)
p0 <- ggumProbability(rep(0, n), theta, alpha, delta, tau)
## Setup the plotting device:
tiff(
    filename = "plots/fig2a.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1, xpd = TRUE)
## Plot the zero response probabilities (and set up the plotting area):
plot(
    x = theta, y = p0,
    type = "l", lty = 2, lwd = 2,
    ylim = c(0, 1.01), ylab = "", yaxt = "n", xaxt = "n", xlab = "", main = ""
)
## Draw a custom x axis:
axis(side = 1, at = -3:3, tick = FALSE, line = -0.75)
mtext(expression(theta - delta), side = 1, line = 1.5)
## Draw a custom y axis:
axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), tick = FALSE, line = -0.75)
mtext("Probability", side = 2, line = 1.5)
## Plot the one response probabilities:
lines(theta, p1, lwd = 2)
## Add a legend
legend(
    "top", inset = c(0, -0.15), bty = "n", horiz = TRUE,
    lty = 2:1, lwd = 2, legend = c("Disagree", "Agree")
)
## Reset the margins
par(opar)
## And close the plotting device
dev.off()

##### Figure 2(b)
## Set example parameter values:
alpha <- 1
delta <- 3
tau   <- c(0, -3)
theta <- seq(-3, 3, 0.01)
## Calculate probability of responses:
n  <- length(theta)
p1 <- ggumProbability(rep(1, n), theta, alpha, delta, tau)
p0 <- ggumProbability(rep(0, n), theta, alpha, delta, tau)
## Setup the plotting device:
tiff(
    filename = "plots/fig2b.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1, xpd = TRUE)
## Plot the zero response probabilities (and set up the plotting area):
plot(
    x = theta, y = p0,
    type = "l", lty = 2, lwd = 2,
    ylim = c(0, 1.01), ylab = "", yaxt = "n", xaxt = "n", xlab = "", main = ""
)
## Draw a custom x axis:
axis(side = 1, at = -3:3, tick = FALSE, line = -0.75, labels = -6:0)
mtext(expression(theta - delta), side = 1, line = 1.5)
## Draw a custom y axis:
axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), tick = FALSE, line = -0.75)
mtext("Probability", side = 2, line = 1.5)
## Plot the one response probabilities:
lines(theta, p1, lwd = 2)
## Add a legend
legend(
    "top", inset = c(0, -0.15), bty = "n", horiz = TRUE,
    lty = 2:1, lwd = 2, legend = c("Disagree", "Agree")
)
## Reset the margins
par(opar)
## And close the plotting device
dev.off()


##### Reproduce Figure 3 -----
## Figure 3 
## Set example parameter values:
discrimination <- 1
location <- 0
thresholds <- c(0, -1, 0, 1)
theta <- seq(-3, 3, 0.01)
## Calculate probability of responses:
p3 <- sapply(theta, gpcm, discrimination, location, thresholds, 3)
p2 <- sapply(theta, gpcm, discrimination, location, thresholds, 2)
p1 <- sapply(theta, gpcm, discrimination, location, thresholds, 1)
p0 <- sapply(theta, gpcm, discrimination, location, thresholds, 0)
## Setup the plotting device:
tiff(
    filename = "plots/fig3.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1)
## Plot the zero response probabilities (and set up the plotting area):
plot(
    x = theta, y = p0,
    type = "l", lty = 4, lwd = 2,
    ylim = c(0, 1.01), ylab = "", yaxt = "n", xaxt = "n", xlab = "", main = ""
)
## Draw a custom x axis:
axis(side = 1, at = -3:3, tick = FALSE, line = -0.75)
axis(
    side = 1, at = c(-1, 1), tick = FALSE, line = 0.5,
    labels = expression(tau[1], -tau[1])
)
mtext(expression(theta - delta), side = 1, line = 1.5)
## Draw a custom y axis:
axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), tick = FALSE, line = -0.75)
mtext("Probability", side = 2, line = 1.5)
## Plot the other response probabilities:
lines(theta, p1, lwd = 2, lty = 3)
lines(theta, p2, lwd = 2, lty = 2)
lines(theta, p3, lwd = 2, lty = 1)
## Add vertical lines at tau_1, 0, -tau_1
abline(v = -1, col = "gray")
abline(v =  0, col = "gray")
abline(v =  1, col = "gray")
## Add a legend
par(xpd = TRUE)
legend(
    "top", inset = c(0, -0.2), bty = "n", lty = c(4, 1, 3:2), lwd = 2, ncol = 2,
    legend = c(
        "Disagree from below", "Disagree from above",
        "Agree from below", "Agree from above"
    )
)
## Reset the margins
par(opar)
## And close the plotting device
dev.off()


##### Reproduce Figures A.1-3 -----
## Appendix A has three figures
## A.1. Demonstrates the effect of changing the alpha parameter
## A.2. Demonstrates the effect of changing the delta parameter
## A.3. Demonstrates the effect of changing the tau parameter

##### Figure A.1
## delta and tau don't change for this one; set them now
delta <- 0
tau   <- c(0, -1)
## panel (a) has alpha = 0.5
alpha <- 0.5
pdf("plots/figA1a.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()
## panel (b) has alpha = 1.0
alpha <- 1.0
pdf("plots/figA1b.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()
## panel (c) has alpha = 2.0
alpha <- 2.0
pdf("plots/figA1c.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()

##### Figure A.2
## alpha and tau don't change for this one; set them now
alpha <- 1
tau   <- c(0, -1)
## panel (a) has delta = -1
delta <- -1
pdf("plots/figA2a.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()
## panel (b) has delta =  0
delta <-  0
pdf("plots/figA2b.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()
## panel (c) has delta =  1
delta <-  1
pdf("plots/figA2c.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()

##### Figure A.3
## alpha and delta don't change for this one; set them now
alpha <- 1
delta <- 0
## panel (a) has tau = (0, -0.5)
tau   <- c(0, -0.5)
pdf("plots/figA3a.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()
## panel (b) has tau = (0, -1.0)
tau   <- c(0, -1.0)
pdf("plots/figA3b.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()
## panel (c) has tau = (0, -2.0)
tau   <- c(0, -2.0)
pdf("plots/figA3c.pdf", width = 6, height = 4)
irf(alpha, delta, tau, main_title = "")
dev.off()


##### Reproduce Figure B.1 -----
## For Figure B.1, we simulate data from 500 respondents to 10 4-option items
set.seed(42)                            ## For reproducibility
sim_data <- ggum_simulation(500, 10, 4) ## From the bggum package
## Extract item parameters for convenience
alpha <- sim_data$alpha
delta <- sim_data$delta
tau <- sim_data$tau
## Now we'll plot profile likelihoods for some bimodal theta parameters
pdf("plots/figB1.pdf", height = 4, width = 12)
## Set layout and plotting parameters
layout(matrix(1:2, nrow = 1))
opar <- par(mar = c(3, 2, 1, 1) + 0.1)
theta_range <- seq(-3.5, 3.5, 0.01)
for ( i in c(241, 445) ) {
    ## Calculate the profile likelihood
    responses <- sim_data$response_matrix[i, ]
    likelihood <- exp(sapply(theta_range, function(th) {
        sum(log(ggumProbability(responses, th, alpha, delta, tau)))
    }))
    ## Plot profile likelihood
    plot(
        x = theta_range, y = likelihood, type = "l",
        main = "", xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    ## Add indicator of true value
    abline(v = sim_data$theta[i], lty = 2, col = "#808080")
    ## Add and label axes
    axis(1, at = -3:3, tick = FALSE, line = -0.75)
    mtext(side = 1, line = 1.50, text = expression("Ideology (" * theta * ")"))
    mtext(side = 2, line = 0.25, text = expression(L(X[i])))
}
par(opar)
dev.off()

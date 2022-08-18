##### Setup -----
## For this demonstration, we need functionality not in the standard release
## of the `bggum` package: We need to be able to select which chain's draws
## we're going to record. The CRAN version of the package does not allow that
## because inference is only valid from the cold chain. However, we can install
## an alternative version from a specific GitHub branch opened to allow it.
## We'll install this in a temporary directory both so as not to disturb our
## regular `bggum` installation and because there will be no reason to hold
## on to the installation of this version of the package
install_dir <- tempdir()
devtools::install_github(
    repo = "duckmayr/bggum",
    ref  = "choose-recorded-chain",
    lib  = install_dir
)
## Now we make sure the official version is not loaded,
## then load the choose-recorded-chain version
if ( "bggum" %in% loadedNamespaces() ) {
    detach(name = "package:bggum", unload = TRUE, character.only = TRUE)
}
library(bggum, lib.loc = install_dir)
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Ensure required directories exist
prepare_directories()


##### Data prep -----
## Deal with reproducibility
set.seed(123)

## Simulate data
n  <- 100 # number of respondents
m  <- 10  # number of items
K  <- 4   # number of options per item
ap <- c(1.5, 1.5,  0.5, 3.0) # alpha (discrimination parameter) distribution
dp <- c(2.0, 2.0, -3.0, 3.0) # delta (location parameter) distribution
sim_data <- ggum_simulation(n, m, K, alpha_params = ap, delta_params = dp)
response_matrix <- sim_data$response_matrix

## Check unanimity issues
# apply(response_matrix, 2, function(x) length(unique(x)))
# # [1] 4 4 4 4 4 4 4 4 4 4 # good


##### Generate posterior samples -----
## Set hyperparameters
sds <- tune_proposals(response_matrix, tune_iterations = 5000)
sapply(sds, mean) # sanity check
# [1] 1.61010 1.10800 0.66000 1.14625

## Record 1000 iterations of cold chain
set.seed(42)
cold_samples <- ggumMC3(
    data = response_matrix,
    sample_iterations = 1000,
    burn_iterations = 0,
    proposal_sds = sds,
    temps = c(1, 0.2),
    swap_interval = 10000
)

## Record 1000 iterations of hottest chain
set.seed(42)
hot_samples <- ggumMC3(
    data = response_matrix,
    sample_iterations = 1000,
    burn_iterations = 0,
    proposal_sds = sds,
    temps = c(1, 0.2),
    swap_interval = 10000,
    recorded_chain = 2
)


##### Reproduce Figure 4 -----
## Setup the plotting device for panel (a): Trace plot of cold chain
tiff(
    filename = "plots/fig4a.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(4, 3, 1, 1) + 0.1, xpd = FALSE)
## Set the palette
pal <- c(cold = "black", hot = "black", truth = "#80808080")
## Plot theta1 trace
i <- 1
trimplot(
    x = NULL, type = "l", ylim = c(-4, 4), xlim = c(1, 1000),
    xlab = "Iteration", ylab = bquote(theta[1]~"draw")
)
abline(h = sim_data$theta[i], lty = 2, lwd = 2, col = pal["truth"])
lines(cold_samples[1:1000, i], col = pal["cold"], lwd = 2)
## Reset the margins
par(opar)
## And close the plotting device
invisible(dev.off())

## Setup the plotting device for panel (b): Trace plot of hot chain
tiff(
    filename = "plots/fig4b.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(4, 3, 1, 1) + 0.1, xpd = FALSE)
## Set the palette
pal <- c(cold = "black", hot = "black", truth = "#80808080")
## Plot theta1 trace
i <- 1
trimplot(
    x = NULL, type = "l", ylim = c(-4, 4), xlim = c(1, 1000),
    xlab = "Iteration", ylab = bquote(theta[1]~"draw")
)
abline(h = sim_data$theta[i], lty = 2, lwd = 2, col = pal["truth"])
lines(hot_samples[1:1000, i],  col = pal["hot"], lwd = 2)
## Reset the margins
par(opar)
## And close the plotting device
invisible(dev.off())

## Setup the plotting device for panel (c): Density plots together
tiff(
    filename = "plots/fig4c.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(4, 3, 1, 1) + 0.1, xpd = FALSE)
## Set the palette
pal <- c(cold = "black", hot = "black", truth = "#80808080")
## Plot theta1 trace
i <- 1
cold_density <- density(cold_samples[251:1000, i], from = -4, to = 4)
hot_density  <- density(hot_samples[251:1000, i],  from = -4, to = 4)
trimplot(
    x = NULL, type = "l", xlim = c(-4, 4), ylim = c(0, 1.31),
    ylab = "Density", xlab = bquote(theta[1])
)
abline(v = sim_data$theta[i], lty = 3, lwd = 3, col = pal["truth"])
lines(cold_density, lwd = 2, col = pal["cold"], lty = 2)
lines(hot_density,  lwd = 2, col = pal["hot"])
## Add a legend
legend(
    "topright", bty = "n", col = pal, lty = c(2, 1, 3), lwd = 2,
    legend = c(
        expression(beta == 1),
        expression(beta == 0.2),
        expression(theta[1])
    )
)
## Reset the margins
par(opar)
## And close the plotting device
invisible(dev.off())

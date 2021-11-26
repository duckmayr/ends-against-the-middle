##### Setup -----
## Load packages
library(bggum)
library(coda)
library(MCMCpack)
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Read in cleaned vote data
votes <- read.csv("data/v23.csv")
## Read in original dataset that has info on the proposals
dat <- read.csv("data/base_ife_eric_feb2021.csv")
## (These data were obtained from the Estevez et al (2008) authors' GitHub:
## https://github.com/emagar/ife-update/tree/master/data)


##### Notes -----
## An advantage of GGUM to the model used in Estevez et al (2008)
## is we can treat abstention as its own category,
## and we can accommodate ends against the middle behavior
## (which we will see a non-trivial amount of in the data).
## As in Estevez et al (2008), we'll analyze the terms separately.


#### Sampling for Term 2 (Woldenberg I) -----
## Reshape into a response matrix
t2r <- which(votes$term == 2)             ## Use only rows for term 2
responses <- as.matrix(votes[t2r, 1:11])  ## First 11 columns are respondents
rownames(responses) <- votes$folio[t2r]   ## Add vote IDs in dimnames
responses <- t(responses)                 ## Reshape so respondents are rows
## Recode votes
responses[responses %in% c(0, 4:6)] <- NA ## Recode all non-votes to NA
responses[responses == 2] <- 0            ## Recode all nays to 0
responses[responses == 1] <- 2            ## Recode all yeas to 2
responses[responses == 3] <- 1            ## Recode all abstention to 1
## Eliminate respondents with no votes (i.e., those only present in term 3)
has_no_votes <- apply(responses, 1, function(r) all(is.na(r)))
responses <- responses[!has_no_votes, ]
## Eliminate unanimous votes
is_unanimous <- apply(responses, 2, function(v) {
    return(length(unique(na.omit(v))) == 1 | all(is.na(v)))
})
t2responses <- responses[ , !is_unanimous]
## Ensure votes are always in {NA, 0, ..., K-1}
m <- ncol(t2responses)
types <- vector(mode = "list", length = m)
set_equal <- function(x, y) setequal(c(NA, x), c(NA, y))
for ( j in 1:m ) {
    if ( set_equal(c(NA, 0, 1), unique(t2responses[ , j])) ) {
        types[[j]] <- c("nay" = 0, "abstain" = 1)
    } else if ( set_equal(c(NA, 0, 2), unique(t2responses[ , j])) ) {
        types[[j]] <- c("nay" = 0, "yea" = 1)
        t2responses[t2responses[ , j] == 2, j] <- 1
    } else if ( set_equal(c(NA, 1, 2), unique(t2responses[ , j])) ) {
        types[[j]] <- c("abstain" = 0, "yea" = 1)
        t2responses[ , j] <- t2responses[ , j] - 1
    } else if ( set_equal(c(NA, 0, 1, 2), unique(t2responses[ , j])) ) {
        types[[j]] <- c("nay" = 0, "abstain" = 1, "yea" = 2)
    } else {
        warning(paste("Check", j))
    }
}
## Double-check
for ( j in 1:m ) {
    k <- max(t2responses[ , j], na.rm = TRUE)
    if ( !set_equal(c(NA, 0:k), t2responses[ , j]) ) {
        warning(paste("Check", j))
    }
}
## Set seed for reproducibility
set.seed(42)
## Tune proposals
sds <- tune_proposals(t2responses, 5000)
## Tune temperatures
temps <- tune_temperatures(t2responses, 6, proposal_sds = sds)
## Sample posterior
samples <- lapply(1:2, function(x) {
    ggumMC3(
        data = t2responses,
        sample_iterations = 50000,
        burn_iterations = 5000,
        proposal_sds = sds,
        temps = temps
    )
})
## Save results to disk
saveRDS(object = samples, file = "output/ife-term2-samples.rds")
## Post-process
for ( i in 1:nrow(responses) ) {
    plot(density(samples[[1]][ , i]), main = paste(rownames(t2responses)[i], i))
} ## 4 should be constraint; Cardenas should be negative per EMR 2008
samples <- lapply(samples, post_process, constraint = 4, expected_sign = "-")
## Diagnostics
conv_stats <- gelman.diag(samples)
summary(conv_stats$psrf[ , 1])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.000   1.001   1.009 
## Analysis
term2estimates <- summary(samples)
saveRDS(term2estimates, file = "output/ife-term2-estimates.rds")
t2theta <- term2estimates$estimates$theta
t2alpha <- term2estimates$estimates$alpha
t2delta <- term2estimates$estimates$delta
t2tau   <- term2estimates$estimates$tau
t2theta_sd <- term2estimates$sds$theta_sds


#### Sampling for Term 3 (Woldenberg II) -----
## Reshape into a response matrix
t3r <- which(votes$term == 3)             ## Use only rows for term 2
responses <- as.matrix(votes[t3r, 1:11])  ## First 11 columns are respondents
rownames(responses) <- votes$folio[t3r]   ## Add vote IDs in dimnames
responses <- t(responses)                 ## Reshape so respondents are rows
## Recode votes
responses[responses %in% c(0, 4:6)] <- NA ## Recode all non-votes to NA
responses[responses == 2] <- 0            ## Recode all nays to 0
responses[responses == 1] <- 2            ## Recode all yeas to 2
responses[responses == 3] <- 1            ## Recode all abstention to 1
## Eliminate respondents with no votes (i.e., those only present in term 2)
has_no_votes <- apply(responses, 1, function(r) all(is.na(r)))
responses <- responses[!has_no_votes, ]
## Eliminate unanimous votes
is_unanimous <- apply(responses, 2, function(v) {
    return(length(unique(na.omit(v))) == 1 | all(is.na(v)))
})
t3responses <- responses[ , !is_unanimous]
## Ensure votes are always in {NA, 0, ..., K-1}
m <- ncol(t3responses)
types <- vector(mode = "list", length = m)
set_equal <- function(x, y) setequal(c(NA, x), c(NA, y))
for ( j in 1:m ) {
    if ( set_equal(c(NA, 0, 1), unique(t3responses[ , j])) ) {
        types[[j]] <- c("nay" = 0, "abstain" = 1)
    } else if ( set_equal(c(NA, 0, 2), unique(t3responses[ , j])) ) {
        types[[j]] <- c("nay" = 0, "yea" = 1)
        t3responses[t3responses[ , j] == 2, j] <- 1
    } else if ( set_equal(c(NA, 1, 2), unique(t3responses[ , j])) ) {
        types[[j]] <- c("abstain" = 0, "yea" = 1)
        t3responses[ , j] <- t3responses[ , j] - 1
    } else if ( set_equal(c(NA, 0, 1, 2), unique(t3responses[ , j])) ) {
        types[[j]] <- c("nay" = 0, "abstain" = 1, "yea" = 2)
    } else {
        warning(paste("Check", j))
    }
}
## Double-check
for ( j in 1:m ) {
    k <- max(t3responses[ , j], na.rm = TRUE)
    if ( !set_equal(c(NA, 0:k), t3responses[ , j]) ) {
        warning(paste("Check", j))
    }
}
## Set seed for reproducibility
set.seed(138)
## Tune proposals
sds <- tune_proposals(t3responses, 5000)
## Tune temperatures
temps <- tune_temperatures(t3responses, 6, proposal_sds = sds)
temps
## Sample posterior
samples <- lapply(1:2, function(x) {
    ggumMC3(
        data = t3responses,
        sample_iterations = 50000,
        burn_iterations = 5000,
        proposal_sds = sds,
        temps = temps
    )
})
## Save results to disk
saveRDS(object = samples, file = "output/ife-term3-samples.rds")
## Post-process
for ( i in 1:nrow(t3responses) ) {
    plot(density(samples[[1]][ , i]), main = paste(rownames(t3responses)[i], i))
} ## 4 should be constraint; Cardenas should be negative per EMR 2008
samples <- lapply(samples, post_process, constraint = 4, expected_sign = "-")
## Diagnostics
conv_stats <- gelman.diag(samples)
summary(conv_stats$psrf[ , 1])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.000   1.001   1.018 
## Analysis
term3estimates <- summary(samples)
saveRDS(term3estimates, file = "output/ife-term3-estimates.rds")
t3theta <- term3estimates$estimates$theta
t3alpha <- term3estimates$estimates$alpha
t3delta <- term3estimates$estimates$delta
t3tau   <- term3estimates$estimates$tau
t3theta_sd <- term3estimates$sds$theta_sds


##### Reproduce Figure K1 -----
## Generate some variables used in both panels
ths <- seq(
    from = round_out(min(t2theta)),
    to = round_out(max(t2theta)),
    by = 0.01
)
N <- length(ths)
xlabel <- expression(theta)
ylabel <- "Pr( Yea )"

## Reproduce panel (a) (folio 871)
j   <- which(colnames(t2responses) == "871")
pal <- c(okabe_ito()[c(6, 5)], "#808080")
pal <- pal[c(1, 3, 2)]
pdf("plots/figK1a.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 2, 1) + 0.1, xpd = TRUE)
## Setup plot and add line for Yea
p <- ggumProbability(rep(2, N), ths, t2alpha[j], t2delta[j], t2tau[[j]])
trimplot(
    ths, y = p, type = "l",
    col = pal[3], lwd = 2,
    ylim = c(0, 1), xlab = xlabel, ylab = ylabel
)
## Add line for Abstain
p <- ggumProbability(rep(1, N), ths, t2alpha[j], t2delta[j], t2tau[[j]])
lines(ths, p, lwd = 2, col = pal[2], lty = 3)
## Add line for Nay
p <- ggumProbability(rep(0, N), ths, t2alpha[j], t2delta[j], t2tau[[j]])
lines(ths, p, lwd = 2, col = pal[1], lty = 2)
p <- ggumProbability(t2responses[ , j], t2theta, t2alpha[j], t2delta[j], t2tau[[j]])
zeros <- which(t2responses[ , j] == 0)
ones  <- which(t2responses[ , j] == 1)
twos  <- which(t2responses[ , j] == 2)
points(t2theta[zeros], p[zeros], pch = 19, col = paste0(pal[1], "80"), cex = 1.5)
points(t2theta[ones],  p[ones],  pch = 19, col = paste0(pal[2], "80"), cex = 1.5)
points(t2theta[twos],  p[twos],  pch = 19, col = paste0(pal[3], "80"), cex = 1.5)
legend("top", inset = c(0, -0.2), lty = c(2, 3, 1), bty = "n", horiz = TRUE,
       legend = c("Nay", "Abstain", "Yea"), col = pal, lwd = 2)
par(opar)
dev.off()

## Reproduce panel (b) (folio 878)
j   <- which(colnames(t2responses) == "878")
pal <- okabe_ito()[c(6, 5)]
pdf("plots/figK1b.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 2, 1) + 0.1, xpd = TRUE)
p <- ggumProbability(rep(1, N), ths, t2alpha[j], t2delta[j], t2tau[[j]])
trimplot(
    ths, y = p, type = "l",
    col = pal[2], lwd = 2,
    ylim = c(0, 1), xlab = xlabel, ylab = ylabel
)
p <- ggumProbability(rep(0, N), ths, t2alpha[j], t2delta[j], t2tau[[j]])
lines(ths, p, lwd = 2, col = pal[1], lty = 2)
p <- ggumProbability(t2responses[ , j], t2theta, t2alpha[j], t2delta[j], t2tau[[j]])
zeros <- which(t2responses[ , j] == 0)
ones  <- which(t2responses[ , j] == 1)
points(t2theta[ones],  p[ones],  pch = 19, col = paste0(pal[2], "80"), cex = 1.5)
points(t2theta[zeros], p[zeros], pch = 19, col = paste0(pal[1], "80"), cex = 1.5)
legend("top", inset = c(0, -0.2), lty = 2:1, bty = "n", horiz = TRUE,
       legend = c("Nay", "Yea"), col = pal, lwd = 2)
par(opar)
dev.off()


##### Reproduce Table K1 -----
## First add the info from Table 6 in EMR (2008)
sponsors <- c(
    "cardenas"   = "PRD",
    "cantu"      = "PT",
    "zebadua"    = "PRD",
    "lujambio"   = "PAN",
    "molinar"    = "PAN",
    "merino"     = "PRI",
    "woldenberg" = "PRI",
    "peschard"   = "PRI",
    "barragan"   = "PRD",
    "luken"      = "PAN",
    "rivera"     = "PRI"
)
term2EMRmean <- c(
    "cardenas"   = -1.79,
    "cantu"      =  0.42,
    "zebadua"    =  0.73,
    "lujambio"   =  0.90,
    "molinar"    =  1.09,
    "merino"     =  1.95,
    "woldenberg" =  2.15,
    "peschard"   =  2.28,
    "barragan"   =  3.25
)
term3EMRmean <- c(
    "cardenas"   = -1.67,
    "barragan"   =  0.40,
    "cantu"      =  1.70,
    "luken"      =  1.98,
    "rivera"     =  3.20,
    "lujambio"   =  3.50,
    "merino"     =  3.60,
    "woldenberg" =  3.70,
    "peschard"   =  3.75
)
term2EMRsd <- c(
    "cardenas"   = 0.44,
    "cantu"      = 0.20,
    "zebadua"    = 0.21,
    "lujambio"   = 0.25,
    "molinar"    = 0.26,
    "merino"     = 0.45,
    "woldenberg" = 0.53,
    "peschard"   = 0.60,
    "barragan"   = 1.03
)
term3EMRsd <- c(
    "cardenas"   = 0.23,
    "barragan"   = 0.12,
    "cantu"      = 0.20,
    "luken"      = 0.24,
    "rivera"     = 0.38,
    "lujambio"   = 0.45,
    "merino"     = 0.44,
    "woldenberg" = 0.47,
    "peschard"   = 0.44
)

## Now put it all together for the comparison table
t2order    <- order(term2theta)
t3order    <- order(term3theta)
t2members  <- rownames(t2responses)[t2order]
t3members  <- rownames(t3responses)[t3order]
comparison <- data.frame(
    Term        = c(rep("Woldenberg I", 9), rep("Woldenberg II", 9)),
    Councilor   = c(t2members, t3members),
    Sponsor     = c(sponsors[t2members], sponsors[t3members]),
    `EMR Mean`  = c(term2EMRmean[t2members], term3EMRmean[t3members]),
    `EMR SD`    = c(term2EMRsd[t2members],   term3EMRsd[t3members]),
    `GGUM Mean` = sprintf("%0.2f", c(t2theta[t2order], t3theta[t3order])),
    `GGUM SD`   = sprintf("%0.2f", c(t2theta_sd[t2order], t3theta_sd[t3order]))
)
comparison

##### Setup -----
## Load required packages
library(MCMCpack)
library(parallel)
library(bggum)
library(pROC)
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Ensure required directories exist
prepare_directories()
## Compile functions to get log likelihoods:
Rcpp::sourceCpp("code/log_likelihood.cpp")


##### Simulate data -----
## Here's a general purpose function to simulate responses from the CJR model
## given its parameters
cjr_simulation <- function(theta, beta, alpha) {
    n <- length(theta)
    m <- length(beta)
    result <- matrix(0, nrow = n, ncol = m)
    zero_prob <- cjrProbability(result, theta, alpha, beta)
    for ( i in 1:n ) {
        for ( j in 1:m ) {
            u <- runif(n = 1, min = 0, max = 1)
            result[i, j] <- as.integer(u > zero_prob[i, j])
        }
    }
    return(result)
}
## Create the simulated data
n <- 100
m <- 400
set.seed(123)
theta <- rnorm(n)
beta  <- rnorm(m)
alpha <- rnorm(m, sd = 0.75)
response_matrix <- cjr_simulation(theta, beta, alpha)
## Ensure there's no unanimous items
# apply(response_matrix, 2, table)


##### Generate posterior samples -----
## First we'll generate samples for the MC3-GGUM model
## Tune proposals and temperatures
set.seed(314)
ggum_sds <- tune_proposals(response_matrix, 5000)
set.seed(2718)
ggum_temps <- tune_temperatures(
    response_matrix,
    n_temps = 6,
    proposal_sds = ggum_sds
)
## Setup cluster to produce chains in parallel
n_cores <- 2
cl <- makeCluster(n_cores, type = "FORK", outfile = "output/monotonic-log.txt")
## Deal with reproducibility
clusterSetRNGStream(cl = cl, iseed = 1119)
## Produce the chains
ggum_chains <- parLapply(cl = cl, X = 1:2, fun = function(x) {
    ggumMC3(
        data = response_matrix,
        sample_iterations = 25000,
        burn_iterations = 2500,
        temps = ggum_temps,
        proposal_sds = ggum_sds
    )
})
## Save the raw output
save(ggum_chains, file = "output/ggum_chains-monotonic-sim.RData")
## and clean up
stopCluster(cl)
## The question is who to use as the constraint.
## I'll use a random respondent whose theta parameter is less than negative one
set.seed(138)
constraint_idx <- sample(x = which(theta < -1), size = 1)
processed_chains <- lapply(
    X = ggum_chains,
    FUN = post_process,
    constraint = constraint_idx,
    expected_sign = "-"
)
## Summarise posterior
ggum_posterior_summary <- summary(processed_chains)
## Save the summary
save(ggum_posterior_summary, file = "output/ggum_summary-monotonic-sim.RData")

## Now we'll generate samples for the CJR model
## I'll give the CJR sampler the most extreme people as constraints
lowest_theta <- which.min(theta)
highest_theta <- which.max(theta)
rownames(response_matrix) <- 1:n
colnames(response_matrix) <- 1:m
cjr_constraints <- list("-", "+")
names(cjr_constraints) <- as.character(c(lowest_theta, highest_theta))
cjr_theta_start <- rep(0, n)
cjr_theta_start[lowest_theta] <- -2
cjr_theta_start[highest_theta] <- 2
## Generate CJR samples
cjr_chains <- lapply(1:2, function(x) {
    MCMCirt1d(
        datamatrix = response_matrix,
        theta.constraints = cjr_constraints,
        burnin = 20000,
        mcmc = 200000,
        verbose = 20000,
        seed = x + 1,
        theta.start = cjr_theta_start,
        alpha.start = 0,
        beta.start = 1,
        AB0 = 1,
        store.item = TRUE,
        store.ability = TRUE
    )
})
## Save the raw output
save(cjr_chains, file = "output/cjr_chains-monotonic-sim.RData")
chain <- do.call(rbind, cjr_chains)
cjr_estimates <- apply(chain, 2, mean)
saveRDS(cjr_estimates, file = "output/cjr_estimates-monotonic-sim.rds")


##### Reproduce Table 1 -----
## Define convenience function to round & put in character representation
round_ <- function(...) sprintf(paste0("%0.2f"), ...)
## Get the log likelihood for and mean theta SD GGUM
ggum_loglik  <- ggum_log_likelihood(response_matrix, ggum_chains[[1]])
ggum_mean_sd <- mean(ggum_posterior_summary$sds$theta)
## Get the log likelihood for and mean theta SD CJR
cjr_loglik  <- cjr_log_likelihood(response_matrix, cjr_chains[[1]])
cjr_mean_sd <- mean(apply(chain[ , 1:n], 2, sd))
results <- data.frame(Model = c("CJR", "GGUM"))
results$`Log Likelihood` <- round(c(cjr_loglik, ggum_loglik))
N <- prod(dim(response_matrix))
results$`L / N` <- round_(results$`Log Likelihood` / N)
results$`Log Likelihood` <- round(results$`Log Likelihood`)
results$`Mean Theta SD` <- round_(c(cjr_mean_sd, ggum_mean_sd))
# print(results)
write.csv(results, file = "tables/Table1.csv", row.names = FALSE)
# round(cor(
#     ggum_posterior_summary$estimates$theta,
#     cjr_estimates[grepl("theta", names(cjr_estimates))]
# ), 3)


##### Reproduce Table D1 -----
## Define some fit stat functions
brier <- function(obs, predprobs) mean((predprobs - obs)^2, na.rm = TRUE)
apre  <- function(obs, predictions) {
    minority_votes <- apply(obs, 2, function(v) min(table(v)))
    classification_errors <- colSums(obs != predictions, na.rm = TRUE)
    return( sum(minority_votes - classification_errors) / sum(minority_votes) )
}
AUC   <- function(obs, predprobs) {
    return(auc(roc(
        response  = c(obs),
        predictor = c(predprobs),
        quiet = TRUE
    )))
}
## Build a matrix of ones to get predicted probabilities
J_mat <- matrix(1, nrow = n, ncol = m)
## Calculate predicted probabilities and classifications for GGUM
theta <- ggum_posterior_summary$estimates$theta
alpha <- ggum_posterior_summary$estimates$alpha
delta <- ggum_posterior_summary$estimates$delta
tau   <- ggum_posterior_summary$estimates$tau
ggum_probs <- ggumProbability(J_mat, theta, alpha, delta, tau)
ggum_class <- ifelse(ggum_probs > 0.5, 1, 0)
## Calculate predicted probabilities and classifications for CJR
cjr_theta <- cjr_estimates[grepl("theta", names(cjr_estimates))]
cjr_beta  <- cjr_estimates[grepl("beta",  names(cjr_estimates))]
cjr_alpha <- cjr_estimates[grepl("alpha", names(cjr_estimates))]
cjr_probs <- cjrProbability(J_mat, cjr_theta, cjr_alpha, cjr_beta)
cjr_class <- ifelse(cjr_probs > 0.5, 1, 0)
## Calculate fit statistics
cjr_correct  <- mean(response_matrix == cjr_class)
ggum_correct <- mean(response_matrix == ggum_class)
cjr_apre   <- apre(response_matrix, cjr_class)
ggum_apre  <- apre(response_matrix, ggum_class)
cjr_auc    <- AUC(response_matrix, cjr_probs)
ggum_auc   <- AUC(response_matrix, ggum_probs)
cjr_brier  <- brier(response_matrix, cjr_class)
ggum_brier <- brier(response_matrix, ggum_class)
## Build and print table
results <- data.frame(Model = c("CJR", "GGUM"))
results$`Proportion Correct` <- round_(c(cjr_correct, ggum_correct))
results$APRE  <- round_(c(cjr_apre,  ggum_apre))
results$AUC   <- round_(c(cjr_auc,   ggum_auc))
results$Brier <- round_(c(cjr_brier, ggum_brier))
results$`Log Likelihood` <- round(c(cjr_loglik, ggum_loglik))
results$`L / N` <- round_(results$`Log Likelihood` / N)
# print(results)
write.csv(results, file = "tables/TableD1.csv", row.names = FALSE)

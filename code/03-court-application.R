##### Setup -----
## Load required packages
library(bggum)     ## For GGUM-MC3 sampler and related tools
library(MCMCpack)  ## For the CJR model sampler
library(pROC)      ## For getting AUC (for appendix fit stats)
library(parallel)  ## For running chains in parallel
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Ensure required directories exist
prepare_directories()


##### Prepare data -----
# Read in the raw data from Martin & Quinn and the SCDB
scdb <- read.csv("data/SCDB_2018_02_justiceCentered_Citation.csv")
load("data/mqData2015.Rda")
# We need to subset everything to the 1704 natural court
scdb <- scdb[scdb$naturalCourt == 1704, ]
# Them we can get our response matrix
# We record the justice names for subsetting mqData (and other things later)
justices <- unique(scdb$justiceName)
# I find it useful to also know the Lexis cite for the cases we deal with
cases <- unique(scdb$lexisCite[scdb$caseId %in% mqData$caseId])
# But we'll use caseId to match cases
ids <- unique(scdb$caseId[scdb$caseId %in% mqData$caseId])
# Then we get the mqData columns just for the justices in the 1704 court
# and the caseId column, which we just use to get rownames then remove
response_matrix <- mqData[mqData$caseId %in% ids, c(justices, "caseId")]
rownames(response_matrix) <- response_matrix$caseId
# We also transpose to get item columns
response_matrix <- t(response_matrix[ , justices])
# Now we start getting the folds for the tenfold cv
# One problem is some cases have only one dissent;
# we lose cases if we let that one dissent vote go missing.
# So, we don't let that happen, but otherwise induce random missingness
# This finds the row and column of un-deletable observations
response_table <- apply(response_matrix, 2, table)
problem_columns <- which(apply(response_table, 2, function(x) any(x == 1)))
problem_rows <- apply(response_matrix[, problem_columns], 2, function(x) which(x == 0))
problems <- cbind(problem_rows, problem_columns)
# Then we turn them into vector-like indices
problems_idx <- apply(problems, 1, function(x) x[1] + 9*(x[2]-1))
# Now we can create the folds
# Set the seed for reproducibility
set.seed(817, sample.kind = "Rounding")
# Get the index numbers that are not problematic
idx <- setdiff(1:length(response_matrix), problems_idx)
# shuffle them
idx <- sample(idx)
# and split them into 10 groups
folds <- cut(1:length(idx), breaks = 10, labels = FALSE)
# Double check that it all works out:
for ( fold in 1:10 ) {
    tmp_data <- response_matrix
    tmp_data[idx[which(folds == fold)]] <- NA
    K <- apply(tmp_data, 2, function(x) length(unique(na.omit(x))))
    if ( any(K != 2) ) {
        print(paste("Fold", fold, "ends up 'unanimous'"))
        break
    }
}
# That's fine.
# Next, the following are just some useful housekeeping variables
n <- nrow(response_matrix)
m <- ncol(response_matrix)
time_map <- mqData$time[match(ids, mqData$caseId)]
time_map <- time_map - min(time_map) + 1


##### Run samplers for each fold -----
# Set up our cluster
number_of_cores <- detectCores() - 2 # or change as appropriate
cl <- makeCluster(number_of_cores, type = "FORK", outfile = "output/log.txt")
# We'll use the following function for each fold:
analyze_fold <- function(fold_number, save = FALSE) {
    # Copy the response matrix and generate the 1/10 missingness
    tmp_data <- response_matrix
    tmp_data[idx[which(folds == fold_number)]] <- NA
    # Set some hyper-parameters
    mq_theta_start <- rep(0, n)
    mq_theta_start[5] <- -2 # Ginsburg
    mq_theta_start[4] <- 2  # Thomas
    mq_theta_constraints <- list(RBGinsburg = "-", CThomas = "+")
    ggum_sds <- tune_proposals(tmp_data, 10000)
    ggum_temps <- tune_temperatures(tmp_data, 6, proposal_sds = ggum_sds)
    # Run the GGUM analysis
    ggum_chain1 <- ggumMC3(data = tmp_data,
                           sample_iterations = 200000,
                           burn_iterations = 20000,
                           proposal_sds = ggum_sds,
                           temps = ggum_temps)
    ggum_chain2 <- ggumMC3(data = tmp_data,
                           sample_iterations = 200000,
                           burn_iterations = 20000,
                           proposal_sds = ggum_sds,
                           temps = ggum_temps)
    # Run the Martin-Quinn analysis
    mq_chain1 <- MCMCdynamicIRT1d(tmp_data,
                                  item.time.map = time_map,
                                  theta.start = mq_theta_start,
                                  mcmc = 200000,
                                  burnin = 20000,
                                  tau2.start = rep(0.1, n),
                                  A0 = 1,
                                  B0 = 1,
                                  beta.start = 1,
                                  alpha.start = 0,
                                  theta.constraints = mq_theta_constraints,
                                  seed = fold_number,
                                  verbose = 20000)
    mq_chain2 <- MCMCdynamicIRT1d(tmp_data,
                                  item.time.map = time_map,
                                  theta.start = mq_theta_start,
                                  mcmc = 200000,
                                  burnin = 20000,
                                  tau2.start = rep(0.1, n),
                                  A0 = 1,
                                  B0 = 1,
                                  beta.start = 1,
                                  alpha.start = 0,
                                  theta.constraints = mq_theta_constraints,
                                  seed = fold_number + 1,
                                  verbose = 20000)
    # Run the CJR analysis
    cjr_chain1 <- MCMCirt1d(datamatrix = response_matrix,
                            theta.constraints = mq_theta_constraints,
                            burnin = 20000,
                            mcmc = 200000,
                            verbose = 0,
                            seed = fold_number,
                            theta.start = mq_theta_start,
                            alpha.start = 0,
                            beta.start = 1,
                            AB0 = 1,
                            store.item = TRUE,
                            store.ability = TRUE)
    cjr_chain2 <- MCMCirt1d(datamatrix = response_matrix,
                            theta.constraints = mq_theta_constraints,
                            burnin = 20000,
                            mcmc = 200000,
                            verbose = 0,
                            seed = fold_number + 1,
                            theta.start = mq_theta_start,
                            alpha.start = 0,
                            beta.start = 1,
                            AB0 = 1,
                            store.item = TRUE,
                            store.ability = TRUE)
    # Save the raw output
    if ( save ) {
        outfile <- paste0("output/fold", fold_number, "-raw-results.RData")
        save(
            ggum_sds, ggum_temps, ggum_chain1, ggum_chain2,
            mq_chain1, mq_chain2,
            cjr_chain1, cjr_chain2,
            file = outfile
        )
    }
    # Post-process GGUM output -- Ginsburg should be negative
    ggum_chain1 <- post_process(ggum_chain1, 5, "-")
    ggum_chain2 <- post_process(ggum_chain2, 5, "-")
    # Check convergence
    mcmclist <- mcmc.list(ggum_chain1, ggum_chain2)
    ggum_convergence <- try(gelman.diag(mcmclist), silent = TRUE)
    mcmclist <- mcmc.list(mq_chain1, mq_chain2)
    mq_convergence <- try(gelman.diag(mcmclist), silent = TRUE)
    mcmclist <- mcmc.list(cjr_chain1, cjr_chain2)
    cjr_convergence <- try(gelman.diag(mcmclist), silent = TRUE)
    # Get estimates
    ggum_summary <- summary(list(ggum_chain1, ggum_chain2))
    ggum_estimates <- ggum_summary$estimates
    ggum_posterior_sds <- ggum_summary$sds$theta_sds
    chain <- rbind(mq_chain1, mq_chain2)
    mq_estimates <- apply(chain, 2, mean)
    mq_posterior_sds <- apply(chain[ , grepl("theta", colnames(chain))], 2, sd)
    chain <- rbind(cjr_chain1, cjr_chain2)
    cjr_estimates <- apply(chain, 2, mean)
    cjr_posterior_sds <- apply(chain[ , grepl("theta", colnames(chain))], 2, sd)
    # Return everything
    return(
        list(
            fold_number = fold_number,
            ggum_temps = ggum_temps, # Not super important, but I like to see
            ggum_estimates = ggum_estimates,
            mq_estimates = mq_estimates,
            cjr_estimates = cjr_estimates,
            ggum_posterior_sds = ggum_posterior_sds,
            mq_posterior_sds = mq_posterior_sds,
            cjr_posterior_sds = cjr_posterior_sds,
            ggum_convergence = ggum_convergence,
            mq_convergence = mq_convergence,
            cjr_convergence = cjr_convergence
        )
    )
}
# First we'll run the samplers on the full data
full_estimates <- analyze_fold(0, save = TRUE)
saveRDS(full_estimates, file = "output/court-estimates.rds")
# Now we deal with some reproducibility issues (since we're going parallel)
clusterSetRNGStream(cl = cl, iseed = 1030)
# run the samplers
fold_estimates <- parLapplyLB(cl = cl, X = 1:10, analyze_fold)
# save the results
saveRDS(fold_estimates, file = "output/fold_estimates.rds")
# and clean up
stopCluster(cl)

##### Reproduce Table H1 -----
# This is the response probability function assumed for M-Q
mq_prob <- function(vote, alpha, beta, theta) {
    # MCMCpack documentation puts it like
    # return(pnorm(beta*theta - alpha, lower.tail = vote))
    # which is equivalent to
    return(pnorm(alpha - beta*theta, lower.tail = !vote))
}
# We'll store the probability of observed responses in these matrices
ggum_probs <- matrix(NA_real_, nrow = n, ncol = m)
mq_probs <- matrix(NA_real_, nrow = n, ncol = m)
# Now we get the probabilities using estimates from each fold
for ( fold in 1:10 ) {
    # Get which observations were removed for this fold
    removed <- idx[which(folds == fold)]
    # Get the results from this fold
    elem <- which(sapply(fold_estimates, '[[', "fold_number") == fold)
    fold_results <- fold_estimates[[elem]]
    # Get the estimates
    thetas <- fold_results$ggum_estimates$theta
    alphas <- fold_results$ggum_estimates$alpha
    deltas <- fold_results$ggum_estimates$delta
    taus <- fold_results$ggum_estimates$tau
    mq_estimates <- fold_results$mq_estimates
    # Find the probability of each actual vote according to the estimates
    for ( vote in removed ) {
        i <- vote %% 9 + 9*(vote %% 9 == 0)
        j <- vote %/% 9 + !(vote %% 9 == 0)
        if ( is.na(response_matrix[i, j]) ) {
            next()
        }
        ggum_probs[i, j] <- ggumProbability(response_matrix[i, j],
                                            thetas[i], alphas[j],
                                            deltas[j], taus[[j]])
        justice <- justices[i]
        tt <- paste0(".t", time_map[j])
        id <- paste0(".", colnames(response_matrix)[j])
        theta <- mq_estimates[paste0("theta.", justice, tt)]
        alpha <- mq_estimates[paste0("alpha", id)]
        beta <- mq_estimates[paste0("beta", id)]
        mq_probs[i, j] <- mq_prob(response_matrix[i, j], alpha, beta, theta)
    }
}
# Now we can get some other useful things, like predicted probabilities,
# classifications, etc.
ggum_pred_prob <- ifelse(response_matrix == 1, ggum_probs, 1 - ggum_probs)
mq_pred_prob <- ifelse(response_matrix == 1, mq_probs, 1 - mq_probs)
ggum_pred <- ifelse(ggum_pred_prob > 0.5, 1, 0)
mq_pred <- ifelse(mq_pred_prob > 0.5, 1, 0)
actual <- c(response_matrix)
naive_pred <- pmax(actual, 1)
save(ggum_pred, ggum_pred_prob, mq_pred, mq_pred_prob,
     file = "data/10fold-predicted-probabilites.Rda")
# Check the confusion matrices
# cat("GGUM confusion matrix:\n")
# caret::confusionMatrix(factor(c(ggum_pred)), factor(actual))
# cat("MQ confusion matrix:\n")
# caret::confusionMatrix(factor(c(mq_pred)), factor(actual))
# Get proportion correctly classified
ggum_correct <- mean(ggum_pred == response_matrix, na.rm = TRUE)
mq_correct <- mean(mq_pred == response_matrix, na.rm = TRUE)
naive_correct <- mean(naive_pred == response_matrix, na.rm = TRUE)
# Brier
ggum_brier <- mean((ggum_pred_prob - response_matrix)^2, na.rm = TRUE)
mq_brier <- mean((mq_pred_prob - response_matrix)^2, na.rm = TRUE)
# APRE
minority_votes <- apply(response_matrix, 2, function(x) sum(!x, na.rm = TRUE))
ggum_errors <- sapply(1:m, function(j) {
    sum(ggum_pred[ , j] != response_matrix[ , j], na.rm = TRUE)
})
mq_errors <- sapply(1:m, function(j) {
    sum(mq_pred[ , j] != response_matrix[ , j], na.rm = TRUE)
})
ggum_apre <- sum(minority_votes - ggum_errors) / sum(minority_votes)
mq_apre <- sum(minority_votes - mq_errors) / sum(minority_votes)
# AUC
ggum_auc <- auc(roc(actual, c(ggum_pred)))
mq_auc <- auc(roc(actual, c(mq_pred)))
# Collect the comparisons and save them
model_comparisons <- data.frame(model = c("GGUM", "MQ", "Naive"),
                                proportion_correct = c(ggum_correct,
                                                       mq_correct,
                                                       naive_correct),
                                brier = c(ggum_brier, mq_brier, NA),
                                apre = c(ggum_apre, mq_apre, NA),
                                auc = c(ggum_auc, mq_auc, NA))
write.csv(model_comparisons, file = "tables/TableH1.csv", row.names = FALSE)


##### Reproduce Figure 7 -----
## Pull out estimates for convenience
theta <- full_estimates$ggum_estimates$theta
alpha <- full_estimates$ggum_estimates$alpha
delta <- full_estimates$ggum_estimates$delta
tau   <- full_estimates$ggum_estimates$tau
## Make a function to get the response probabilities for the Martin-Quinn model
mq_prob <- function(vote, alpha, beta, theta) {
    return(pnorm(alpha - beta*theta, lower.tail = !vote))
}
## Set the theta range to use for the Martin-Quinn models
th_range <- seq(-3, 3, 0.01)
## Set palette
pal   <- c("#808080", "black")
## Store justices' initials to identify their position on the plots
j_id  <- c("JR", "AS", "AK", "CT", "RG", "SB", "SA", "SS", "EK")

## Plot GGUM IRF for Wynne
tiff(
    filename = "plots/fig7a.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
item <- 166
custom_irf(
    alpha[item], delta[item], tau[[item]], from = -2.25, to = 2.25,
    option_names = c("Dissent", "Majority"), line_types = 2:1, color = TRUE,
    color_palette = pal, main_title = "", x_axis_at = -2:2, lwd = 2
)
resps <- response_matrix[ , item]
probs <- ggumProbability(resps, theta, alpha[item], delta[item], tau[[item]])
points(theta, probs, pch = 19, col = pal[resps + 1], cex = 1.5)
#            JG     AS    AM     CT     RG     SB     SA     SS     EK
xadj <- c(-0.08, -0.08, 0.08,  0.08,  0.08, -0.08,  0.00, -0.20, -0.08)
yadj <- c(-0.08,  0.08, 0.08, -0.08,  0.08,  0.08, -0.08,  0.00, -0.08)
text(theta + xadj, probs + yadj, labels = j_id, col = pal[resps + 1])
par(opar)
invisible(dev.off())

## Plot M-Q IRF for Wynne
tiff(
    filename = "plots/fig7b.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
item <- 166
mq_estimates <- full_estimates$mq_estimates
theta <- mq_estimates[grepl(paste0("t", time_map[item]), names(mq_estimates))]
alpha <- mq_estimates[paste("alpha", ids[item], sep = ".")]
beta  <- mq_estimates[paste("beta", ids[item], sep = ".")]
resps <- response_matrix[ , item]
probs <- sapply(1:9, function(i) mq_prob(resps[i], alpha, beta, theta[i]))
ones  <- mq_prob(1, alpha, beta, th_range)
zeros <- mq_prob(0, alpha, beta, th_range)
plot(
    x = th_range, y = zeros,
    type = "l", lty = 2, col = pal[1], lwd = 2,
    ylim = c(0, 1.01), xaxt = "n", yaxt = "n",
    main = "", xlab = "", ylab = ""
)
axis(
    side = 2,
    at = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    tick = FALSE, line = -0.75
)
title(ylab = expression(P[ij](k)), line = 1.5)
axis(1, at = -3:3, tick = FALSE, line = -0.75)
mtext(expression("Ideology (" * theta * ")"), side = 1, line = 1.5)
lines(th_range, ones, col = pal[2], lwd = 2)
points(theta, probs, pch = 19, col = pal[resps+1], cex = 1.5)
#            JG     AS    AM     CT     RG     SB     SA     SS     EK
xadj <- c(-0.00, -0.00, 0.00,  0.00,  0.00, -0.00,  0.00, -0.00, -0.00)
yadj <- c( 0.08, -0.08, 0.08, -0.08, -0.08,  0.08,  0.08,  0.08, -0.08)
text(theta + xadj, probs + yadj, labels = j_id, col = pal[resps+1])
par(xpd = TRUE)
legend(
    "topleft", inset = c(0, -0.1), legend = c("Dissent", "Majority"),
    lty = 2:1, horiz = TRUE, col = pal, bty = "n", lwd = 2
)
par(opar)
invisible(dev.off())


##### Reproduce Figure 8 -----
## Pull out estimates for convenience
theta <- full_estimates$ggum_estimates$theta
alpha <- full_estimates$ggum_estimates$alpha
delta <- full_estimates$ggum_estimates$delta
tau   <- full_estimates$ggum_estimates$tau
## Make a function to get the response probabilities for the Martin-Quinn model
mq_prob <- function(vote, alpha, beta, theta) {
    return(pnorm(alpha - beta*theta, lower.tail = !vote))
}
## Set the theta range to use for the Martin-Quinn models
th_range <- seq(-3, 3, 0.01)
## Set palette
pal   <- c("#808080", "black")
## Store justices' initials to identify their position on the plots
j_id  <- c("JR", "AS", "AK", "CT", "RG", "SB", "SA", "SS", "EK")

## Plot GGUM IRF for Arizona
tiff(
    filename = "plots/fig8a.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
item <- 83
custom_irf(
    alpha[item], delta[item], tau[[item]], from = -2.25, to = 2.25,
    option_names = c("Dissent", "Majority"), line_types = 2:1, color = TRUE,
    color_palette = pal, main_title = "", x_axis_at = -2:2, lwd = 2
)
resps <- response_matrix[ , item]
probs <- ggumProbability(resps, theta, alpha[item], delta[item], tau[[item]])
points(theta, probs, pch = 19, col = pal[resps + 1], cex = 1.5)
#            JG     AS    AM     CT     RG    SB    SA    SS     EK
xadj <- c(-0.20, -0.12, 0.08, -0.04, -0.00,  0.00,  0.08,  0.00, -0.08)
yadj <- c(-0.00,  0.08, 0.08,  0.08, -0.08, -0.08, -0.08, -0.08, -0.08)
text(theta + xadj, probs + yadj, labels = j_id, col = pal[resps + 1])
par(opar)
invisible(dev.off())

## Plot M-Q IRF for Arizona
tiff(
    filename = "plots/fig8b.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
item <- 83
mq_estimates <- full_estimates$mq_estimates
theta <- mq_estimates[grepl(paste0("t", time_map[item]), names(mq_estimates))]
alpha <- mq_estimates[paste("alpha", ids[item], sep = ".")]
beta  <- mq_estimates[paste("beta", ids[item], sep = ".")]
resps <- response_matrix[ , item]
probs <- sapply(1:9, function(i) mq_prob(resps[i], alpha, beta, theta[i]))
ones  <- mq_prob(1, alpha, beta, th_range)
zeros <- mq_prob(0, alpha, beta, th_range)
plot(
    x = th_range, y = zeros,
    type = "l", lty = 2, col = pal[1], lwd = 2,
    ylim = c(0, 1.01), xaxt = "n", yaxt = "n",
    main = "", xlab = "", ylab = ""
)
axis(
    side = 2,
    at = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    tick = FALSE, line = -0.75
)
title(ylab = expression(P[ij](k)), line = 1.5)
axis(1, at = -3:3, tick = FALSE, line = -0.75)
mtext(expression("Ideology (" * theta * ")"), side = 1, line = 1.5)
lines(th_range, ones, col = pal[2], lwd = 2)
points(theta, probs, pch = 19, col = pal[resps+1], cex = 1.5)
#            JG     AS    AM     CT     RG    SB    SA    SS
xadj <- c( 0.20, -0.08, 0.08, -0.08, -0.00,  0.00,  0.08,  0.00)
yadj <- c( 0.08,  0.08, 0.08,  0.08, -0.08, -0.08, -0.08, -0.08)
text(theta[-9] + xadj, probs[-9] + yadj, labels = j_id[-9], col = pal[resps+1])
par(xpd = TRUE)
legend(
    "topleft", inset = c(0, -0.1), legend = c("Dissent", "Majority"),
    lty = 2:1, horiz = TRUE, col = pal, bty = "n", lwd = 2
)
par(opar)
invisible(dev.off())

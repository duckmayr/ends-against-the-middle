setwd("~/repos/ends-against-the-middle/")
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
number_of_cores <- 2 # or change as appropriate
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
# Now we deal with some reproducibility issues (since we're going parallel)
clusterSetRNGStream(cl = cl, iseed = 1030)
# run the samplers
fold_estimates <- parLapplyLB(cl = cl, X = 1:10, analyze_fold)
# save the results
saveRDS(fold_estimates, file = "output/fold_estimates.rds")
# and clean up
stopCluster(cl)

## Notify of finishing
write(x = "done", file = "~/Dropbox/done.txt")


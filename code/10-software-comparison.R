##### Setup -----
## Load required packages
library(dplyr)
library(tidyr)
library(bggum)
library(coda)
library(GGUM)
library(parallel)


##### Set up replication cells -----
n_respondents <- rep(c(100, 500, 1000), 4)
n_items <- rep(c(10, 20), each = 6)
n_options <- rep(c(2, 4), each = 3, times = 2)
cells <- cbind(n_respondents, n_items, n_options)
omit_rows <- c(4, 10) # Only 100 respondents with 4 options may be problematic
cells <- cells[-omit_rows, ]


##### Run analysis for each replicate in each cell -----
## Set up our cluster
n_cores <- 5 ## or change as appropriate
cl <- makeCluster(n_cores, type = "FORK", outfile = "output/mml-log.txt")
## Deal with reproducibility
clusterSetRNGStream(cl = cl, iseed = 327)
## Run the replicates
results <- parLapplyLB(cl = cl, X = 1:100, fun = function(x) {
    cell <-  (x - 1) %/% 10 + 1
    replication <- (x - 1) %% 10 + 1
    n <- cells[cell, 1]
    m <- cells[cell, 2]
    K <- cells[cell, 3]
    sim_data <- ggum_simulation(n, m, K, delta_params = c(2, 2, -3, 3),
                                alpha_params = c(1.5, 1.5, 0, 3))
    Z <- apply(sim_data$response_matrix, 2, function(x) length(unique(x)))
    if ( any(Z < K) ) {
        return(paste("Not all options were chosen for at least one question",
                     "for cell", cell, " replicate ", replication))
    }
    sds <- tune_proposals(sim_data$response_matrix, 5000)
    temps <- tune_temperatures(sim_data$response_matrix, n_temps = 5)
    ggum_chain1 <- ggumMC3(sim_data$response_matrix,
                           burn_iterations = 5000,
                           sample_iterations = 20000,
                           temps = temps, proposal_sds = sds)
    ggum_chain2 <- ggumMC3(sim_data$response_matrix,
                           burn_iterations = 5000,
                           sample_iterations = 20000,
                           temps = temps, proposal_sds = sds)
    aic <- which.min(sim_data$theta)
    processed_chains <- lapply(list(ggum_chain1, ggum_chain2), post_process,
                               constraint = aic, expected_sign = "-")
    ggum_summary <- summary(processed_chains)
    mc3_estimates <- ggum_summary$estimates
    ggum_mml <- GGUM(data = sim_data$response_matrix, C = K - 1)
    ggum_eap_thetas <- Theta.EAP(ggum_mml)
    out_name <- paste0("output/cell", cell, "r", replication, "results.RData")
    save(sim_data, ggum_chain1, ggum_chain2, ggum_summary,
         ggum_mml, ggum_eap_thetas, mc3_estimates, file = out_name)
    rmse <- function(x, y) return(sqrt(mean((x - y)^2, na.rm = TRUE)))
    parameter_names <- c("theta", "alpha", "delta", "tau")
    result <- data.frame(cell = rep(cell, 8),
                         replicate = rep(replication, 8),
                         parameter = rep(parameter_names, 2),
                         model = c(rep("mc3", 4), rep("mml/eap", 4)),
                         rmse = rep(NA_real_, 8))
    result[1, 5] <- rmse(sim_data$theta, mc3_estimates$theta)
    result[5, 5] <- rmse(sim_data$theta, ggum_eap_thetas[ , 2])
    result[2, 5] <- rmse(sim_data$alpha, mc3_estimates$alpha)
    result[6, 5] <- rmse(sim_data$alpha, ggum_mml$alpha)
    result[3, 5] <- rmse(sim_data$delta, mc3_estimates$delta)
    result[7, 5] <- rmse(sim_data$delta, ggum_mml$delta)
    tau_true_values <- unlist(sim_data$tau)
    tau_true_values <- tau_true_values[which(tau_true_values != 0)]
    mc3_taus <- unlist(mc3_estimates$tau)
    mc3_taus <- mc3_taus[which(mc3_taus != 0)]
    mml_taus <- c(t(ggum_mml$taus[ , 1:(K-1)]))
    result[4, 5] <- rmse(tau_true_values, mc3_taus)
    result[8, 5] <- rmse(tau_true_values, mml_taus)
    saveRDS(result, file = sub("RData", "rds", out_name))
})


##### Reproduce Table C.2 -----
params <- c("theta", "alpha", "delta", "tau")
files <- list.files(path = "output", pattern = "cell.*rds$", full.names = TRUE)
results <- do.call("rbind", lapply(files, readRDS))
results <- results %>%
    mutate(
        parameter = factor(parameter, levels = params),
        model = factor(model, levels = c("mml/eap", "mc3"))
    ) %>% 
    group_by(model, parameter) %>%
    summarise(val = mean(as.numeric(rmse))) %>% 
    pivot_wider(names_from = model, values_from = val)
# results
results <- format(results, digits = 2)
write.csv(results, file = "tables/TableC2.csv", row.names = FALSE)

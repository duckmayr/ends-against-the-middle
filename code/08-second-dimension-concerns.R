##### Setup -----
## Since the last script had to use a variant of the `bggum` package,
## we ensure that version is no longer loaded
if ( "bggum" %in% loadedNamespaces() ) {
    detach(name = "package:bggum", unload = TRUE, character.only = TRUE)
}
## Now we load the packages we need
library(wnominate)
library(MCMCpack)
library(parallel)
library(bggum)
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Ensure required directories exist
prepare_directories()


##### Simulation data prep -----
## Function giving probability of one response under multi-dimensional 2PL
irt2pl_prob <- function(theta, alpha, delta) {
    # theta is n x d
    # alpha is m x d
    # delta is m x 1
    result <- theta %*% t(alpha)
    result <- result + matrix(delta, nrow = n, ncol = m, byrow = TRUE)
    result <- plogis(result)
    return(result)
}

## Function that simulates multidimensional IRT responses using above function
irt2d_sim <- function(n, m, d, upweight_first_dim = TRUE) {
    theta <- round(matrix(rnorm(n*d), ncol = d), 3)
    alpha <- round(matrix(rnorm(m*d), ncol = d), 3)
    delta <- round(matrix(rnorm(m), ncol = 1), 3)
    if ( upweight_first_dim ) {
        if ( is.numeric(upweight_first_dim) ) {
            weight <- upweight_first_dim
        } else {
            weight <- 2.0
        }
        ## Make the first dimension worth more
        alpha[ , 1] <- weight * alpha[ , 1]
    }
    probs <- irt2pl_prob(theta, alpha, delta)
    result <- sapply(1:m, function(j) {
        sapply(1:n, function(i) {
            u <- runif(1, min = 0, max = 1)
            return("if"(u < probs[i, j], 0, 1))
        })
    })
    class(result) <- "irt2d_sim"
    attr(result, "theta") <- theta
    attr(result, "alpha") <- alpha
    attr(result, "delta") <- delta
    return(result)
}

## Simulate the data
set.seed(123)
n <- 100
m <- 400
d <- 2
sim_data <- irt2d_sim(n, m, d)
theta1d <- attr(sim_data, "theta")[ , 1]
theta2d <- attr(sim_data, "theta")[ , 2]


##### Fit models for simulation -----
## MC3-GGUM
set.seed(123)
proposal_sds <- tune_proposals(sim_data, tune_iterations = 5000)
temps <- tune_temperatures(sim_data, n_temps = 6, proposal_sds = proposal_sds)
number_of_cores <- 2
cl <- makeCluster(number_of_cores, type = "FORK",
                  outfile = "output/2d-sim-mc3-log.txt")
clusterSetRNGStream(cl = cl, iseed = 1119) # parallel reproducibility
chains <- parLapply(cl = cl, X = 1:2, fun = function(x) {
    ggumMC3(sim_data, sample_iterations = 50000, burn_iterations = 5000,
            temps = temps, proposal_sds = proposal_sds)
})
saveRDS(chains, file = "output/2d-sim-ggum-chains.rds")
raw_chains <- chains
aic <- which.max(theta1d)
chains <- lapply(chains, post_process, constraint = aic, expected_sign = "+")
mc3ggum_summary <- summary(chains)
ggum_theta <- mc3ggum_summary$estimates$theta
saveRDS(mc3ggum_summary, file = "output/2d-sim-mc3ggum-summary.rds")

## 2D IRT
sim_alpha   <- attr(sim_data, "alpha")
constraints <- list(
    list(2, "+"), list(2, "-"), list(2, 0),
    list(3, "+"), list(3, "-"), list(3, 0)
)
names(constraints) <- c(
    paste0("V", which.max(sim_alpha[ , 1])),
    paste0("V", which.min(sim_alpha[ , 1])),
    paste0("V", which.min(abs(sim_alpha[ , 1]))),
    paste0("V", which.max(sim_alpha[ , 2])),
    paste0("V", which.min(sim_alpha[ , 2])),
    paste0("V", which.min(abs(sim_alpha[ , 2])))
)
cjr_chain   <- MCMCirtKd(
    datamatrix = sim_data,
    dimensions = 2,
    burnin = 5000,
    mcmc = 50000,
    B0 = 0.25,
    seed = 1234,
    item.constraints = constraints,
    verbose = 5000,
    store.item = TRUE
)
saveRDS(cjr_chain, file = "output/2d-sim-cjr-chain.rds")
cjr_summary <- summary(cjr_chain)
saveRDS(cjr_summary, file = "output/2d-sim-cjr-summary.rds")
parnames  <- rownames(cjr_summary$statistics)
dim1_rows <- which(grepl(pattern = "theta\\.[0-9]+\\.1", x = parnames))
cjr_dim1  <- -cjr_summary$statistics[dim1_rows, "Mean"]
dim2_rows <- which(grepl(pattern = "theta\\.[0-9]+\\.2", x = parnames))
cjr_dim2  <- -cjr_summary$statistics[dim2_rows, "Mean"]
# head(cjr_chain[,"beta.V160.1"])
cjr_means <- cjr_summary$statistics[ , "Mean"]
cjr_theta <- cbind(
    cjr_means[grepl("theta.+1$", names(cjr_means))],
    cjr_means[grepl("theta.+2$", names(cjr_means))]
)
cjr2d_d1  <- cjr_theta[ , 1]
cjr2d_d2  <- cjr_theta[ , 2]
cjr_alpha <- cjr_means[grepl("alpha", names(cjr_means))]
cjr_beta1 <- cjr_means[grepl("beta.+1$", names(cjr_means))]
cjr_beta2 <- cjr_means[grepl("beta.+2$", names(cjr_means))]
j1 <- which(grepl("beta.V130.2",  names(cjr_beta2)))
j2 <- which(grepl("beta.V125.1", names(cjr_beta1)))
cjr_beta1 <- c(cjr_beta1[1:(j1-1)], 0, cjr_beta1[j1:(length(cjr_alpha)-1)])
cjr_beta2 <- c(cjr_beta2[1:(j2-1)], 0, cjr_beta2[j2:(length(cjr_alpha)-1)])
cjr_beta  <- cbind(cjr_beta1, cjr_beta2)

## 1D IRT
theta.constraints <- list("+", "-")
names(theta.constraints) <- c(which.max(theta1d), which.min(theta1d))
cjr_chain   <- MCMCirt1d(
    datamatrix = sim_data,
    burnin = 5000,
    mcmc = 50000,
    seed = 1234,
    theta.constraints = theta.constraints,
    verbose = 5000,
    store.item = TRUE
)
saveRDS(cjr_chain, file = "output/2d-sim-1d-cjr-chain.rds")
cjr_1d_summary <- summary(cjr_chain)
saveRDS(cjr_1d_summary, file = "output/2d-sim-1d-cjr-summary.rds")

## W-NOMINATE
nominate1d <- wnominate(rollcall(unclass(sim_data)), dims = 1, polarity = aic)
saveRDS(nominate1d, file = "data/2d-sim-nominate1d.rds")
nom1d_d1 <- nominate1d$legislators$coord1D

aic <- c(aic, which.max(attr(sim_data, "theta")[ , 2]))
nominate2d <- wnominate(rollcall(unclass(sim_data)), dims = 2, polarity = aic)
saveRDS(nominate2d, file = "data/2d-sim-nominate2d.rds")
nom2d_d1 <- nominate2d$legislators$coord1D
nom2d_d2 <- nominate2d$legislators$coord2D


##### Reproduce Figure F1 -----
true_theta <- attr(sim_data, "theta")
parmat <- cbind(
    theta1   = true_theta[ , 1],
    GGUM     = ggum_theta,
    nom1d_d1 = nom1d_d1,
    nom2d_d1 = nom2d_d1,
    cjr2d_d1 = -cjr2d_d1,
    theta2   = true_theta[ , 2],
    nom2d_d2 = nom2d_d2,
    cjr2d_d2 = -cjr2d_d2
)
cormat <- cor(parmat)
cormat[lower.tri(cormat)] <- NA
cordat <- cormat %>%
    data.frame() %>%
    mutate(model2 = colnames(parmat)) %>%
    pivot_longer(!model2, names_to = "model1", values_to = "cor")

## Change model label variables so they're ordered how I want
model_names_ordered <- c(
    "theta1", "GGUM", "nom1d_d1", "nom2d_d1", "cjr2d_d1",
    "theta2", "nom2d_d2", "cjr2d_d2"
)
cordat <- cordat %>%
    mutate(
        model1 = factor(model1, levels = model_names_ordered),
        model2 = factor(model2, levels = model_names_ordered)
    )

## Plot correlation matrix
varlabs <- c(
    expression(theta[1]), "GGUM", "1D W-N D1", "2D W-N D1", "2D CJR D1",
    expression(theta[2]), "2D W-N D2", "2D CJR D2"
)
plt = ggplot(data = cordat, mapping = aes(x = model1, y = model2, fill = cor)) +
    geom_tile(color = "white", na.rm = TRUE) +
    scale_x_discrete(labels = varlabs) +
    scale_y_discrete(labels = varlabs) +
    scale_fill_gradient2(low = okabe_ito()[6],
                         mid = "white",
                         high = okabe_ito()[5],
                         na.value = "white",
                         midpoint = 0,
                         limits = c(-1, 1)) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.4, 0.6),
        legend.direction = "horizontal") +
    guides(fill = guide_colorbar(barwidth = 7,
                                 barheight = 1,
                                 title = "Correlation",
                                 title.position = "top",
                                 title.hjust = 0.5))
ggsave(
    plot = plt, filename = "plots/figF1.pdf",
    device = "pdf", height = 7, width = 7
)


##### Reproduce Figure F2 -----
pdf("plots/figF2.pdf", height = 6, width = 6)
opar <- par(mar = c(0, 0, 0, 0) + 0.1)
pairs(parmat, labels = varlabs, pch = 19, col = paste0(okabe_ito()[5], "80"))
par(opar)
invisible(dev.off())


##### Reproduce Figure F3 -----
pdf("plots/figF3.pdf", height = 4, width = 4)
trimplot(
    x = nom2d_d2, y = abs(ggum_theta),
    xlab = "W-NOMINATE 2nd Dim.",
    ylab = bquote(group("|", " GGUM "*theta*" ", "|")),
    pch = 19, col = "#80808080"
)
invisible(dev.off())


##### Reproduce Table F1 -----
## Define some fit stat functions
brier <- function(obs, predprobs) mean((predprobs - obs)^2, na.rm = TRUE)
apre  <- function(responses, predictions) {
    minority_votes <- apply(responses, 2, function(v) min(table(v)))
    classification_errors <- colSums(responses != predictions, na.rm = TRUE)
    return( sum(minority_votes - classification_errors) / sum(minority_votes) )
}
## Setup object to hold results for table
tab <- data.frame(
    Model = c("GGUM", "1D CJR", "2D CJR", "1D W-NOMINATE", "2D W-NOMINATE"),
    `Proportion Correct` = "", APRE = "", AUC = "", Brier = "",
    check.names = FALSE
)
## It will be useful to have a matrix of ones to get predicted probabilities
J <- matrix(1, nrow = nrow(sim_data), ncol = ncol(sim_data))
## Now get the fit statistics from each model
i <- which(tab$Model == "GGUM")
theta <- mc3ggum_summary$estimates$theta
alpha <- mc3ggum_summary$estimates$alpha
delta <- mc3ggum_summary$estimates$delta
tau   <- mc3ggum_summary$estimates$tau
P <- ggumProbability(J, theta, alpha, delta, tau)
C <- ifelse(P > 0.5, 1, 0)
tab$`Proportion Correct`[i] <- sprintf("%0.2f", mean(C == sim_data))
tab$APRE[i]  <- sprintf("%0.2f", apre(sim_data, C))
tab$AUC[i]   <- sprintf("%0.2f", auc(roc(c(sim_data), c(P))))
tab$Brier[i] <- sprintf("%0.2f", brier(sim_data, P))
i <- which(tab$Model == "1D CJR")
stats <- cjr_1d_summary$statistics
theta <- stats[grepl("theta", rownames(stats)), "Mean"]
alpha <- stats[grepl("alpha", rownames(stats)), "Mean"]
beta  <- stats[grepl("beta",  rownames(stats)), "Mean"]
P <- cjrProbability(J, theta, alpha, beta)
C <- ifelse(P > 0.5, 1, 0)
tab$`Proportion Correct`[i] <- sprintf("%0.2f", mean(C == sim_data))
tab$APRE[i]  <- sprintf("%0.2f", apre(sim_data, C))
tab$AUC[i]   <- sprintf("%0.2f", auc(roc(c(sim_data), c(P))))
tab$Brier[i] <- sprintf("%0.2f", brier(sim_data, P))
i <- which(tab$Model == "2D CJR")
P <- cjrProbability(J, cjr_theta, cjr_alpha, cjr_beta)
C <- ifelse(P > 0.5, 1, 0)
tab$`Proportion Correct`[i] <- sprintf("%0.2f", mean(C == sim_data))
tab$APRE[i]  <- sprintf("%0.2f", apre(sim_data, C))
tab$AUC[i]   <- sprintf("%0.2f", auc(roc(c(sim_data), c(P))))
tab$Brier[i] <- sprintf("%0.2f", brier(sim_data, P))
i <- which(tab$Model == "1D W-NOMINATE")
fits <- nominate1d$fits
tab$`Proportion Correct`[i] <- sprintf("%0.2f", fits["correctclass1D"] / 100)
tab$APRE[i] <- sprintf("%0.2f", fits["apre1D"])
i <- which(tab$Model == "2D W-NOMINATE")
fits <- nominate2d$fits
tab$`Proportion Correct`[i] <- sprintf("%0.2f", fits["correctclass2D"] / 100)
tab$APRE[i] <- sprintf("%0.2f", fits["apre2D"])
## Display/save the table
# tab
write.csv(tab, file = "tables/TableF1.csv", row.names = FALSE)


##### 92nd Senate data prep -----
## We analyze the second session of the 92nd Senate, so we first need to get
## the roll numbers of the relevant votes
rollnumbers <- read.csv("data/S092_rollcalls.csv") %>%
    filter(date >= "1972-01-01" & date <= "1972-12-31") %>%
    pull(rollnumber)
## Now we can construct the response matrix in the same manner as in
## 02-house-application.R (see there for further details on the process)
responses <- read.csv("data/S092_votes.csv") %>%
    filter(rollnumber %in% rollnumbers) %>%
    mutate(response = case_when(
        cast_code == 1 ~ 1L,
        cast_code == 6 ~ 0L,
        TRUE ~ NA_integer_
    )) %>%
    select(icpsr, rollnumber, response) %>%
    pivot_wider(names_from = "rollnumber", values_from = "response") %>%
    column_to_rownames("icpsr") %>%
    as.matrix()
lopsided  <- which(apply(responses, 2, is_lopsided, cutoff = 0.025))
responses <- responses[ , -lopsided]
few_votes <- which(apply(responses, 1, has_few_votes))
responses <- responses[-few_votes, ]
unanimous <- which(apply(responses, 2, is_unanimous))
responses <- responses[ , -unanimous]


##### Fit model for 92nd Senate -----
## Set the discrimination parameter prior wider than the default
alpha_prior <- c(1.5, 1.5, 0.25, 8.0)
## Set the seed for reproducibility
set.seed(138)
## Tune hyperparameters
sds <- tune_proposals(responses, 5000, alpha_prior_params = alpha_prior)
temps <- tune_temperatures(
    data = responses,
    n_temps = 6,
    proposal_sds = sds,
    alpha_prior_params = alpha_prior
)
## Set up the cluster
n_cores <- 2
cl <- makeCluster(n_cores, type = "FORK", outfile = "output/S92-log.txt")
## Deal with reproducibility
clusterSetRNGStream(cl = cl, iseed = 1119)
## Produce the chains
chains0 <- parLapply(cl = cl, X = 1:2, fun = function(x) {
    ggumMC3(
        data = responses,
        sample_iterations = 50000,
        burn_iterations = 5000,
        temps = temps,
        proposal_sds = sds,
        alpha_prior_params = alpha_prior
    )
})
# Save the raw output
save(chains0, file = "output/S92-chains.RData")
## Post-process; we'll use Ted Kennedy (ICPSR # 10808) for the constraint
aic <- which(rownames(responses) == "10808")
chains1 <- lapply(chains0, post_process, constraint = aic, expected_sign = "-")
posterior_summary <- summary(chains1)
saveRDS(posterior_summary, file = "output/S92-summary.rds")
theta <- posterior_summary$estimates$theta


##### Reproduce Figure 5 -----
iters <- nrow(chains1[[1]])
idx <- seq(from = 1, to = iters, by = floor(iters / 1000))
tiff(
    filename = "plots/fig5.tif",
    width = 12, height = 8, units = "in",
    res = 1000, compression = "lzw"
)
layout(matrix(1:4, ncol = 2))
lims  <- range(chains0[[1]][ , 21])
draws <- chains0[[1]][ , 21]
dens  <- density(draws)
trimplot(
    x = idx, y = draws[idx], type = "l", col = "#333333dd", ylim = lims,
    xlab = "Iteration", ylab = bquote("Ideology ("*theta*")"),
    main = "Barry Goldwater Theta Draws"
)
trimplot(
    x = dens, type = "l", xlim = lims, lwd = 2,
    ylab = "Density", xlab = bquote("Ideology ("*theta*")"),
    main = "Barry Goldwater Theta Density"
)
draws <- chains1[[1]][ , 21]
dens  <- density(draws)
trimplot(
    x = idx, y = draws[idx], type = "l", col = "#333333dd", ylim = lims,
    xlab = "Iteration", ylab = bquote("Ideology ("*theta*")"),
    main = "Barry Goldwater Theta Draws (Post-processed)"
)
trimplot(
    x = dens, type = "l", xlim = lims, lwd = 2,
    ylab = "Density", xlab = bquote("Ideology ("*theta*")"),
    main = "Barry Goldwater Theta Density (Post-processed)"
)
invisible(dev.off())


##### Reproduce Figure F4 -----
## Determine which senators were "Southern Democrats"
the_south <- c(
    "FL", "GA", "MD", "NC", "SC", "VA", "WV", "AL",
    "KY", "MS", "TN", "AR", "LA", "OK", "TX"
)
member_data <- read.csv("data/S092_members.csv") %>%
    mutate(southern_dem = state_abbrev %in% the_south & party_code == 100) %>%
    select(icpsr, bioname, nominate_dim1, party_code, southern_dem)
names(theta) <- rownames(responses)
ggum_df <- data.frame(theta = theta, icpsr = as.integer(names(theta)))
plot_data <- merge(ggum_df, member_data)
## Run 1-D W-NOMINATE on the data
aic <- which(rownames(responses) == "3658")
nominate1d <- wnominate(rollcall(responses), dims = 1, polarity = aic)
plot_data$nom1d_dim1 <- nominate1d$legislators$coord1D
## Plot the results
pdf("plots/figF4.pdf", height = 4, width = 6)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(
    x = plot_data$theta, y = plot_data$nom1d_dim1,
    col = ifelse(plot_data$southern_dem, "#ff0000", "#33333380"),
    pch = ifelse(plot_data$southern_dem, 19, 1),
    xlab = expression("GGUM Ideology (" * theta * ")"),
    ylab = "NOMINATE Dimension 1"
)
legend("topleft", legend = c("Southern Democrat", "Other"), bty = "n",
       col = c("#ff0000", "#333333"), pch = c(19, 1))
par(opar)
invisible(dev.off())


##### Reproduce Table F2 -----
## Run 2-D W-NOMINATE on the data
aic <- c(aic, which(rownames(responses) == "12100"))
nominate2d <- wnominate(rollcall(responses), dims = 2, polarity = aic)
tab <- data.frame(
    Model = c("GGUM", "1D W-NOMINATE", "2D W-NOMINATE"),
    `Proportion Correct` = "", APRE = "",
    check.names = FALSE
)
## It will be useful to have a matrix of ones to get predicted probabilities
J <- matrix(1, nrow = nrow(responses), ncol = ncol(responses))
## Now get the fit statistics from each model
i <- which(tab$Model == "GGUM")
theta <- posterior_summary$estimates$theta
alpha <- posterior_summary$estimates$alpha
delta <- posterior_summary$estimates$delta
tau   <- posterior_summary$estimates$tau
P <- ggumProbability(J, theta, alpha, delta, tau)
C <- ifelse(P > 0.5, 1, 0)
tab$`Proportion Correct`[i] <- sprintf("%0.2f", mean(C == responses, na.rm = T))
tab$APRE[i]  <- sprintf("%0.2f", apre(responses, C))
i <- which(tab$Model == "1D W-NOMINATE")
fits <- nominate1d$fits
tab$`Proportion Correct`[i] <- sprintf("%0.2f", fits["correctclass1D"] / 100)
tab$APRE[i] <- sprintf("%0.2f", fits["apre1D"])
i <- which(tab$Model == "2D W-NOMINATE")
fits <- nominate2d$fits
tab$`Proportion Correct`[i] <- sprintf("%0.2f", fits["correctclass2D"] / 100)
tab$APRE[i] <- sprintf("%0.2f", fits["apre2D"])
## Display/save the table
# tab
write.csv(tab, file = "tables/TableF2.csv", row.names = FALSE)

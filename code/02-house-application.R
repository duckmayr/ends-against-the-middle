##### Setup -----
## Load required packages
library(bggum)     ## For GGUM-MC3 sampler and related tools
library(MCMCpack)  ## For the CJR model sampler
library(wnominate) ## To estimate W-NOMINATE
library(oc)        ## For Poole's Optimal Classification
library(pROC)      ## For getting AUC (for appendix fit stats)
library(dplyr)     ## For data manipulation and summarization
library(tidyr)     ## For data reshaping
library(tibble)    ## For extra data manipulation functionality not in dplyr
library(parallel)  ## For running chains in parallel
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Ensure required directories exist
prepare_directories()


##### Prepare data -----
## Read in raw data from Voteview (downloaded March 28, 2021)
member_data <- read.csv("data/H116_members.csv") ## Data on MCs
vote_data   <- read.csv("data/H116_votes.csv")   ## Data on members' votes
## Create a dichotomous response matrix from the data on members' votes
responses   <- vote_data %>%        ## First we dichotomize responses,
    mutate(                         ## with a "yea" vote being a "1",
        response = case_when(       ## "nay" being a "0", & everything
            cast_code == 1 ~ 1L,    ## else counted as missing. In theory
            cast_code == 6 ~ 0L,    ## we could leverage GGUM's polytomous
            TRUE ~ NA_integer_      ## nature & treat abstention, etc., as
        )                           ## substantively interesting responses,
    ) %>%                           ## but we leave that for the future.
    select(                         ## Now we need to move the data to "wide"
        icpsr, rollnumber, response ## form, where the rows represent MCs
    ) %>%                           ## and the columns represent roll call
    pivot_wider(                    ## votes. After eliminating the now
        names_from = "rollnumber",  ## superfluous columns of vote_data,
        values_from = "response"    ## we do this via the relatively new
    ) %>%                           ## pivot_wider() function (see its
    column_to_rownames("icpsr") %>% ## vignette for details). Finally we fix
    as.matrix()                     ## row names and convert to matrix.
## Eliminate Amash
responses <- responses[-which(rownames(responses) %in% c("21143", "91143")), ]
## Eliminate unanimous and lopsided votes and legislators with little data
lopsided  <- which(apply(responses, 2, is_lopsided))   ## find lopsided votes
responses <- responses[ , -lopsided]                   ## and remove them
few_votes <- which(apply(responses, 1, has_few_votes)) ## find MCs w > 90%
responses <- responses[-few_votes, ]                   ## missing & remove them
unanimous <- which(apply(responses, 2, is_unanimous))  ## find unanimous votes
responses <- responses[ , -unanimous]                  ## and remove them


##### Obtain and summarize MC3-GGUM posterior samples -----
## Before any sampling, we define our prior over alpha (the descrimination
## parameter); while the default settings should work just fine for delta
## (the location parameters) and tau (the option threshold parameters),
## as there are some sharp cutpoints on these bills we could have censoring
## without a wider prior on alpha; we double the default range.
alpha_prior <- c(1.5, 1.5, 0.25, 8.00)
## Set seed for reproducibility
set.seed(42)
## Tune proposal densities
sds <- tune_proposals(responses, 5000, alpha_prior_params = alpha_prior)
## Tune temperature schedule
temps <- tune_temperatures(responses, 6, proposal_sds = sds,
                           alpha_prior_params = alpha_prior)
## Get posterior samples; we want at least two chains to assess convergence,
## so we'll obtain those chains in parallel to reduce computation time.
cl <- makeCluster(2, type = "FORK", outfile = "output/H116-log.txt")
clusterSetRNGStream(cl = cl, iseed = 314)
ggum_chains <- parLapplyLB(cl = cl, X = 1:2, fun = function(x) {
    ggumMC3(data = responses,
            sample_iterations = 100000,
            burn_iterations = 10000,
            proposal_sds = sds,
            temps = temps,
            alpha_prior_params = alpha_prior)
})
stopCluster(cl)
saveRDS(ggum_chains, file = "output/H116-chains.rds")
## Post process to deal with reflection
aic <- which(rownames(responses) == "21726") ## Use Jayapal to identify
processed_ggum_chains <- lapply(ggum_chains, post_process, constraint = aic,
                                expected_sign = "-")
## Summarise posterior
ggum_posterior_summary <- summary(processed_ggum_chains)
saveRDS(ggum_posterior_summary, file = "output/H116-ggum-post-summary.rds")


##### Obtain and summarize CJR posterior samples -----
## Set up constraints, AOC and Jim Jordan
constraints <- list("21949" = "-", "20738" = "+")
## Generate a chain (using default priors)
cjr_chain   <- MCMCirt1d(
    datamatrix = responses,
    theta.constraints = constraints,
    burnin = 10000,
    mcmc = 100000,
    verbose = 10000,
    seed = 1234,
    store.item = TRUE
)
## Save output to disk (if needed for later inspection)
saveRDS(cjr_chain, file = "output/cjr-house-chain.rds")
## Summarise posterior
cjr_posterior_summary <- summary(cjr_chain)
saveRDS(cjr_posterior_summary, file = "output/cjr-post-summary.rds")


##### Reproduce Figure 10 -----
ggum_ideo <- ggum_posterior_summary$estimates$theta
cjr_means <- cjr_posterior_summary$statistics[ , "Mean"]
cjr_ideo  <- cjr_means[which(grepl("theta", names(cjr_means)))]
## Find the Squad
squad_names <- "OCASIO|TLAIB|PRESSLEY|OMAR"
squad_ids   <- member_data$icpsr[grepl(squad_names, member_data$bioname)]
icpsrs      <- rownames(responses)
in_squad    <- grepl(paste(squad_ids, collapse = "|"), icpsrs)
## Set point types and colors based on Squad & party membership
member_idx  <- match(icpsrs, member_data$icpsr)
is_dem      <- member_data$party_code[member_idx] == 100
point_types <- ifelse(in_squad, 17, ifelse(is_dem, 15, 16))
squad_col   <- "#000000"
dem_col     <- "#80808080"
rep_col     <- "#80808080"
point_cols  <- ifelse(in_squad, squad_col, ifelse(is_dem, dem_col, rep_col))
## Setup the plotting device:
tiff(
    filename = "plots/fig10.tif",
    width = 6, height = 4, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1)
## Plot the comparison
trimplot(
    x = cjr_ideo, y = ggum_ideo,
    xlab = "CJR Ideology", ylab = "GGUM Ideology",
    pch = point_types, col = point_cols
)
## Add a legend
legend(
    "topleft", pch = c(16, 15, 17), col = c(rep_col, dem_col, squad_col),
    bty = "n", legend = c("Republicans", "Democrats", "The Squad")
)
## Reset the margins
par(opar)
## And close the plotting device
invisible(dev.off())
## Look at outliers
comparison <- data.frame(
    icpsr = icpsrs,
    name  = member_data$bioname[member_idx],
    ggum_ideo = ggum_ideo,
    cjr_ideo  = cjr_ideo
)
comparison$diff <- comparison$ggum_ideo - comparison$cjr_ideo
reps <- comparison[comparison$ggum_ideo >  2.90 & comparison$cjr_ideo <  2.10, ]
dems <- comparison[comparison$ggum_ideo < -2.15 & comparison$cjr_ideo > -2.10, ]
reps[order(abs(reps$diff), decreasing = TRUE), ]
dems[order(abs(dems$diff), decreasing = TRUE), ]


##### Reproduce Figure 9 -----
## Extract parameters for convenience
theta <- ggum_posterior_summary$estimates$theta
alpha <- ggum_posterior_summary$estimates$alpha
delta <- ggum_posterior_summary$estimates$delta
tau   <- ggum_posterior_summary$estimates$tau
ths <- seq(from = round_out(min(theta)), to = round_out(max(theta)), by = 0.01)
N <- length(ths)

## Define axis labels, AOC index, & plot dims, which don't change between plots
xlabel <- expression(theta)
ylabel <- "Pr( Yea )"
aoc    <- which(rownames(responses) == "21949")
width  <- 5
height <- 3

## Plot IRF for H.J. Res. 31 vote (panel (a))
item <- which(colnames(responses) == "86")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
yat <- c(0, 0.5, 1)
yticklabs <- as.character(yat)
## Setup the plotting device:
tiff(
    filename = "plots/fig9a.tif",
    width = width, height = height, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1)
## Plot the IRF
plot(
    x = ths, y = probs, type = "l",
    xaxt = "n", yaxt = "n", xlab = "", ylab = ""
)
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75, at = yat, labels = yticklabs)
axis(side = 1, at = theta[aoc], labels = "AOC")
mtext(side = 1, text = xlabel, line =  1.50)
mtext(side = 2, text = ylabel, line =  1.50)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
## Reset the margins
par(opar)
## And close the plotting device
invisible(dev.off())

## Plot IRF for H.R. 2740 (panel (b))
item <- which(colnames(responses) == "366")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
yat <- c(0, 0.5, 1)
yticklabs <- as.character(yat)
## Setup the plotting device:
tiff(
    filename = "plots/fig9b.tif",
    width = width, height = height, units = "in",
    res = 1000, compression = "lzw"
)
## Edit the margins
opar <- par(mar = c(3, 3, 2.2, 1) + 0.1)
## Plot the IRF
plot(
    x = ths, y = probs, type = "l",
    xaxt = "n", yaxt = "n", xlab = "", ylab = ""
)
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75, at = yat, labels = yticklabs)
axis(side = 1, at = theta[aoc], labels = "AOC")
mtext(side = 1, text = xlabel, line =  1.50)
mtext(side = 2, text = ylabel, line =  1.50)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
## Reset the margins
par(opar)
## And close the plotting device
invisible(dev.off())


##### Reproduce Figures I.1-I.4 -----
## The remaining plots are for the appendix on non-monotonic votes

## Plot IRF for Figure I.1: Defense Funding
## (H.R. 2500 (National Defense Authorization Act))
item <- which(colnames(responses) == "472")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
pdf("plots/figI1.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())

## Plot IRFs for Figure I.2: Humanitarian Aid for Immigrants
## (Plot IRFs for H.R. 3401 (Border security appropriations))
## House version -- panel (a)
item <- which(colnames(responses) == "413")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
pdf("plots/figI2a.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())
## Senate version -- panel (b)
item <- which(colnames(responses) == "428")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
pdf("plots/figI2b.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())

## Plot IRF for Figure I.3: A 2 State Solution to the Israel-Palestine Conflict
## (H. Res. 326 (Two state solution to Israel-Palestine conflict))
item <- which(colnames(responses) == "651")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
pdf("plots/figI3.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())

## Plot IRF for I.4: The HEROES Act
## (H. Res. 866 (HEROES Act rule))
item <- which(colnames(responses) == "805")
probs <- ggumProbability(rep(1, N), ths, alpha[item], delta[item], tau[[item]])
pdf("plots/figI4.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(theta[responses[ , item] == 1], side = 3)
rug(theta[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())


##### Determine how many roll calls have non-monotonic IRFs -----
prop_nonmonotonic(ggum_posterior_summary)
# [1] 0.1678487


##### Reproduce Figure E.2 -----

## Run W-NOMINATE on 116th data
rc  <- rollcall(responses)
i1  <- which.max(member_data$nominate_dim1)
i2  <- which.max(member_data$nominate_dim2)
aic <- match(member_data[c(i1, i2), "icpsr"], rownames(responses))
nominate2d <- wnominate(rc, dims = 2, lop = 0.01, polarity = aic)
saveRDS(nominate2d, file = "output/H116-nominate2d.rds")

## Plot results
squad_col   <- "#80008080"
dem_col     <- "#00008080"
rep_col     <- "#80000080"
point_cols  <- ifelse(in_squad, squad_col, ifelse(is_dem, dem_col, rep_col))
party <- ifelse(in_squad, "Squad", ifelse(is_dem, "Dem", "Rep"))
pdf("plots/figE2.pdf")
opar <- par(mar = c(4, 4, 1, 1) + 0.1)
plot_coords(nominate2d, plotBy = party)
par(opar)
invisible(dev.off())


##### Reproduce Figure E.1 -----

## Run 2D CJR on 116th data
dim1j <- which.max(with(nominate2d$rollcalls, midpoint1D - spread1D))
dim2j <- which.max(with(nominate2d$rollcalls, midpoint2D - spread2D))
colnames(responses)[c(dim1j, dim2j)]
constraints <- list(
    "754" = list(2, "+"), ## Roll call 754 must have positive slope on 1st dim
    "754" = list(3,   0), ## Roll call 754 must load onto only 1st dim
    "43"  = list(3, "+"), ## Roll call  43 must have positive slope on 2nd dim
    "43"  = list(2,   0)  ## Roll call  43 must load onto only 2nd dim
)
cjr_chain   <- MCMCirtKd(
    datamatrix = responses,
    dimensions = 2,
    item.constraints = constraints,
    B0 = 0.25,
    burnin = 10000,
    mcmc = 100000,
    verbose = 10000,
    seed = 1234,
    store.item = TRUE
)
cjr2d_posterior_summary <- summary(cjr_chain)
saveRDS(cjr2d_posterior_summary, file = "output/cjr-2d-post-summary.rds")
parnames  <- rownames(cjr2d_posterior_summary$statistics)
dim1_rows <- which(grepl(pattern = "theta\\.[0-9]+\\.1", x = parnames))
cjr_dim1  <- -cjr2d_posterior_summary$statistics[dim1_rows, "Mean"]

## Run GGUM on 115th data
members115 <- read.csv("data/H115_members.csv")
votes115 <- read.csv("data/H115_votes.csv")
responses115 <- votes115 %>%        ## First we dichotomize responses,
    mutate(                         ## with a "yea" vote being a "1",
        response = case_when(       ## "nay" being a "0", & everything
            cast_code == 1 ~ 1L,    ## else counted as missing. In theory
            cast_code == 6 ~ 0L,    ## we could leverage GGUM's polytomous
            TRUE ~ NA_integer_      ## nature & treat abstention, etc., as
        )                           ## substantively interesting responses,
    ) %>%                           ## but we leave that for the future.
    select(                         ## Now we need to move the data to "wide"
        icpsr, rollnumber, response ## form, where the rows represent MCs
    ) %>%                           ## and the columns represent roll call
    pivot_wider(                    ## votes. After eliminating the now
        names_from = "rollnumber",  ## superfluous columns of vote_data,
        values_from = "response"    ## we do this via the relatively new
    ) %>%                           ## pivot_wider() function (see its
    column_to_rownames("icpsr") %>% ## vignette for details). Finally we fix
    as.matrix()                     ## row names and convert to matrix.
## Eliminate president
responses115 <- responses115[-which(grepl("^999", rownames(responses115))), ]
## Eliminate unanimous votes
unanimous <- which(apply(responses115, 2, is_unanimous))
responses115 <- responses115[ , -unanimous]
## Save ICPSRs
icpsrs115 <- rownames(responses115)
## Run MC3-GGUM, post-process, and summarise
set.seed(314)
chain115 <- ggumMC3(
    data = responses115,
    sample_iterations = 100000,
    burn_iterations = 10000,
    n_temps = 6,
    alpha_prior_params = alpha_prior
)
aic <- as.character(members115$icpsr[which.min(members115$nominate_dim1)])
aic <- which(rownames(responses115) == aic)
processed_ggum_chains <- lapply(ggum_chains, post_process, constraint = aic,
                                expected_sign = "-")
ggum_post <- summary(processed_ggum_chains)
saveRDS(ggum_post, file = "output/H115-ggum-post-summary.rds")


## Run W-NOMINATE on 115th data
h115rc      <- rollcall(responses115, legis.names = icpsrs115)
constraints <- c(which(icpsrs115 == 21721), which(icpsrs115 == 29127))
nom2d115    <- wnominate(h115rc, polarity = constraints)

## Panel (a): 2D CJR Dim. 1 vs. GGUM
pdf("plots/figE1a.pdf", width = 6, height = 4)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(cjr_dim1, ggum_ideo, xlab = "2D CJR Dim. 1", ylab = "GGUM Ideology",
         pch = point_types, col = point_cols)
legend("topleft", pch = 15:17, col = c(dem_col, rep_col, squad_col),
       bty = "n", legend = c("Republicans", "Democrats", "The Squad"))
par(opar)
invisible(dev.off())

## Panel (b): 2D DW-NOMINATE Dim. 1 vs. GGUM
dw_nom_ideo <- with(member_data, nominate_dim1[match(icpsrs, icpsr)])
pdf("plots/figE1b.pdf", width = 6, height = 4)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(
    x = dw_nom_ideo, y = ggum_ideo,
    xlab = "2D DW-NOMINATE Dim. 1", ylab = "GGUM Ideology",
    pch = point_types, col = point_cols
)
legend(
    "topleft", pch = c(16, 15, 17), col = c(rep_col, dem_col, squad_col),
    bty = "n", legend = c("Republicans", "Democrats", "The Squad")
)
par(opar)
invisible(dev.off())

## Panel (c): 2D W-NOMINATE Dim. 1 vs. GGUM
w_nom_ideo <- nominate2d$legislators$coord1D
pdf("plots/figE1c.pdf", width = 6, height = 4)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(
    x = w_nom_ideo, y = ggum_ideo,
    xlab = "2D W-NOMINATE Dim. 1", ylab = "GGUM Ideology",
    pch = point_types, col = point_cols
)
legend(
    "topleft", pch = c(16, 15, 17), col = c(rep_col, dem_col, squad_col),
    bty = "n", legend = c("Republicans", "Democrats", "The Squad")
)
par(opar)
invisible(dev.off())

## Panel (d): 2D W-NOMINATE Dim. 1 vs. GGUM, 115th Congress
MCs_used    <- rownames(nom2d115$legislators)
nom2d_dim1  <- setNames(nom2d115$legislators$coord1D, MCs_used)
ggum_theta  <- ggum_post$estimates$theta[match(MCs_used, icpsrs115)]
lc_names    <- "AMASH|GOSAR|GRIFFITH|JONES, W|MASSIE"
lc_ids      <- members115$icpsr[grepl(lc_names, members115$bioname)]
in_caucus   <- grepl(paste(lc_ids, collapse = "|"), icpsrs115)
member_idx  <- match(icpsrs115, members115$icpsr)
is_dem      <- members115$party_code[member_idx] == 100
point_types <- ifelse(in_caucus, 17, ifelse(is_dem, 15, 16))
lc_col      <- "black"
dem_col     <- "#00008080"
rep_col     <- "#80000080"
point_cols  <- ifelse(in_caucus, lc_col, ifelse(is_dem, dem_col, rep_col))
pdf("plots/figE1d.pdf", width = 6, height = 4)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(
    x = nom2d_dim1, y = ggum_theta,
    xlab = "2D W-NOMINATE Dim. 1", ylab = "GGUM Ideology",
    pch = point_types, col = point_cols
)
legend(
    "topleft", pch = c(16, 15, 17), col = c(rep_col, dem_col, lc_col),
    bty = "n", legend = c("Republicans", "Democrats", "Liberty Caucus")
)
par(opar)
invisible(dev.off())


##### Reproduce Figures E3-5 -----

## Define function to get NOMINATE response probabilities
utility <- function(x, z, beta, w) beta * exp(-0.5 * w * (x - z)^2)
nominate_irf <- function(nomObject, items) {
    if ( !inherits(nomObject, "nomObject") ) {
        stop("nomObject must be of class 'nomObject'")
    }
    if ( nomObject$dimensions > 1 ) {
        stop("nominate_irf currently only accommodates 1D models")
    }
    ths  <- seq(-1, 1, 0.01)
    beta <- nomObject$beta
    w    <- nomObject$weights
    y    <- with(nomObject$rollcalls, midpoint1D[item] - spread1D[item])
    n    <- with(nomObject$rollcalls, midpoint1D[item] + spread1D[item])
    return(plogis(utility(ths, y, beta, w) - utility(ths, n, beta, w)))
}

## Run 1D W-NOMINATE
nominate1d <- wnominate(rc, dims = 1, lop = 0.01, polarity = aic[1])
saveRDS(nominate1d, file = "output/H116-nominate1d.rds")
nom1d_dim1 <- nominate1d$legislators$coord1D

## Define axis labels, AOC index, & plot dims, which don't change between plots
xlabel <- expression(theta)
ylabel <- "Pr( Yea )"
aoc    <- which(rownames(responses) == "21949")
width  <- 5
height <- 3

## Plot Figure E3(b)
## (NOMINATE IRF for H.J. Res. 31 vote)
## (panel (a) is the same as fig9a)
ths   <- seq(-1, 1, 0.01)
item  <- which(colnames(responses) == "86")
probs <- nominate_irf(nominate1d, item)
pdf("plots/figE3b.pdf", height = height, width = width)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(nom1d_dim1[responses[ , item] == 1], side = 3)
rug(nom1d_dim1[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())

## Plot Figure E4(b)
## (NOMINATE IRF for H.R. 2740)
## (panel (a) is the same as fig9b)
item  <- which(colnames(responses) == "366")
probs <- nominate_irf(nominate1d, item)
pdf("plots/figE4b.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(nom1d_dim1[responses[ , item] == 1], side = 3)
rug(nom1d_dim1[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())

## Plot Figure E5(b)
## (NOMINATE IRF for H. Res. 326)
## (panel (a) is the same as figI3)
item  <- which(colnames(responses) == "651")
probs <- nominate_irf(nominate1d, item)
pdf("plots/figE5b.pdf", height = 3, width = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
trimplot(ths, y = probs, type = "l", xlab = xlabel, ylab = ylabel)
rug(nom1d_dim1[responses[ , item] == 1], side = 3)
rug(nom1d_dim1[responses[ , item] == 0], side = 1)
par(opar)
invisible(dev.off())


##### Reproduce Table E.1 -----
## Define some fit stat functions
brier <- function(obs, predprobs) mean((predprobs - obs)^2, na.rm = TRUE)
apre  <- function(responses, predictions) {
    minority_votes <- apply(responses, 2, function(v) min(table(v)))
    classification_errors <- colSums(responses != predictions, na.rm = TRUE)
    return( sum(minority_votes - classification_errors) / sum(minority_votes) )
}
## Create results storage object
fitstats <- data.frame(
    Model = c("GGUM", outer(c("1D", "2D"), c("CJR", "W-NOMINATE", "OC"), paste)),
    Correct = 0, APRE = 0, Brier = 0, AUC = 0
)
## Get GGUM fit statistics
J_mat <- matrix(1, nrow = nrow(responses), ncol = ncol(responses))
probs <- ggumProbability(J_mat, theta, alpha, delta, tau)
preds <- ifelse(probs > 0.5, 1, 0)
idx   <- which(fitstats$Model == "GGUM")
fitstats$Correct[idx] <- mean(preds == responses, na.rm = TRUE)
fitstats$APRE[idx] <- apre(responses, preds)
fitstats$Brier[idx] <- brier(responses, probs)
fitstats$AUC[idx] <- auc(roc(response = c(responses), predictor = c(preds)))
## Get CJR fit statistics
cjr_means <- cjr_posterior_summary$statistics[ , "Mean"]
cjr_theta <- cjr_means[grepl("theta", names(cjr_means))]
cjr_beta  <- cjr_means[grepl("beta",  names(cjr_means))]
cjr_alpha <- cjr_means[grepl("alpha", names(cjr_means))]
probs <- cjrProbability(J_mat, cjr_theta, cjr_alpha, cjr_beta)
preds <- ifelse(probs > 0.5, 1, 0)
idx   <- which(fitstats$Model == "1D CJR")
fitstats$Correct[idx] <- mean(preds == responses, na.rm = TRUE)
fitstats$APRE[idx] <- apre(responses, preds)
fitstats$Brier[idx] <- brier(responses, probs)
fitstats$AUC[idx] <- auc(roc(response = c(responses), predictor = c(preds)))
## Get 2D CJR fit statistics
cjr_means <- cjr2d_posterior_summary$statistics[ , "Mean"]
cjr_theta <- cbind(
    cjr_means[grepl("theta.+1$", names(cjr_means))],
    cjr_means[grepl("theta.+2$", names(cjr_means))]
)
cjr_beta1 <- cjr_means[grepl("beta.+1$", names(cjr_means))]
cjr_beta2 <- cjr_means[grepl("beta.+2$", names(cjr_means))]
j1 <- which(grepl("beta.43.2",  names(cjr_beta2)))
j2 <- which(grepl("beta.754.1", names(cjr_beta1)))
cjr_beta1 <- c( cjr_beta1[1:(j1-1)], 0, cjr_beta1[j1:(length(cjr_alpha)-1)])
cjr_beta2 <- c(cjr_beta2[1:(j2-1)], 0, cjr_beta2[j2:(length(cjr_alpha)-1)])
cjr_beta  <- cbind(cjr_beta1, cjr_beta2)
cjr_alpha <- cjr_means[grepl("alpha", names(cjr_means))]
probs <- cjrProbability(J_mat, cjr_theta, cjr_alpha, cjr_beta)
preds <- ifelse(probs > 0.5, 1, 0)
idx   <- which(fitstats$Model == "2D CJR")
fitstats$Correct[idx] <- mean(preds == responses, na.rm = TRUE)
fitstats$APRE[idx] <- apre(responses, preds)
fitstats$Brier[idx] <- brier(responses, probs)
fitstats$AUC[idx] <- auc(roc(response = c(responses), predictor = c(preds)))
## Get 1D W-NOMINATE fit statistics
idx   <- which(fitstats$Model == "1D W-NOMINATE")
fitstats$Correct[idx] <- nominate1d$fits["correctclass1D"] / 100
fitstats$APRE[idx] <- nominate1d$fits["apre1D"]
fitstats$Brier[idx] <- NA
fitstats$AUC[idx] <- NA
## Get 2D W-NOMINATE fit statistics
idx   <- which(fitstats$Model == "2D W-NOMINATE")
fitstats$Correct[idx] <- nominate2d$fits["correctclass2D"] / 100
fitstats$APRE[idx] <- nominate2d$fits["apre2D"]
fitstats$Brier[idx] <- NA
fitstats$AUC[idx] <- NA
## Get 1D OC fit statistics
aic1  <- as.character(with(members115, icpsr[which.max(nominate_dim1)]))
aic1  <- which(rownames(responses115) == aic1)
oc1d  <- oc(rc, dims = 1, lop = 0.01, polarity = aic1)
idx   <- which(fitstats$Model == "1D OC")
fitstats$Correct[idx] <- oc1d$fits[1]
fitstats$APRE[idx] <- oc1d$fits[2]
fitstats$Brier[idx] <- NA
fitstats$AUC[idx] <- NA
## Get 2D OC fit statistics
aic2  <- as.character(with(members115, icpsr[which.max(nominate_dim2)]))
aic2  <- which(rownames(responses115) == aic2)
aic   <- c(aic1, aic2)
oc2d  <- oc(rc, dims = 2, lop = 0.01, polarity = aic)
idx   <- which(fitstats$Model == "2D OC")
fitstats$Correct[idx] <- oc2d$fits[1]
fitstats$APRE[idx] <- oc2d$fits[2]
fitstats$Brier[idx] <- NA
fitstats$AUC[idx] <- NA
## View results
# print(fitstats, digits = 2)
#           Model Correct APRE Brier  AUC
# 1          GGUM    0.96 0.89 0.029 0.96
# 2        1D CJR    0.96 0.88 0.031 0.95
# 3        2D CJR    0.96 0.89 0.027 0.96
# 4 1D W-NOMINATE    0.96 0.88    NA   NA
# 5 2D W-NOMINATE    0.95 0.85    NA   NA
# 6         1D OC    0.97 0.91    NA   NA
# 7         2D OC    0.97 0.92    NA   NA
## Save results
write.csv(fitstats, file = "tables/TableE1.csv", row.names = FALSE)

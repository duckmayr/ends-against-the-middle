##### Setup -----
## Load required packages
library(bggum)     ## For GGUM-MC3 sampler and related tools
library(dplyr)     ## For data manipulation and summarization
library(tidyr)     ## For data reshaping
library(tibble)    ## For extra data manipulation functionality not in dplyr
## Source helper functions (documented in code/util.R)
source("code/00-util.R")
## Ensure required directories exist
prepare_directories()
## Compile functions to get log likelihoods:
Rcpp::sourceCpp("code/log_likelihood.cpp")


##### Compare log likelihood and theta s.d.s for the House -----
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

## Load the posterior samples and summaries
ggum_chains <- readRDS(file = "output/H116-chains.rds")
cjr_chain   <- readRDS(file = "output/cjr-house-chain.rds")
ggum_post   <- readRDS(file = "output/H116-ggum-post-summary.rds")
cjr_post    <- readRDS(file = "output/cjr-post-summary.rds")

## Get the log likelihood and mean theta SD for CJR
cjr_house_loglik  <- cjr_log_likelihood(responses, cjr_chain)
cjr_house_mean_sd <- mean(with(
    cjr_post, statistics[grepl("theta", rownames(statistics)), "SD"]
))
## Get the log likelihood and mean theta SD for GGUM
ggum_house_loglik  <- ggum_log_likelihood(responses, ggum_chains[[1]])
ggum_house_mean_sd <- mean(ggum_post$sds$theta_sds)


##### Compare log likelihood and theta s.d.s for the Court -----
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
# Next, the following are just some useful housekeeping variables
n <- nrow(response_matrix)
m <- ncol(response_matrix)
time_map <- mqData$time[match(ids, mqData$caseId)]
time_map <- time_map - min(time_map) + 1

## Load posterior samples
load("output/fold0-raw-results.RData")

## Get the log likelihood for and mean theta SD CJR
cjr_court_loglik  <- cjr_log_likelihood(response_matrix, cjr_chain1)
cjr_court_mean_sd <- mean(apply(cjr_chain1[ , 1:9], 2, sd))
## Get the log likelihood for and mean theta SD GGUM
ggum_court_loglik  <- ggum_log_likelihood(response_matrix, ggum_chain1)
ggum_court_mean_sd <- mean(apply(ggum_chain1[ , 1:9], 2, sd))
## Get the log likelihood for and mean theta SD M-Q
mq_court_loglik  <- mq_log_likelihood(response_matrix, mq_chain1, time_map)
mq_court_mean_sd <- mean(apply(mq_chain1[ , 1:54], 2, sd))


##### Reproduce Table 2 -----
results <- data.frame(Model = c("GGUM", "CJR", "MQ", "GGUM", "CJR"),
                      Data  = c(rep("Court", 3), rep("House", 2)),
                      check.names = FALSE)
results$`Log Likelihood` <- c(
    ggum_court_loglik, cjr_court_loglik, mq_court_loglik,
    ggum_house_loglik, cjr_house_loglik
)
N <- c(rep(sum(!is.na(response_matrix)), 3), rep(sum(!is.na(responses)), 2))
results$`L / N` <- sprintf("%0.2f", results$`Log Likelihood` / N)
results$`Log Likelihood` <- round(results$`Log Likelihood`)
results$`Mean Theta SD` <- sprintf("%0.2f", c(
    ggum_court_mean_sd, cjr_court_mean_sd, mq_court_mean_sd,
    ggum_house_mean_sd, cjr_house_mean_sd
))
# print(results)
write.csv(results, file = "tables/Table2.csv", row.names = FALSE)

##### Setup -----
## Load required packages
library(tibble)
library(dplyr)
library(tidyr)
library(bggum)
## Source helper functions (documented in code/util.R)
source("code/00-util.R")


##### Run MC3-GGUM & count non-monotonic items for each Congress -----
non_monotonicity <- expand.grid(Chamber = c("H", "S"), Congress = 110:115)
non_monotonicity$Proportion <- 0
for ( i in 110:115 ) {
    for ( j in c("H", "S") ) {
        ## Read in data on 
        votes <- read.csv(paste0("data/", j, i, "_votes.csv"))
        responses <- votes %>%
            filter(!grepl("^999", icpsr)) %>%
            mutate(response = case_when(
                cast_code == 1 ~ 1L,
                cast_code == 6 ~ 0L,
                TRUE ~ NA_integer_ 
            )) %>%
            select( icpsr, rollnumber, response) %>%
            pivot_wider(names_from = "rollnumber", values_from = "response") %>%
            column_to_rownames("icpsr") %>%
            as.matrix()
        ## Eliminate unanimous votes
        unanimous_votes <- which(apply(responses, 2, is_unanimous))
        if ( length(unanimous_votes) == 0 ) unanimous_votes <- 1e6
        responses <- responses[ , -unanimous_votes]
        ## Generate samples
        set.seed(1e6 * exp(1))
        chain <- ggumMC3(
            data = responses,
            sample_iterations = 100000,
            burn_iterations = 10000,
            temp_n_draws = 50,
            n_temps = 6,
            alpha_prior_params = c(1.5, 1.5, 0.25, 8.00)
        )
        ## Post-process
        members <- read.csv(paste0("data/", j, i, "_members.csv"))
        aic <- as.character(members$icpsr[which.min(members$nominate_dim1)])
        aic <- which(rownames(responses) == aic)
        processed_ggum_chain <- post_process(chain, aic, "-")
        ## Save posterior summary
        ggum_post <- summary(processed_ggum_chain)
        outfile <- paste0("output/", j, i, "-ggum-post-summary.rds")
        saveRDS(ggum_post, file = outfile)
        ## Check amount of non-monotonic items
        idx <- with(non_monotonicity, which(Chamber == j & Congress == i))
        non_monotonicity$Proportion[idx] <- prop_nonmonotonic(ggum_post)
    }
}


##### Reproduce Figure J.1 -----
## Setup the plotting device:
pdf("plots/figJ1.pdf", height = 4, width = 6)
## Edit the margins
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
## Plot the comparison
point_types <- ifelse(non_monotonicity$Chamber == "H", 19, 17)
point_cols  <- ifelse(non_monotonicity$Chamber == "H", "#d55e00", "#0072b2")
trimplot(
    x = non_monotonicity$Congress, y = non_monotonicity$Proportion,
    xlab = "Congress", ylab = "Proportion non-monotonic",
    pch = point_types, col = point_cols
)
## Add a legend
legend(
    "bottomright", pch = c(19, 17), col = c("#d55e00", "#0072b2"),
    bty = "n", legend = c("House", "Senate")
)
## Reset the margins
par(opar)
## And close the plotting device
dev.off()

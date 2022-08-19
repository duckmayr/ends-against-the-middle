##### Setup -----
## Load required packages
library(ggplot2) ## For plotting
library(bggum)   ## For GGUM analysis
library(ltm)     ## For GRM analysis
library(coda)    ## To assess convergence


##### Prepare data -----
## Read in the raw responses
responses <- read.csv("data/lucid_data.csv")
## It will be useful to record the question wording for later
questions <- responses[1, ]
QUESTIONS <- gsub("\n", " ", questions[ , grepl("^IMM", colnames(questions))])
## Eliminate respondents who do not pass the attention checks
feelings  <- which(names(responses) == "SCREENER_FEELINGS")
correct1  <- responses[ , feelings] == "Proud,None of the above"
interest  <- which(names(responses) == "SCREENER_INTEREST")
correct2  <- responses[ , interest] == "Extremely interested,Not interested at all"
colors    <- which(names(responses) == "SCREENER_COLORS")
correct3  <- responses[ , colors] == "Red,Green"
attentive <- correct1 & correct2 & correct3
responses <- responses[attentive, ]
## Get the ideology and party ID variables to validate
vinfo     <- responses[ , c(27:30, 32)]
party_id  <- (vinfo[ , "PID"] == "Republican") - (vinfo[ , "PID"] == "Democrat")
party_id  <- party_id + (vinfo[ , "R_STRENGTH"] == "A strong Republican")
party_id  <- party_id - (vinfo[ , "D_STRENGTH"] == "A strong Democrat")
party_id  <- party_id + (vinfo[ , "LEANERS"] == "Republican")
party_id  <- party_id - (vinfo[ , "LEANERS"] == "Democratic")
ideo_levs <- c("Very liberal", "Somewhat liberal", "Slightly liberal",
               "Moderate; middle of the road", "Slightly conservative",
               "Somewhat conservative", "Very conservative")
ideo      <- as.integer(factor(vinfo[ , "IDEO"], levels = ideo_levs))
vinfo     <- data.frame(party_id = party_id, ideo = ideo)
## Get just the responses to the immigration battery
responses <- responses[ , grepl("^IMM", colnames(responses))]
## Code the responses in {NA, 0, 1, 2, 3, 4}
q_options <- c(paste(c("Strongly", "Somewhat"), "disagree"),
               "Neither disagree nor agree",
               paste(c("Strongly", "Somewhat"), "agree"))
responses <- apply(responses, 2, match, q_options) - 1
## Eliminate respondents who straight-line
str8lines <- apply(responses, 1, function(x) all(x %in% 0:2) | all(x %in% 2:4))
responses <- responses[!str8lines, ]
vinfo     <- vinfo[!str8lines, ]
## Ensure all options have been chosen by some respondents
apply(responses, 2, table)


##### Sampling -----
## Tune hyperparameters
set.seed(123)
sds <- tune_proposals(responses, 5000)
sapply(sds, mean)
temps <- tune_temperatures(responses, 6, proposal_sds = sds)
temps
save(sds, temps, file = "output/lucid-hypers.RData")
## Draw posterior samples
chain1 <- ggumMC3(responses, proposal_sds = sds, temps = temps)
## Save output
saveRDS(chain1, "output/lucid-chain1.rds")
## Use "All undocumented immigrants currently living in the U.S. should be
## required to return to their home country" to post-process results
n <- nrow(responses)
m <- ncol(responses)
processed_chain <- post_process(chain1, n + m + 1, "+")
posterior_summary <- summary(processed_chain)
saveRDS(posterior_summary, "output/lucid-posterior.rds")
theta <- posterior_summary$estimates$theta
alpha <- posterior_summary$estimates$alpha
delta <- posterior_summary$estimates$delta
tau   <- posterior_summary$estimates$tau


##### Assess convergence -----
set.seed(456)
chain2 <- ggumMC3(responses, proposal_sds = sds, temps = temps)
processed_chain2 <- post_process(chain2, n + m + 1, "+")
convergence <- gelman.diag(list(processed_chain, processed_chain2))
summary(convergence$psrf[ , 1])


##### Compare with GRM -----
## Fit the GRM model and get the ability scores
grm_fit  <- grm(responses + 1)
grm_ideo <- factor.scores(grm_fit, resp.patterns = responses)
## Look at correlation between GGUM theta, GRM theta, and self ideo
## (referenced in Section 5.1)
cor(
    cbind(
        self_ideo = vinfo$ideo,
        grm_ideo = grm_ideo$score.dat$z1,
        ggum_ideo = theta
    ),
    use = "pairwise"
)


##### Reproduce Figure 6 -----
irf_probs <- function(alpha, delta, tau, theta_range = c(-2.75, 2.5)) {
    theta <- seq(from = theta_range[1], to = theta_range[2], by = 0.01)
    K     <- length(tau)
    n     <- length(theta)
    res   <- matrix(NA_real_, nrow = length(theta), ncol = K)
    for ( k in 1:K ) {
        res[ , k] <- ggumProbability(rep(k-1, n), theta, alpha, delta, tau)
    }
    return(res)
}
response_options <- c(
    "Strongly Disagree", "Disagree", "Neither", "Agree", "Strongly Agree"
)
option_colors <- colorRampPalette(colors = c("#808080", "black"))(5)
grm_item9  <- plot(grm_fit, zrange = c(-2.25, 2.25), items = 9)
grm_item2  <- plot(grm_fit, zrange = c(-2.25, 2.25), items = 2)
grm_item4  <- plot(grm_fit, zrange = c(-2.25, 2.25), items = 4)
theta_plot <- seq(-2.75, 2.5, 0.01)
ggum_item9 <- irf_probs(alpha[9], delta[9], tau[[9]])
ggum_item2 <- irf_probs(alpha[2], delta[2], tau[[2]])
ggum_item4 <- irf_probs(alpha[4], delta[4], tau[[4]])
## Plot Item 2:
dat <- data.frame(
    Model  = c(rep("GGUM", n1 * 5), rep("GRM", n2 * 5)),
    Option = factor(c(rep(opts, each = n1), rep(opts, each = n2)),
                    levels = opts),
    theta  = c(rep(theta_plot, 5), rep(grm_item2$z, 5)),
    p      = c(c(ggum_item2), c(grm_item2$pr$IMM_2))
)
map <- aes(x = theta, y = p, color = Option, linetype = Option, group = Option)
plt <- ggplot(data = dat, mapping = map) +
    geom_line() +
    facet_wrap(~ Model, scales = "free_x") +
    scale_color_manual(values = option_colors) +
    scale_linetype_manual(values = 5:1) +
    xlab(expression(theta)) +
    ylab("Probability") +
    theme_bw() + 
    theme(
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.box.spacing = unit(0, "in")
    )
ggsave(
    filename = "plots/fig6a.tif", plot = plt, device = "tiff",
    width = 6.5, height = 2.5, units = "in",
    dpi = 1000, compression = "lzw"
)
## Plot Item 4:
dat <- data.frame(
    Model  = c(rep("GGUM", n1 * 5), rep("GRM", n2 * 5)),
    Option = factor(c(rep(opts, each = n1), rep(opts, each = n2)),
                    levels = opts),
    theta  = c(rep(theta_plot, 5), rep(grm_item4$z, 5)),
    p      = c(c(ggum_item4), c(grm_item4$pr$IMM_4))
)
map <- aes(x = theta, y = p, color = Option, linetype = Option, group = Option)
plt <- ggplot(data = dat, mapping = map) +
    geom_line() +
    facet_wrap(~ Model, scales = "free_x") +
    scale_color_manual(values = option_colors) +
    scale_linetype_manual(values = 5:1) +
    xlab(expression(theta)) +
    ylab("Probability") +
    theme_bw() + 
    theme(
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.box.spacing = unit(0, "in")
    )
ggsave(
    filename = "plots/fig6b.tif", plot = plt, device = "tiff",
    width = 6.5, height = 2.5, units = "in",
    dpi = 1000, compression = "lzw"
)
## Plot Item 9:
n1 <- length(theta_plot)
n2 <- length(grm_item9$z)
opts <- c("Strongly disagree", "Disagree", "Neither", "Agree", "Strongly agree")
dat <- data.frame(
    Model  = c(rep("GGUM", n1 * 5), rep("GRM", n2 * 5)),
    Option = factor(c(rep(opts, each = n1), rep(opts, each = n2)),
                    levels = opts),
    theta  = c(rep(theta_plot, 5), rep(grm_item9$z, 5)),
    p      = c(c(ggum_item9), c(grm_item9$pr$IMM_9))
)
map <- aes(x = theta, y = p, color = Option, linetype = Option, group = Option)
plt <- ggplot(data = dat, mapping = map) +
    geom_line() +
    facet_wrap(~ Model, scales = "free_x") +
    scale_color_manual(values = option_colors) +
    scale_linetype_manual(values = 5:1) +
    xlab(expression(theta)) +
    ylab("Probability") +
    theme_bw() + 
    theme(
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.box.spacing = unit(0, "in")
    )
ggsave(
    filename = "plots/fig6c.tif", plot = plt, device = "tiff",
    width = 6.5, height = 2.5, units = "in",
    dpi = 1000, compression = "lzw"
)


##### Fit statistics (mentioned in text) -----
ggum_preds <- array(0, dim = c(n, m, 5))
for ( k in 0:4 ) {
    tmp <- matrix(k, nrow = n, ncol = m)
    ggum_preds[ , , k+1] <- ggumProbability(tmp, theta, alpha, delta, tau)
}
grm_preds  <- fitted(
    grm_fit, resp.patterns = responses + 1,
    type = "conditional-probabilities"
)
## Need to convert from list to array
list_to_array <- function(list_in) {
    dims <- c(nrow(list_in[[1]]), ncol(list_in[[1]]), length(list_in))
    return(array(as.numeric(unlist(list_in)), dim = dims))
}
grm_preds  <- aperm(list_to_array(grm_preds), c(1, 3, 2))
ggum_classifications <- apply(ggum_preds, 1:2, which.max)
grm_classifications  <- apply(grm_preds, 1:2, which.max)
actual <- factor(c(responses + 1))
# caret::confusionMatrix(factor(c(grm_classifications)),  actual)
# caret::confusionMatrix(factor(c(ggum_classifications)), actual)

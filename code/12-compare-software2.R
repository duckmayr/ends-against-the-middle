##### Setup -----
library(bggum)

##### Simulate Data -----
## Define function to generate draws from the four-parameter Beta distribution
r4beta = function(n, shape1, shape2, a, b){
    return((b - a) * rbeta(n, shape1, shape2) + a)
}
## Set the number of respondents and items
n = 1000
m = 10
## Draw item and respondent parameters
set.seed(1)
alphas = r4beta(m, 1.5, 1.5, 0.4, 2.0)
deltas = r4beta(m, 2.0, 2.0, -3.0, 3.0)
taus = sapply(1:m, function(j) {
    return(c(0, r4beta(3, 2.0, 2.0, -2.0, 0.0)))
}, simplify = FALSE)
thetas = rnorm(n)
## Simulate responses
simulate_response = function(theta, alpha, delta, tau) {
    probs = ggumProbability(1:length(tau) - 1, theta, alpha, delta, tau)
    return(which(rmultinom(1, 1, probs) == 1) - 1)
}
simulate_response_matrix = function(thetas, alphas, deltas, taus) {
    return(t(sapply(1:length(thetas), function(i) {
        sapply(1:length(alphas), function(j) {
            simulate_response(thetas[i], alphas[j], deltas[j], taus[[j]])
        })
    })))
}
responses = simulate_response_matrix(thetas, alphas, deltas, taus)
## Order the items (to give the Windows program every benefit of the doubt)
ordering  = order(deltas)
alphas    = alphas[ordering]
deltas    = deltas[ordering]
taus      = taus[ordering]
responses = responses[ , ordering]
## Write out the data in a format the Windows program can use
collapse = function(x, sep = "") paste(x, collapse = sep)
dat = sapply(1:n, function(i) collapse(sprintf("  %d", responses[i, ])))
writeLines(dat, con = "data/sim-data.txt")

##### Generate MC3-GGUM posterior samples -----
set.seed(2)
sds = tune_proposals(responses, 5000, thetas = thetas,
                     alphas = alphas, deltas = deltas, taus = taus)
iters = 1000
delta_init1 = c(-2.73, -1.03, -0.08, -0.03, 0.65, 0.98, 1.33, 1.51, 1.76, 2.21)
delta_init2 = rep(4.5, 10)
temps = c(1, 0.95, 0.9, 0.86, 0.82, 0.78)
set.seed(3)
chain1 = ggumMC3(
    data = responses,
    sample_iterations = iters,
    burn_iterations = 0,
    proposal_sds = sds,
    delta_init = delta_init1,
    temps = temps,
    swap_interval = 1
)
saveRDS(chain1, "output/winsim-bggum-chain1.rds")
set.seed(3)
chain2 = ggumMC3(
    data = responses,
    sample_iterations = iters,
    burn_iterations = 0,
    proposal_sds = sds,
    delta_init = delta_init2,
    temps = temps,
    swap_interval = 1
)
saveRDS(chain2, "output/winsim-bggum-chain2.rds")

##### Reproduce Figure C.1 -----
chain1 = read.table("data/sim01-item-draws.txt")
chain2 = read.table("data/sim03-item-draws.txt")
idx    = seq(1, nrow(chain1), floor(nrow(chain1) / 1000))
delta0 = c(-2.73, -1.03, -0.08, -0.03, 0.65, 0.98, 1.33, 1.51, 1.76, 2.21)
delta1 = apply(chain1[ , 11:20], 2, mean)
delta2 = apply(chain2[ , 11:20], 2, mean)
pdf("plots/figC1.pdf", height = 6, width = 6)
opar <- par(mar = c(3, 3, 1, 1) + 0.1)
plot(delta0, delta1, xlim = c(-4, 4), ylim = c(-4, 4), pch = 19, xaxt = "n",
     yaxt = "n", col = "#80000080", xlab = "", ylab = "")
points(delta0, delta2, pch = 19, col = "#00008080")
axis(1, at = seq(-4, 4, 2), tick = FALSE, line = -0.75)
axis(2, at = seq(-4, 4, 2), tick = FALSE, line = -0.75)
mtext("Estimates", side = 2, line = 1.5)
mtext("True Values", side = 1, line = 1.5)
legend("topleft", pch = 19, col = c("#80000080", "#00008080"), bty = "n",
       legend = c("Given True Ordering", "Given Alternate Ordering"))
par(opar)
dev.off()

##### Reproduce Figure C.2 -----
pal  = c("#0072b280", "#d55e0080")
one  = m + 1
two  = m + 2
xlab = expression(paste(delta[1], " Draws"))
ylab = expression(paste(delta[2], " Draws"))
xat  = c(-4, -2.73, -2, 0, 2, 4)
yat  = c(-4, -2, -1.03, 0, 2, 4)
lgnd = c(
    expression(paste("Chain 1 (", delta^0, " = truth)")),
    expression(paste("Chain 2 (", delta^0, " = ", 4.5, ")"))
)
xticklabs = c("-4", expression(delta[1]), "-2", "0", "2", "4")
yticklabs = c("-4", "-2", expression(delta[2]), "0", "2", "4")

chain1 = readRDS("output/winsim-bggum-chain1.rds")
chain2 = readRDS("output/winsim-bggum-chain2.rds")
pdf("plots/figC2a.pdf", height = 11, width = 11)
idx  = seq(from = 1, to = nrow(chain1), by = floor(nrow(chain1) / 1000))
opar = par(mar = c(7, 7, 2, 2) + 0.1)
plot(
    x = chain1[idx, n + one], y = chain1[idx, n + two],
    xlim = c(-5, 5), ylim = c(-5, 5), xlab = "", ylab = "",
    col = pal[1], pch = 20,
    cex.lab = 3, cex = 2, mgp = c(4, 2, 0),
    yaxt = "n", xaxt = "n"
)
points(
    x = chain2[idx, n + one], y = chain2[idx, n + two],
    col = pal[2], pch = 20, type = "o", cex = 2
)
abline(v = -2.73, lty = 2, col = "#80808080")
abline(h = -1.03, lty = 2, col = "#80808080")
axis(side = 1, at = xat, cex.axis = 3, mgp = c(4.5, 2, 0), labels = xticklabs)
axis(side = 2, at = yat, cex.axis = 3, mgp = c(4, 1.5, 0), labels = yticklabs)
mtext(text = xlab, side = 1, cex = 3, line = 4.75)
mtext(text = ylab, side = 2, cex = 3, line = 4)
legend(-2.7, 5.5, pch = 20, col = pal, bty = "n", legend = lgnd, cex = 3)
par(opar)
dev.off()

chain1 = read.table("data/sim01-item-draws.txt")
chain2 = read.table("data/sim02-item-draws.txt")
pdf("plots/figC2b.pdf", height = 11, width = 11)
idx  = seq(from = 1, to = nrow(chain1), by = floor(nrow(chain1) / 1000))
opar = par(mar = c(7, 7, 2, 2) + 0.1)
plot(
    x = chain1[idx, one], y = chain1[idx, two],
    xlim = c(-5, 5), ylim = c(-5, 5), xlab = "", ylab = "",
    col = pal[1], pch = 20,
    cex.lab = 3, cex = 2, mgp = c(4, 2, 0),
    yaxt = "n", xaxt = "n"
)
points(
    x = chain2[idx, one], y = chain2[idx, two],
    col = pal[2], pch = 20, type = "o", cex = 2
)
abline(v = -2.73, lty = 2, col = "#80808080")
abline(h = -1.03, lty = 2, col = "#80808080")
axis(side = 1, at = xat, cex.axis = 3, mgp = c(4.5, 2, 0), labels = xticklabs)
axis(side = 2, at = yat, cex.axis = 3, mgp = c(4, 1.5, 0), labels = yticklabs)
mtext(text = xlab, side = 1, cex = 3, line = 4.75)
mtext(text = ylab, side = 2, cex = 3, line = 4)
legend(-2.7, 5.5, pch = 20, col = pal, bty = "n", legend = lgnd, cex = 3)
par(opar)
dev.off()

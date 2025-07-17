library(MASS)
library(microbenchmark)
set.seed(0625)
# Design
n_vals <- c(200, 500, 1000)
k_vals <- c(2, 5, 10)
rho <- 0.65
rs_target <- (6 / pi) * asin(rho / 2)
epsilon <- 0.002
m_flex <- 20
reps <- 1000  # increase for full run

# Helpers
make_cs <- function(k, val) {
  mat <- matrix(val, k, k)
  diag(mat) <- 1
  mat
}
bias_metric <- function(cormat, target) max(abs(cormat - target))

# Output container
summary_table <- data.frame()

# Loop
for (n in n_vals) {
  for (k in k_vals) {
    target <- make_cs(k, rs_target)
    
    ic_deltas <- numeric(reps)
    flex_deltas <- numeric(reps)
    ic_times <- numeric(reps)
    flex_times <- numeric(reps)
    flex_redraws <- integer(reps)
    
    for (rep in 1:reps) {
      x <- matrix(rnorm(n * k), nrow = n)
      Sigma <- make_cs(k, rs_target)
      
      # IC timing and run
      ic_benchmark <- microbenchmark(
        {
          z_ic <- mvrnorm(n, rep(0, k), Sigma)
          rank_z_ic <- apply(z_ic, 2, rank, ties.method = "first")
          x_ic <- matrix(NA, n, k)
          for (j in 1:k) x_ic[ , j] <- sort(x[ , j])[rank_z_ic[ , j]]
          ic_deltas[rep] <- bias_metric(cor(x_ic, method = "spearman"), target)
        },
        times = 1
      )
      ic_times[rep] <- ic_benchmark$time / 1e6  # convert ns → ms
      
      # flexIC timing and run
      flex_benchmark <- microbenchmark(
        {
          best_bias <- Inf
          flex_redraw <- 0
          for (m in 1:m_flex) {
            z <- mvrnorm(n, rep(0, k), Sigma)
            rank_z <- apply(z, 2, rank, ties.method = "first")
            x_flex <- matrix(NA, n, k)
            for (j in 1:k) x_flex[ , j] <- sort(x[ , j])[rank_z[ , j]]
            this_bias <- bias_metric(cor(x_flex, method = "spearman"), target)
            flex_redraw <- m
            if (this_bias < best_bias) best_bias <- this_bias
            if (best_bias <= epsilon) break
          }
          flex_deltas[rep] <- best_bias
          flex_redraws[rep] <- flex_redraw
        },
        times = 1
      )
      flex_times[rep] <- flex_benchmark$time / 1e6  # convert ns → ms
    }
    
    # Summary row
    summary_table <- rbind(summary_table, data.frame(
      n = n,
      k = k,
      epsilon = epsilon,
      IC.ms = formatC(median(ic_times), format = "f", digits = 2),
      FLEX.ms = formatC(median(flex_times), format = "f", digits = 2),
      Time.Difference.ms = formatC(median(flex_times - ic_times), format = "f", digits = 2),
      redraws = as.integer(median(flex_redraws)),
      IC.max.Delta = round(mean(ic_deltas), 4),
      FLEX.max.Delta = round(mean(flex_deltas), 4),
      Bias.Reduction.pct = round(100 * (mean(ic_deltas) - mean(flex_deltas)) / mean(ic_deltas), 1)
    ))
  }
}

# Output
print(summary_table, row.names = FALSE)

## Sim results analysis

devtools::load_all()
library(dplyr)
library(ggplot2)

#### Simulation 1 ----

# Standard normal error, single beta
sim_1 <- readRDS("inst/sim_res/2024-03-26_beta_1.rds")

sim_1_df <- sim_res_to_df(sim_1)
mean(sim_1_df$pct_truncated)
mean(sim_1_df$pct_censored)

mean(sim_1_df$naive_beta)
sd(sim_1_df$naive_beta)
mean(sim_1_df$beta)
sd(sim_1_df$beta)

sim_1_df %>%
  ggplot() +
  aes(x = beta) +
  geom_density() +
  theme_bw() +
  xlab("Beta Estimate") +
  ylab("Probability Density") +
  geom_vline(
    xintercept = mean(sim_1_df$beta),
    color = "darkgreen"
  ) +
  geom_vline(
    xintercept = 1,
    lty = "dashed",
    color = "darkblue"
  ) +
  geom_vline(
    xintercept = mean(sim_1_df$naive_beta),
    lty = "dashed",
    color = "red"
  ) +
  annotate(
    geom = "segment",
    x = 0.5, y = 4,
    xend = .81, yend = 3
  ) +
  annotate(
    geom = "label",
    x = 0.5, y = 4,
    label = "Naive Beta Average\n0.819"
  ) +
  annotate(
    geom = "segment",
    x = 1.3, y = 4,
    xend = 1.02, yend = 3
  ) +
  annotate(
    geom = "label",
    x = 1.3, y = 4,
    label = "Our Estimate Average\n1.011"
  )

mean(sim_1_df$norm_score)
median(sim_1_df$norm_score)

sim_1_df %>%
  filter(norm_score < 1000) %>%
  ggplot() +
  aes(x = beta, y = beta_score) +
  geom_point() +
  theme_bw() +
  xlab("Beta Estimate") +
  ylab("Beta Component of Score Vector")

## Plots of the hazard

x <- seq(-3, 3, by = 0.01)
plot(x, rep(0, length(x)), type = "l", ylim = c(0, 10))
for (i in 1:100) {
  y <- predict_hazard(ltrc_obj = sim_1[[i]], vals = x)
  lines(x, y, col = "red")
}
lines(x, true_hazard(x), col = "black", lty = "dashed")

## get avg hazard:
haz_mat <- matrix(NA, nrow = length(sim_1), ncol = length(x))
for (i in 1:length(sim_1)) {
  y <- predict_hazard(ltrc_obj = sim_1[[i]], vals = x, within_bounds = TRUE)
  haz_mat[i,] <- y
}
new_y <- colSums(haz_mat)

lines(x, y, col = "blue")

ggplot() +
  geom_line(
    aes(x = x, y = true_hazard(x)),
    color = "black",
    lty = "dashed"
  ) +
  geom_line(
    aes(x = x, y = y),
    color = "blue"
  ) +
  theme_bw() +
  xlab("Epsilon") +
  ylab("Hazard Function") +
  ggtitle("Average Estimated Hazard Function", "Estimated (blue) vs. Truth (black)")



#### Simulation 2 ----

# Standard normal error, single beta
sim_2 <- readRDS("inst/sim_res/2024-03-27_beta_2_neg1.rds")

sim_2_df <- sim_res_to_df(sim_2)

mean(sim_2_df$naive_beta_1)
sd(sim_2_df$naive_beta_1)
mean(sim_2_df$beta_1)
sd(sim_2_df$beta_1)

mean(sim_2_df$naive_beta_2)
sd(sim_2_df$naive_beta_2)
mean(sim_2_df$beta_2)
sd(sim_2_df$beta_2)

mean(sim_2_df$norm_score)
median(sim_2_df$norm_score)


#### Simulation 3 ----
sim_3 <- readRDS("inst/sim_res/2024-03-28_beta_1.rds")
sim_3_df <- sim_res_to_df(sim_3)

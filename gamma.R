# Gamma factor
rm(list=ls())
library(tidyverse)
library(tidyquant)
library(scales)

sample <- read_csv("sample data 2.csv", 
        col_names = FALSE,
        col_types = list(
          X1 = col_character(),
          X2 = col_date(format = "%d-%b-%y"),
          X3 = col_double()
          )
        )
colnames(sample) <- c("ticker", "date", "price")

# pretty chart
sample |>
  ggplot(aes(
    x = date,
    y = price,
    color = ticker
  )) +
  geom_line() +
  labs(
    x = NULL,
    y = NULL,
    title = "Stock prices of index constituents"
  ) +
  theme(legend.position = "none")

all_returns <- sample |>
  group_by(ticker) |>
  mutate(ret = price / lag(price) - 1) |>
  select(ticker, date, ret) |>
  drop_na(ret)

#summary stats
summary_returns <- all_returns |>
  group_by(ticker) |>
  summarize(across(
    ret,
    list(
      daily_mean = mean,
      daily_sd = sd,
      daily_min = min,
      daily_max = max
    ),
    .names = "{.fn}"
  )) |>
  print(n = Inf)

# only use tickers with complete data
sample <- sample |>
  group_by(ticker) |>
  mutate(n = n()) |>
  ungroup() |>
  filter(n == max(n)) |>
  select(-n)

# monthly returns
returns <- sample |>
  mutate(month = floor_date(date, "month")) |>
  group_by(ticker, month) |>
  summarize(price = last(price), .groups = "drop_last") |>
  mutate(ret = price / lag(price) - 1) |>
  drop_na(ret) |>
  select(-price)

returns_matrix <- returns |>
  pivot_wider(
    names_from = ticker,
    values_from = ret
  ) |>
  select(-month)

Sigma <- cov(returns_matrix)
mu <- colMeans(returns_matrix)

N <- ncol(returns_matrix)
iota <- rep(1, N)
mvp_weights <- solve(Sigma) %*% iota
mvp_weights <- mvp_weights / sum(mvp_weights)

tibble(
  average_ret = as.numeric(t(mvp_weights) %*% mu),
  volatility = as.numeric(sqrt(t(mvp_weights) %*% Sigma %*% mvp_weights))
)

benchmark_multiple <- 30         # only bc bm return ~1%
mu_bar <- benchmark_multiple * t(mvp_weights) %*% mu
C <- as.numeric(t(iota) %*% solve(Sigma) %*% iota)
D <- as.numeric(t(iota) %*% solve(Sigma) %*% mu)
E <- as.numeric(t(mu) %*% solve(Sigma) %*% mu)
lambda_tilde <- as.numeric(2 * (mu_bar - D / C) / (E - D^2 / C))
efp_weights <- mvp_weights +
  lambda_tilde / 2 * (solve(Sigma) %*% mu - D * mvp_weights)

length_year <- 12   # Diff from online
a <- seq(from = -0.4, to = 1.9, by = 0.01)
res <- tibble(
  a = a,
  mu = NA,
  sd = NA
)

for (i in seq_along(a)) {
  w <- (1 - a[i]) * mvp_weights + (a[i]) * efp_weights
  res$mu[i] <- length_year * t(w) %*% mu   
  res$sd[i] <- sqrt(length_year) * sqrt(t(w) %*% Sigma %*% w)
}

res |>
  ggplot(aes(x = sd, y = mu)) +
  geom_point() +
  geom_point(
    data = res |> filter(a %in% c(0, 1)),
    size = 4
  ) +
  geom_point(
    data = tibble(
      mu = length_year * mu,       
      sd = sqrt(length_year) * sqrt(diag(Sigma))
    ),
    aes(y = mu, x = sd), size = 1
  ) +
  labs(
    x = "Annualized standard deviation",
    y = "Annualized expected return",
    title = "Efficient frontier for group constituents"
  ) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent)

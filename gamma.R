# Gamma factor
library(tidyverse)
library(tidyquant)
library(scales)

# download the exchange listings from 
# https://www.nasdaq.com/market-activity/stocks/screener?exchange=nyse
# and change the name
# NYSE <- tq_exchange("NYSE") doesn't work
# AMEX <- jsonlite::fromJSON('https://api.nasdaq.com/api/screener/stocks?tableonly=true&limit=25&offset=0&exchange=AMEX&download=true')
# head(AMEX$data$rows) doesn't work either

nyse <- read_csv("nyse.csv", col_names = TRUE)
amex <- read_csv("amex.csv", col_names = TRUE)
nasdaq <- read_csv("nasdaq.csv", col_names = TRUE)

nyse$Exchange <- "NYSE"
amex$Exchange <- "AMEX"
nasdaq$Exchange <- "NASDAQ"

tickers <- bind_rows(nyse, nasdaq, amex)

zip <- tickers |> 
  
index_prices <- tq_get(ticker,
                       get = "stock.prices",
                       from = "2000-01-01",
                       to = "2022-12-31"
)

# Warning message:
# There was 1 warning in `dplyr::mutate()`.
# ℹ In argument: `data.. = purrr::map(...)`.
# Caused by warning:
#  ! x = '-', get = 'stock.prices': Error in getSymbols.yahoo(Symbols = "-", env = <environment>, verbose = FALSE, : Unable to import “-”.
#  HTTP error 404.
#  Removing -.

index_prices |>
  ggplot(aes(
    x = date,
    y = adjusted,
    color = symbol
  )) +
  geom_line() +
  labs(
    x = NULL,
    y = NULL,
    color = NULL,
    title = "Stock prices of DOW index constituents"
  ) +
  theme(legend.position = "none")

all_returns <- index_prices |>
  group_by(symbol) |>
  mutate(ret = adjusted / lag(adjusted) - 1) |>
  select(symbol, date, ret) |>
  drop_na(ret)

all_returns |>
  group_by(symbol) |>
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

#   ticker <- tq_index("SP500")

volume <- index_prices |>
  group_by(date) |>
  summarize(volume = sum(volume * close))

volume |>
  ggplot(aes(x = date, y = volume)) +
  geom_line() +
  labs(
    x = NULL, y = NULL,
    title = "Aggregate daily trading volume of DOW index constitutens"
  ) +
  scale_y_continuous(labels = unit_format(unit = "B", scale = 1e-9))

volume |>
  ggplot(aes(x = lag(volume), y = volume)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1),
              linetype = "dashed"
  ) +
  labs(
    x = "Previous day aggregate trading volume",
    y = "Aggregate trading volume",
    title = "Persistence in daily trading volume of DOW index constituents"
  ) + 
  scale_x_continuous(labels = unit_format(unit = "B", scale = 1e-9)) +
  scale_y_continuous(labels = unit_format(unit = "B", scale = 1e-9))

index_prices <- index_prices |>
  group_by(symbol) |>
  mutate(n = n()) |>
  ungroup() |>
  filter(n == max(n)) |>
  select(-n)

returns <- index_prices |>
  mutate(month = floor_date(date, "month")) |>
  group_by(symbol, month) |>
  summarize(price = last(adjusted), .groups = "drop_last") |>
  mutate(ret = price / lag(price) - 1) |>
  drop_na(ret) |>
  select(-price)

returns_matrix <- returns |>
  pivot_wider(
    names_from = symbol,
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

benchmark_multiple <- 3
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
    title = "Efficient frontier for DOW index constituents"
  ) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent)
# Create simulated data

set.seed(123)

# simulate initial mass of litter
i <- rnorm(609, 4111, 156)
# log initial mass
log_i <- log(i)

# time in years
t <- rep(seq(0, 0.7, .035), 29)

# create clusters in data
group <- rep(1:29, each = 21)

# negative exponential
k <- rnorm(609, 1.3, .4)
log_mean_k <- log_i - (k*t)
ne <- data.frame(t = t, log_i = log_i, log_mean = log_mean_k, group_id = group)
save(ne, file = "data/sim_neg_exp.rda")

# weibull
a <- rnorm(609, 0.6, .1)
b <- rnorm(609, 1.2, .2)
log_mean_w <- log_i - ((t/b)^a)
w <- data.frame(t = t, log_i = log_i, log_mean = log_mean_w, group_id = group)
save(w, file = "data/sim_weibull.rda")

source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra)
library(speedglm)
library(foreach); library(doParallel)
library(robustbase)

# Specify data-generating distribution ------------------------------------

# Set of parameters 1
# lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# beta <- 0.25
# gamma <- 0.125
# kappa <- 1 / (2 * (gamma + 1 - beta))
# Set of parameters 2
lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
kappa <- 5 / 4
beta <- 2 - kappa
gamma <- 1 - kappa / 2
# Set of parameters 3
# lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
# beta <- 7/8; gamma <- 1/16
true_rate <- 1 / (2 * (gamma + 1 - beta))

# Define Estimation tasks
eta <- 2
ns <- sort(10^seq(from = 3, to = 7.5, by = 0.5), decreasing = T)
candidate_rates <- c(0.8 * true_rate, true_rate, 1.2 * true_rate)
estimation_tasks <- expand.grid(n = ns, candidate_rate = candidate_rates)

# Subsamples' tree
nb_children <- c(10, 5, 2, rep(1,6)) # nb of children at each level (at a given depth, each node has the
#same number of children)
phi <- function(i){ # indices of level i are [phi(i-1) + 1, phi(i)]
  if(i == 0){
    return(0) 
  }else if(i == 1){
    return(1)
  }else{
    phi(i-1) + prod(nb_children[1:(i-1)])
  }
}
levels.upper_bounds <- sapply(1:(length(nb_children)+1), Vectorize(phi))
levels.lower_bounds <- sapply(1:(length(nb_children)+1), Vectorize(function(i) phi(i-1) + 1))

get_level <- function(i){
  min(which(levels.upper_bounds >= i))
}

get_children <- function(i){
  l <- get_level(i)
  r <- i - levels.lower_bounds[l] + 1
  (phi(l) + nb_children[l] * (r-1) + 1):(phi(l) + nb_children[l] * r)
}

get_parent <- function(i){
  l <- get_level(i)
  nb_siblings <- nb_children[l-1]
  r <- i - levels.lower_bounds[l] + 1 # Position of i in its level
  r_parent <- floor((r - 1)/ nb_siblings) + 1 # Position of its parent in its own level
  levels.lower_bounds[l-1] + r_parent - 1
}

# Generate the subsampling tree
subsamples <- list()
subsamples[[1]] <- 1:max(ns)
for(i in 2:max(levels.upper_bounds)){
  subsamples[[i]] <- sample(subsamples[[get_parent(i)]], ns[get_level(i)], replace = F)
}

# Compute sigma_n and fin_diff

compute_sigma_n_and_fin_diff <- function(indices, delta, Delta){
  replicate_data <- observed_data[indices, ]
  
  TMLE_delta <- TMLE_EY1_speedglm(replicate_data, delta)
  TMLE_delta_plus_Delta <- TMLE_EY1_speedglm(replicate_data, delta + Delta)
  sigma_n <- sqrt(TMLE_delta$var_IC)
  fin_diff <- (TMLE_delta_plus_Delta$Psi_n - TMLE_delta$Psi_n) / Delta
  
  c(fin_diff, sigma_n)
}

# Generate one sample
observed_data.list <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, max(ns))
observed_data <- data.frame(L0 = observed_data.list$L0,
                            A0 = observed_data.list$A0,
                            L1 = observed_data.list$L1)

cl <- makeCluster(getOption("cl.cores", 64), outfile = '')
registerDoParallel(cl)

results <- foreach(job = 1:nrow(estimation_tasks), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- estimation_tasks[job, ]$n
                     candidate_rate <- estimation_tasks[job, ]$candidate_rate
                     cat("Beginning of n = ", n, "\n")
                     delta <- n^(-candidate_rate / eta)
                     Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
                     
                     current_level <- which(ns == n)
                     
                     subsamples.results <- 
                       sapply(levels.lower_bounds[current_level]:levels.upper_bounds[current_level],
                              function(i) compute_sigma_n_and_fin_diff(indices = subsamples[[i]],
                                                                       delta = delta,
                                                                       Delta = Delta))
                     fin_diff <- median(subsamples.results[1, ])
                     sigma_n <- median(subsamples.results[2, ])
                     
                     cat("n = ", n, ", delta = ", delta, ", sigma_n = ", sigma_n, "\n")
                     c(n = n, candidate_rate = candidate_rate, delta = delta, sigma_n = sigma_n, 
                       fin_diff = fin_diff)
                   }

pivot <- (abs(results[, "fin_diff"]) * results[, "delta"] / results[, "sigma_n"])^eta * sqrt(n)
results<- cbind(results, pivot = pivot)
log_sigma_n_log_delta <- ggplot(as.data.frame(results), aes(x = log(delta) / log(10),
                                                            y = log(sigma_n) / log(10),
                                                            coulour = factor(candidate_rate))) +
  geom_point() + geom_line() +
  geom_abline(intercept = -0.7, slope = -gamma) +
  geom_abline(intercept = -0.8, slope = -gamma)

log_sigma_n_log_n <- ggplot(as.data.frame(results), aes(x = log(n) / log(10),
                                                        y = log(sigma_n) / log(10),
                                                        coulour = factor(candidate_rate))) +
  geom_point() + geom_line() +
  geom_abline(intercept = -0.5, slope = gamma * true_rate / eta) +
  geom_abline(intercept = -0.7, slope = gamma * true_rate / eta) +
  geom_abline(intercept = -0.9, slope = gamma * true_rate / eta)

log_fin_diff_log_delta <- ggplot(as.data.frame(results), aes(x = log(delta) / log(10),
                                                             y = log(abs(fin_diff)) / log(10),
                                                             coulour = factor(candidate_rate))) +
  geom_point() + geom_line() +
  geom_abline(intercept = -2, slope = -beta)

log_fin_diff_log_n <- ggplot(as.data.frame(results), aes(x = log(n) / log(10),
                                                         y = log(abs(fin_diff)) / log(10),
                                                         coulour = factor(candidate_rate))) +
  geom_point() + geom_line() + 
  geom_abline(intercept = -2, slope = beta * true_rate / eta) +
  geom_abline(intercept = -2.3, slope = beta * true_rate / eta) + 
  geom_abline(intercept = -1.7, slope = beta * true_rate / eta)

log_pivot_log_n <- ggplot(as.data.frame(results), aes(x = log(n) / log(10),
                                                      y = log(pivot) / log(10),
                                                      coulour = factor(candidate_rate))) +
  geom_point() + geom_line()

log_half_pivot_log_n <- ggplot(as.data.frame(results), aes(x = log(n) / log(10),
                                                           y = log((abs(fin_diff) * delta / sigma_n)^eta * sqrt(n)) / log(10),
                                                           coulour = factor(candidate_rate))) +
  geom_point() + geom_line() + 
  geom_abline(intercept = -2, slope = 0) + 
  geom_abline(intercept = -2.5, slope = 0) +
  geom_abline(intercept = -3, slope = 0)

grid.arrange(nrow = 2, ncol = 2, log_pivot_log_n, log_sigma_n_log_n,
             log_half_pivot_log_n, log_fin_diff_log_n)

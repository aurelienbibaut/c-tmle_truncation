beta <- 0.1
gamma <- 0.7 / 2

R <- function(delta, n) delta^(1 - beta) + n^(-0.5) * delta^(-gamma)
R_prime <- function(delta, n) (1 - beta) * delta^(-beta) - gamma * n^(-0.5) * delta^(-gamma - 1)

sigma_R_prime <- function(delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  Delta * (delta^(-beta - 1) + n^(-0.5) * delta^(-gamma - 2)) +
    Delta^(-1) * n^(-0.5) * delta^(-gamma)
  # n^(-0.25) * delta^(-(beta  + gamma + 1)/2)
}

n <- 100000
deltas <- seq(from = 1e-5, to = 1e-1, length = 100)
plot(deltas, sapply(deltas, Vectorize(function(delta) R_prime(delta, n))), ylim = c(-4, 2))
lines(deltas, sapply(deltas, Vectorize(function(delta) (1 - beta) * delta^(-beta))), col = "green")
lines(deltas, sapply(deltas, Vectorize(function(delta) -gamma * n^(-0.5) * delta^(-gamma - 1))), col = "red")
abline(0, 0)

finite_diff_estimator <- function(delta, n) R_prime(delta, n) + sigma_R_prime(delta, n)
# 
delta_n <- vector(); bias_prime_delta_n_star <- vector(); sigma_prime_delta_n_star <- vector()
sigma_R_prime_delta_n_plus <- vector()
epsilon <- 0.2
eta <- 2
# 

LHS_eta <- vector(); LHS_eta_a <- vector(); LHS_eta_b <- vector()

rate_delta_n <- -1 / (2 * (gamma + 1 - beta))
cat("delta_n's rate: ", rate_delta_n, "\n")

ns <- 10^(seq(from = 2, to = 10, length = 50))
delta_n_star <- vector()

for(n in ns){
  delta_n <- c(delta_n,
               uniroot(Vectorize(function(delta) R_prime(delta, n)), interval = c(1/n, 1e-1),
                       extendInt = "yes")$root)
  
  delta_n_plus <- delta_n[length(delta_n)] * n^epsilon
  bias_prime_delta_n_star <- c(bias_prime_delta_n_star,
                               delta_n_plus^(-beta))
  sigma_prime_delta_n_star <- c(sigma_prime_delta_n_star,
                                n^(-0.5) * delta_n_plus^(-gamma - 1))
  
  sigma_R_prime_delta_n_plus <- c(sigma_R_prime_delta_n_plus,
                                  sigma_R_prime(delta_n_plus, n))
  
  delta_n_star_tilde <- n^(-(1 / (2 * eta * (gamma + 1 - beta))))
  
  delta_n_b <- delta_n_star_tilde * n^(+epsilon)
  delta_n_a <- delta_n_star_tilde * n^(-epsilon)
  
  
  LHS_eta<- c(LHS_eta,
               (finite_diff_estimator(delta_n_star_tilde, n) * delta_n_star_tilde^(gamma + 1))^eta / n^(-0.5))
  LHS_eta_a <- c(LHS_eta_a,
                 (finite_diff_estimator(delta_n_a, n) * delta_n_a^(gamma + 1))^eta / n^(-0.5))
  LHS_eta_b <- c(LHS_eta_b,
                 (finite_diff_estimator(delta_n_b, n) * delta_n_b^(gamma + 1))^eta / n^(-0.5))
  # LHS_eta <- c(LHS_eta,
  #              (R_prime(delta_n_star, n) * delta_n_star^(gamma + 1))^eta)
  
  delta_n_star <- c(delta_n_star,
                    try(uniroot(Vectorize(function(delta) (finite_diff_estimator(delta, n) / delta^(-gamma - 1))^eta - n^(-0.5)),
                            interval = c(delta_n_a, delta_n_b),
                            extendInt = "yes")$root))
}

par(mfrow = c(2, 2), mar=c(4.1,4.5,4.9,2.1))
plot(log(ns), log(delta_n))
abline(log(delta_n)[1] + 1/(2*(gamma + 1 - beta)) * log(10^2), -1/(2*(gamma + 1 - beta)))

plot(log(ns), log(bias_prime_delta_n_star), main = "log b'(delta_n^+)")
plot(log(ns), log(sigma_prime_delta_n_star), main = "log n^-0.5 * sigma0(delta_n^+)")

plot(log(ns), log(sigma_R_prime_delta_n_plus), main = "sigma_{R'(delta_n^+)}")

plot(log(ns), log(LHS_eta), type = "l", col = "green", ylim = c(min(c(log(LHS_eta), 
                                                                      log(LHS_eta_a),
                                                                      log(LHS_eta_b))),
                                                                max(c(log(LHS_eta), 
                                                                      log(LHS_eta_a),
                                                                      log(LHS_eta_b)))),
     xlab = expression(log(n)),
     ylab = '',
     main = expression(g(tilde(delta)[n])* hat(" = ") * eta * log* bgroup("(", 
                                                                         frac(widehat(Delta * R[n](tilde(delta)[n])), sigma[0](tilde(delta)[n])*{tilde(delta)[n]}^{-1})
                                                                          ,
                                                                         ")") + 1/2 * log(n)))
abline(0, 0)
lines(log(ns), log(LHS_eta_a), col = "blue")
lines(log(ns), log(LHS_eta_b), col = "red")
legend('bottomleft', c(expression(g(delta[n]^"*" * n^{-epsilon})),
                       expression(g(delta[n]^"*")), 
                       expression(g(delta[n]^"*" * n^epsilon)),
                       expression(0)), 
       col = c("blue", "green", "red", "black"), lty = c(1,1,1, 1))


plot(log(ns), log(LHS_eta) - 0.5 * log(ns), type = "l", col = "green", ylim = c(min(c(log(LHS_eta) - 0.5 * log(ns), 
                                                                                    log(LHS_eta_a) - 0.5 * log(ns),
                                                                                    log(LHS_eta_b) - 0.5 * log(ns))),
                                                                                max(c(log(LHS_eta) - 0.5 * log(ns), 
                                                                                      log(LHS_eta_a) - 0.5 * log(ns),
                                                                                      log(LHS_eta_b)) - 0.5 * log(ns))),
     xlab = expression(log(n)),
     ylab = '',
     main = expression(LHS(tilde(delta)[n]) * hat(" = ") * eta * log* bgroup("(", 
                                                                         frac(widehat(Delta * R[n](tilde(delta)[n])), sigma[0](tilde(delta)[n])*{tilde(delta)[n]}^{-1}),
                                                                         ")" )))

lines(log(ns), log(LHS_eta_a) - 0.5 * log(ns), col = "blue")
lines(log(ns), log(LHS_eta_b) - 0.5 * log(ns), col = "red")
abline(0, -0.5)
legend('bottomleft', c(expression(LHS(delta[n]^"*" * n^{-epsilon})),
                       expression(LHS(delta[n]^"*")), 
                       expression(LHS(delta[n]^"*" * n^epsilon)),
                       expression(-0.5 * log(n))), 
       col = c("blue", "green", "red", "black"), lty = c(1,1,1, 1))



plot(log(ns), log(delta_n_star),
     xlab = expression(log(n)),
     ylab = expression(log(delta[n]^"*")))
abline(0, -1 / (2 * eta * (gamma + 1 - beta)))
legend('bottomleft', c(expression(log(delta[n]^"*")),
  expression(frac(-log(n), 2 * eta * (gamma + 1 - beta)))), lty = c(NA,1), pch = (c(1, NA)))

mtext(expression(eta == 2), side = 3, line = -25, outer = TRUE, cex = 2)

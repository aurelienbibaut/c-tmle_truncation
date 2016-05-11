# X <- outer(seq(from = -9e-4, to = 9e-4, length = 7), c(0,1,2,3), "^")
# y <- X %*% c(-8, 1,-1, 1)
# # y <- y + mean(y) * 1e-1 * rnorm(length(y))
# plot(y)
# lm(y ~ X - 1)
# 
# weights <- c(1e-2, 1e-1, 5e-1, 1, 5e-1, 1e-1, 1e-5)
# W <- diag(weights)
# 
# wfit <- lm.wfit(X, y, weights)$coefficients
# (lm.wfit(X, y, weights))$coefficients
# ginv(t(X)%*% W %*% X) %*% t(X) %*% W %*% y
# A <- apply(diag(nrow(X)), 2, function(C) lm.wfit(X, C, weights)$coefficients)
# wfit$qr$qr %*% t(X) %*% W %*% y
# solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

# Functions to be used later
logit <- function(x){
  log(x/(1-x))
}

expit<-function(x){
  result<-exp(x)/(1+exp(x))
  result[is.nan(result)]<-1
  result
}

n <- 1e3
random_offset <- runif(n, min = -2, max = 2)
W <- runif(n, min = -2, max = 2)
Qbar <- expit(random_offset + 1e2 * W)
Y <- rbinom(n, 1, Qbar)

print(glm(Y ~ W - 1, family = binomial, offset = random_offset)$coefficients[1], subset = 2*(1:500))

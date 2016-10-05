# Treatment rule(s)
alwaysTreated0 <- function(L0){
  1
}

# Functions to be used later
logit <- function(x){
  log(x/(1-x))
}
  
expit <- function(x){
  result <- exp(x)/(1+exp(x))
  result[is.nan(result)] <- 1
  result
}

g_to_g_delta<-function(delta, g){
  (g<delta) * delta + (g>=delta) * g
}
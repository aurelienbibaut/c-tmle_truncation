xs <- seq(from = 0, to = 10, by = 0.01)
f <- function(x){
  if(x <= 3){
    return(0.1 * x)
  }else if(3 < x & x <= 5){
    return(f(3) - 2 * (x - 3))
  }else{
    return (f(5) + 0.7 * (x - 5))
  }
}
fs <- sapply(xs, Vectorize(f))
ys <- fs + 0.2 * rnorm(length(fs))
plot(xs, ys)

library('segmented')
lm.out <- lm(ys ~ xs)
segmented.out <- segmented(lm.out, seg.Z=~xs, psi = list(xs = NA))
plot(segmented.out)

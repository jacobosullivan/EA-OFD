## Tokeshi method
shape <- 1 # set to 1 for unimodal left, 2 for unimodal right, 3 for bimodal

S <- 100
N <- 100
l <- 10
h <- 0.1

if (shape==1) {
  x <- sample(1:N, S, replace=T, prob=exp(-(1:N)/l))  
} else if (shape==2) {
  x <- sample(1:N, S, replace=T, prob=exp((1:N)/l))  
} else if (shape==3) {
  x <- c(sample(1:N, S/2, replace=T, prob=exp(-(1:N)/l)),
         sample(1:N, S/2, replace=T, prob=exp((1:N)/l)))
}

x <- x/max(x)
x <- as.numeric((cut(x, breaks=seq(0,1,by=h))))
hist(x)

x <- table(x)

P <- c()
for (i in 1:length(x)) {
  p <- 0
  for (j in x[i]:S) {
    bc <- factorial(S) / (factorial(j) * factorial(S-j))
    p <- p + (bc * (h^j) * ((1-h)^(S-j)))
  }
  P[i] <- p
}
P
P[1] < 0.05
P[length(P)] < 0.05


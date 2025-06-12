bnpmediation <- function(obj1, obj0, q, NN = 100, n1, n0, extra.thin = 5) {
  library(mnormt)
  Len.MCMC <- 1:dim(obj0)[1]
  if (extra.thin != 0) {Len.MCMC <- Len.MCMC[seq(1, length(Len.MCMC), extra.thin)]}
  Y1 <- obj1[, 1]
  Y0 <- obj0[, 1]
  Y11 <- Y1[Len.MCMC]
  Y00 <- Y0[Len.MCMC]
  mat.given.ij <- function(x, y) ifelse(x <= y, (q - 1) * (x - 1) + y - x * (x - 1) / 2, (q - 1) * (y - 1) + x - y * (y - 1) / 2)
  mat <- function(q) outer(1:q, 1:q, mat.given.ij)
  pb <- txtProgressBar(min = 0, max = length(Len.MCMC), style = 3)
  Y10 <- NULL
  index <- 0
  for (j in Len.MCMC) {
    index <- index + 1; mu2 <- obj0[j, 2:q]
    sigma22 <- matrix(obj0[j, seq(2 * q + 1, q + (q * (q + 1) / 2))][mat(q - 1)], q - 1, q - 1, byrow = TRUE)
    joint0 <- do.call("rbind", replicate(NN, data.frame(sapply(1:n0, function(x) {
      sigma22_x <- sigma22
      rmnorm(1, mu2, sigma22_x)
    }))))
    b01 <- NULL
    Weight.num0 <- matrix(nrow = 1, ncol = n0 * NN)  
    B0 <- matrix(nrow = 1, ncol = n0 * NN) 
    mu1 <- obj1[j, 1] 
    mu2 <- obj1[j, 2:q]
    sigma1 <- obj1[j, q + 1]
    sigma12 <- obj1[j, (q + 2):(2 * q)] 
    sigma22 <- matrix(obj1[j, (2 * q + 1):(q + (q * (q + 1) / 2))][mat(q - 1)], q - 1, q - 1, byrow = TRUE)
    Weight.num0[1, 1:(n0 * NN)] <- dmnorm(joint0, mu2, sigma22)
    b01[1] <- mu1 - sigma12 %*% solve(sigma22) %*% t(t(mu2))
    B0[1, 1:(n0 * NN)] <- sigma12 %*% solve(sigma22) %*% t(joint0)
    Weight <- apply(Weight.num0, 2, function(x) x / sum(x))
    test <- Weight * (b01 + B0)
    Y10[index] <- mean(apply(test, 2, sum))
    setTxtProgressBar(pb, index)
  }
  z <- list(Y11=Y11, Y00=Y00, Y10=Y10, 
            ENIE=mean(Y11-Y10), ENDE=mean(Y10-Y00), 
            ETE=mean(Y11-Y00), 
            TE.c.i=c(sort(Y11-Y00)[length(Len.MCMC)*0.025],
                     sort(Y11-Y00)[length(Len.MCMC)*0.975]),
            IE.c.i=c(sort(Y11-Y10)[length(Len.MCMC)*0.025],
                     sort(Y11-Y10)[length(Len.MCMC)*0.975]),
            DE.c.i=c(sort(Y10-Y00)[length(Len.MCMC)*0.025],
                     sort(Y10-Y00)[length(Len.MCMC)*0.975]))  
  z$call <- match.call()
  class(z) <- "bnpmediation"
  return(z)
}
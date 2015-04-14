
LS <-function(x,y)
{
   solve(t(x)%*%x)%*%t(x)%*%y; 
}


ConvexOpt <-function(x,y, rho=0.00000001, threshold=0.001, iterations=1500)
{
  # regularized matrix version of least squares
  # x = (A'*A)*A'*b
  
  # rho, the wiggle factor
  
  # A = x^t*x
  a  <- t(x)%*%x
  
  # b = x^t*y  
  b <- t(x)%*%y
  
  #identity matrix
  I <- diag(1,nrow(a))
  
  Aprime <- solve(a+(rho*I))
  
  u <- rep(0, nrow(b))
  z <- rep(0, nrow(b))
  x <- rep(0, nrow(b))
  
  RMS <- c()
  for (i in 1:iterations) {
    x <- Aprime%*%(b + rho*(z - u))
    
    #   for each element in x + u
    #   z = (x + u);
    for (j in seq_along(x)) {
      temp <- x[j] + u[j]
      if (temp > 0) {
        z[j] <- temp
      } else {
        z[j] <- 0
      }    
    }
    #   for each element in x - z
    #   u = u + (x - z);
    for (k in seq_along(u)) {
      u[k] <- u[k] + (x[k] - z[k])
    }
    
    sumSquares <- c() # you may have to move this above the first for loop
    #   calculate RMS error
    for (ii in seq_along(x)) {
      sumSquares <- c(sumSquares, (x[ii] - z[ii])*(x[ii] - z[ii]))
    }
    sumSquares <- sum(sumSquares)
    
    RMS <- sqrt(sumSquares/length(x))
    if (RMS < threshold) break  
  }
  
  return(t(t(z)))
}


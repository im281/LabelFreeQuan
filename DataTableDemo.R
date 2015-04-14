
#One thing we get is some very nice and powerful syntax. 
#Consider some simple data of replicate time series:
time <- rep(1:10, 10)
replicate <- sort(time)
value <- rnorm(100)
df <- data.frame(replicate, time, value)


#To apply a function to each set of replicates, instead of
sapply(1:max(df$replicate), function(i) 
  mean( df[df$replicate == i,]$value) 
)

#We can use
dt <- data.table(df)
dt[, mean(value), by="replicate"]

# Note that we could have passed multiple arguments to the function, 
# f(time, value), or functions of the arguments, mean(value*time), 
# etc. While this can be much faster data-frames to begin with (see below), 
# when the function is much more computationally intensive than “mean”, it 
# might be sensible to parallelize this command instead:

grpsize = ceiling(1e7/26^2) 
DF <- data.frame(x=rep(factor(LETTERS), each = 26 * grpsize), 
                 y=rep(factor(letters), each = grpsize), 
                 v=runif(grpsize * 26 ^ 2))
DT <- data.table(DF)
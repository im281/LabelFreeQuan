#Missing Values
x <- c(1,2,NA,10,3)
is.na(x)
x <- c(1,2,NaN,NA,4)
is.na(x)
is.nan(x)

#data types 

#data ftame
x <- data.frame(foo = 1:4, bar = c(T,T,F,F))
x
nrow(x)
ncol(x)

#Names Attribute
x <- 1:3
names(x)
names(x) <- c("foo","bar","norf")
x
names(x)
x <- list(a=1,b=2,c=3)
x

m <- matrix(1:4,nrow =2,ncol =2)
dimnames(m) <-list(c("a","b"),c("d","e"))
m

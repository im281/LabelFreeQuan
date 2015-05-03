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

#assign new column name
colnames(X)[2] <- "superduper"

#drop column named Replicate
mq$Replicate <- null

m <- matrix(1:4,nrow =2,ncol =2)
dimnames(m) <-list(c("a","b"),c("d","e"))
m

#mean of a data table column with missing values removed
mean(dt$nitrate,na.rm = TRUE)

#put 3 digits example 001 
test <-  sprintf("%03d",1)

airquality <- test
airquality[is.na(airquality$nitrate),]
airquality <- airquality[!is.na(airquality$Ozone),]


# 4 figures arranged in 2 rows and 2 columns
attach(mtcars)
par(mfrow=c(2,2))
plot(wt,mpg, main="Scatterplot of wt vs. mpg")
plot(wt,disp, main="Scatterplot of wt vs disp")
hist(wt, main="Histogram of wt")
boxplot(wt, main="Boxplot of wt")

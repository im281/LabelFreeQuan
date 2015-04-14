source('Functions.R')

# Read in the data
#pep <- jnb_PrepareDataForAnalysis('File-11-12.csv')
pep <- jnb_PrepareDataForAnalysis('File5-6-TN.csv')
# Perform the test
pVals <- jnb_PerformTheTest(pep)



# Examine the results
pVals

#resultFileNames <- c('130124_dilA_1_01.pdResult', '130124_dilA_2_01.pdResult')

#1.2 #########################################################################################
resultFileNames1.2 <- c('130124_dilA_1-2.pdResult')
dt1.2 <- GetAllPsmsTablesFromDatabase(resultFileNames1.2)
unique(dt1.2$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt1.2$Condition <- rep('Unknown', nrow(dt1.2))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt1.2[grep('_1_', dt1.2$SpectrumFileName)]$Condition <- 'Dilution01'
dt1.2[grep('_2_', dt1.2$SpectrumFileName)]$Condition <- 'Dilution02'

# tmp <- dt[grep('_1_', dt$SpectrumFileName)]

proteinIDTypes <- c('sol', 'A0A', 'hum')

res1.2 <- MyCoolTtest(dt1.2, proteinIDTypes)

tbl <- PrepDataFromRoc(res1.2)
RocResults <- MyRocR(tbl)
RocResults[['AUC']]
plot(RocResults[[1]], colorize=T)
#head(RocResults[[1]]@y.values)
#head(RocResults[[1]]@x.values)


#8.9#########################################################################################
resultFileNames8.9 <- c('130124_dilA_8-9.pdResult')
dt8.9 <- GetAllPsmsTablesFromDatabase(resultFileNames8.9)
unique(dt8.9$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt8.9$Condition <- rep('Unknown', nrow(dt8.9))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt8.9[grep('_8_', dt8.9$SpectrumFileName)]$Condition <- 'Dilution01'
dt8.9[grep('_9_', dt8.9$SpectrumFileName)]$Condition <- 'Dilution02'

# tmp <- dt[grep('_1_', dt$SpectrumFileName)]

proteinIDTypes <- c('sol', 'A0A', 'hum')

res8.9 <- MyCoolTtest(dt8.9, proteinIDTypes)

preppedRes8.9 <- PrepDataFromRoc(res8.9)
RocResults8.9 <- MyRocR(preppedRes8.9)
RocResults8.9[['AUC']]
plot(RocResults8.9[[1]], colorize=T)


##########################################################################################
#5.6
resultFileNames5.6 <- c('130124_dilA_5-6.pdResult')
dt5.6 <- GetAllPsmsTablesFromDatabase(resultFileNames5.6)
unique(dt5.6$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt5.6$Condition <- rep('Unknown', nrow(dt5.6))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt5.6[grep('_5_', dt5.6$SpectrumFileName)]$Condition <- 'Dilution01'
dt5.6[grep('_6_', dt5.6$SpectrumFileName)]$Condition <- 'Dilution02'

# tmp <- dt[grep('_1_', dt$SpectrumFileName)]

proteinIDTypes <- c('sol', 'A0A', 'hum')

res5.6 <- MyCoolTtest(dt5.6, proteinIDTypes)

preppedRes5.6 <- PrepDataFromRoc(res5.6)
RocResults5.6 <- MyRocR(preppedRes5.6)
RocResults5.6[['AUC']]
plot(RocResults5.6[[1]], colorize=T)
##########################################################################################

#Minora5.6
resultFileNames5.6 <- c('130124_dilA_5_6-Minora.pdResult')
dt5.6 <- GetAllPsmsTablesFromDatabase(resultFileNames8.9)
unique(dt5.6$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt5.6$Condition <- rep('Unknown', nrow(dt5.6))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt5.6[grep('_5_', dt5.6$SpectrumFileName)]$Condition <- 'Dilution01'
dt5.6[grep('_6_', dt5.6$SpectrumFileName)]$Condition <- 'Dilution02'

# tmp <- dt[grep('_1_', dt$SpectrumFileName)]

proteinIDTypes <- c('sol', 'A0A', 'hum')

res5.6 <- MinoraTtest(dt5.6, proteinIDTypes)

preppedRes5.6 <- PrepDataFromRoc(res5.6)
RocResults5.6 <- MyRocR(preppedRes5.6)
RocResults5.6[['AUC']]
plot(RocResults5.6[[1]], colorize=T)



#Minora-8.9#########################################################################################
resultFileNames8.9 <- c('130124_dilA_8_9-Minora.pdResult')
dt8.9 <- GetAllPsmsTablesFromDatabase(resultFileNames8.9)
unique(dt8.9$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt8.9$Condition <- rep('Unknown', nrow(dt8.9))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt8.9[grep('_8_', dt8.9$SpectrumFileName)]$Condition <- 'Dilution01'
dt8.9[grep('_9_', dt8.9$SpectrumFileName)]$Condition <- 'Dilution02'

# tmp <- dt[grep('_1_', dt$SpectrumFileName)]

proteinIDTypes <- c('sol', 'A0A', 'hum')

res8.9 <- MinoraTtest(dt8.9, proteinIDTypes)

preppedRes8.9 <- PrepDataFromRoc(res8.9)
RocResults8.9 <- MyRocR(preppedRes8.9)
RocResults8.9[['AUC']]
plot(RocResults8.9[[1]], colorize=T)


##########################################################################################

#Minora1.2
resultFileNames1.2 <- c('130124_dilA_1_2-Minora.pdResult')
dt1.2 <- GetAllPsmsTablesFromDatabase(resultFileNames1.2)
unique(dt1.2$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt1.2$Condition <- rep('Unknown', nrow(dt1.2))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt1.2[grep('_1_', dt1.2$SpectrumFileName)]$Condition <- 'Dilution01'
dt1.2[grep('_2_', dt1.2$SpectrumFileName)]$Condition <- 'Dilution02'

# tmp <- dt[grep('_1_', dt$SpectrumFileName)]

proteinIDTypes <- c('sol', 'A0A', 'hum')

res1.2 <- MinoraTtest(dt1.2, proteinIDTypes)

preppedRes1.2 <- PrepDataFromRoc(res1.2)
RocResults1.2 <- MyRocR(preppedRes1.2)
RocResults1.2[['AUC']]
plot(RocResults1.2[[1]], colorize=T)
 
#Log2 ratios
############################################################################################

#1.2 #########################################################################################
testFile <- c('MaxQuant-Fraction14.pdResult')
dt1 <- GetAllPsmsTablesFromDatabase(testFile)
unique(dt1$SpectrumFileName)

# Create a new column, 'Condition', and populate this column with the string unknown
# This column will be used to add conditions based on file names.
dt1$Condition <- rep('Unknown', nrow(dt1))
dt1$Replicate <- rep('Unknown', nrow(dt1))

# Subset the dt table for each file type, and assign the appropriate dilution condition.
dt1[grep('_black_', dt1$SpectrumFileName)]$Replicate <- 'Replicate01'
dt1[grep('_green_', dt1$SpectrumFileName)]$Replicate <- 'replicate02'
dt1[grep('_red_', dt1$SpectrumFileName)]$Replicate <- 'replicate03'

dt1[grep('_Ecoli10_', dt1$SpectrumFileName)]$Condition <- '1:3'
dt1[grep('_Ecoli30_', dt1$SpectrumFileName)]$Condition <- '1:1'

proteinIDTypes <- c('Homo sapiens','Escherichia coli')
#proteinIDTypes <- c('Homo sapiens','Escherichia coli')

minArea <- min(dt1$PeakArea, na.rm=TRUE)
fileNames <- unique(dt1$SpectrumFileName)

a <- dt1[grep('23944', dt1$ParentProteinAccessions)]

a <- dt1[grep('Homo sapiens', dt1$ParentProteinDescriptions)]

b <- a[which(a$XCorr_SequestHT > .67)]
c <- subset(b,XCorr_SequestHT > 0.69)


b <- a[which(dt1$Sequence == 'LVLVGDGGTGK')]
c <- subset(dt1,XCorr_SequestHT > 2.5)
maxArea <- max(c$PeakArea, na.rm=TRUE)

c <- subset(dt1,PercolatorPEP < 0.001 & PeakArea > 0)

source('Functions.R')
myresult <- Log2Ratio(c, proteinIDTypes)
human <- data.frame(myresult[1])
ecoli <- data.frame(myresult[2])
hist(human[,3],breaks = 200, main = "Log2Ratio human 1:1",
     xlab = "Log2Ratio",col.lab="red",col = "green")
hist(ecoli[,3],breaks = 200, main = "Log2Ratio ecoli 1:3",
     xlab = "Log2Ratio",col.lab="red",col = "red")


h <- human[,3]
e <- ecoli[,3]

#par(mfrow=c(1,2))
hist(h,breaks = 200,main = "",
     xlab = "Log2Ratio",col.lab="red",col = "green")
#par(new=TRUE)
hist(e,breaks = 200, main = "",
     xlab = "Log2Ratio",col.lab="red",col = "red")

h$Species <- 'human'
e$Species <- 'ecoli'

combined <- rbind(h,e)
ggplot(combined)#,aes(length, fill = veg)) + geom_density(alpha = 0.2)
hist(combined,breaks = 200)

write.csv(c,"MaxQuantData.csv")

a=rnorm(1000, 3, 1)
b=rnorm(1000, 6, 1)
hist(a, xlim=c(0,10), col="red")
hist(b, add=T, col=rgb(0, 1, 0, 0.5) )

#example density plot
######################################################

carrots <- data.frame(length = rnorm(100000, 6, 2))
cukes <- data.frame(length = rnorm(50000, 7, 2.5))

#Now, combine your two dataframes into one.  First make a new column in each.
carrots$veg <- 'carrot'
cukes$veg <- 'cuke'

#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes)

#now make your lovely plot
ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)
#######################################################

ec <-data.frame(e)
hu <-data.frame(h)

ec$Species <- 'Ecoli'
hu$Species <- 'Human'
names(ec) <- c("value","species")
names(hu) <- c("value","species")
combined <-rbind(ec,hu)

length = nrow(combined)

ggplot(combined, aes(length, fill = Species)) + geom_density(alpha = 0.2)




#Two objects: 100 numbers each, randomly distributed with a standard deviation of 1. One has mean of 10, the other with mean of 13. 
a=rnorm(100,mean=10,sd=1)
b=rnorm(100,mean=13,sd=1)

# print the two histograms side-by-side
par(mfrow=c(1,2))
hist(a,main="")
hist(b,main="")

hist(a,xlim=c(5,18),ylim=c(0,30),breaks=10,col=rgb(1,1,0,0.7),main="",xlab="number")
par(new=TRUE)
hist(b,xlim=c(5,18),ylim=c(0,30),breaks=10,col=rgb(0,1,1,0.4),main="",xlab="",ylab="")



par(mfrow=c(1,2))
hist(h,breaks = 200,xlim=c(-10,10),ylim=c(0,60),main = "",xlab = "Log2Ratio",col.lab="red",col="green")
par(new=TRUE)
hist(e,breaks = 200, xlim=c(-10,10),ylim=c(0,60),main = "",xlab = "Log2Ratio",col.lab="red",col = "red")

#Isotope Annotations
################130124_dilA_1_01-(76).pdResult
source('Functions.R')
r <- c('PotatoAnnotations-TN.pdResult')
d <- AssignDilutionReplicates(r)
WriteSqlite(r,"DilutionTable",d)
chro <- GetAnnotationsFromSqlite(r)
############################################


tp <- data.table(read.csv("PotatoAnnotations-TP.csv"))
tn <- data.table(read.csv("PotatoAnnotations-TN.csv"))


chro$Dilution <- rep('Unknown', nrow(chro))
chro[grep('81', chro$FileId)]$Dilution <- 'Dilution-01'
chro[grep('130', chro$FileId)]$Dilution <- 'Dilution-05'
chro[grep('150', chro$FileId)]$Dilution <- 'Dilution-10'


#AAGASAQVLGQEGK
#IANNLPSDQDAIK
r <- tp[TheoreticalSequence == 'IANNLPSDQDAIK']
r <-SnTtest(r)
box(r)

bp <- data.frame(unlist(r[1]),unlist(r[2]),unlist(r[3]),unlist(r[4]),unlist(r[5]))
boxplot(bp)

unique(chro$FileId) 

colnames(chro)

#get maximum peak area per condition
chroT <- chro[,list(MaxSN=max(SignaltoNoise)),
         by=list(FileId,TheoreticalSequence,RetentionTimeAtApex)]






ma <- data.table(read.csv('MA-PotatoFeatures.csv'))
ma <- ma[-1,]
ma$Time <- as.numeric(as.character(ma$Time))
ma$m.z <- as.numeric(as.character(ma$m.z))
ma$X130124_dilA_1_01.raw <- as.numeric(as.character(ma$X130124_dilA_1_01.raw))

df <- read.csv('MA-PotatoFeatures.csv',stringsAsFactors = FALSE)
ma1 <- data.table(read.csv('MA-PotatoFeatures.csv'))

tt <- data.table(read.csv("130124_dilA_1_01_tt.csv"))
x <- Xtopia(tpd1,tt)

source('Functions.R')
source('Xtopia.R')
t <- MassAnalyzer(ma,tp1)
table <- MyBS(ma,503.755)


mat <- MA(tpd1,ma)

#select columns and rows
View(data.frame(ma)[,c(1,7,8,10)])


subset(ma, select=colnames(ma)[2:10])[1:10]

#So, use the subset function to select the columns, and then take the first 10 rows.

#There’s also other ways, but this is the way that I’ve been doing it. 
#I think I read a blog somewhere that said this is the fastest, been optimized for the data.table.


test <- GetUniquePeptides(x)

ma <-data.frame(ma)

select <- c("Time", "m/z", "Mass",  "130124_dilA_1_01.raw")
f1r1 <-ma[,select]
View(f1r1)

colnames(ma)


str(ma)

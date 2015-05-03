data <- fread("MA-Potato-f1-allreps-TP.csv")
udata1 <- data[Dilution != 'NA']
data1 <- fread("MA-Potato-f2-allreps-TP.csv")
udata2 <- data[Dilution != 'NA']

mape1 <- GetUniquePeptides(udata)
test <- rbind(mapeps,mapeps2)

Run <-function(){
  
  for(i in 1:length(peps)){
    
    tbl <-data.table()
    
    pep <- udata[PeptideSequence == peps[i]]
    
    max <- pep[which(pep$Area == max(pep$Area))]
    
    rbind(tbl,max)
  }
  
  return (tbl)
}




data <- data.table(read.csv("TT-Potato-f1-allreps-TP.csv"))
udata1 <- data[Dilution != 'NA']
udata1$Classlabel<- rep('TruePositive', nrow(udata1))

data <- data.table(read.csv("TT-Potato-f2-allreps-TP.csv"))
udata2 <- data[Dilution != 'NA']
udata2$Classlabel<- rep('TruePositive', nrow(udata2))


finaltp <- rbind(udata1,udata2)
source('RocAnalysis.R')
res <- QuanAnova(finaltp,1,2)
final <- tttp10#[ClassLabel != 'unknown']
RocResults <- MyRocR(final)
RocResults[['AUC']]
plot(RocResults[[1]]@x.values[[1]], RocResults[[1]]@y.values[[1]], 
     type='l', xlab=RocResults[[1]]@x.name, ylab=RocResults[[1]]@y.name)

View(tttp10)


test <- finaltp[PeptideSequence == 'APVGLEVQVTPNK']
View(test)
a1 <- test[Dilution == 'Dilution9']
a2 <- test[Dilution == 'Dilution10']

boxplot(a1$Area,a2$Area)
View(a1$Area)
class(a1)

a1 <- c(6.52E+06,  34443.4,	1766.37,	83993.1)
a2 <- c(3.84E+06,  28251.1,	19967.2,	33578.8)
boxplot(a1,a2)


source('MassAnalyzer.R')
#read all annotions
tp <- data.table(read.csv('PotatoAnnotations-TP.csv'))
tp$Classlabel<- rep('TruePositive', nrow(tp))
mafeatures <- fread('MA-potato.csv')
matp4test<- MATest(tp,mafeatures,4)


matp4Results <- MATest(tp,mafeatures,4)
write.csv(matp4Results,"MA-Potato-f4-allreps-TP.csv") 

matp5Results <- MATest(tp,mafeatures,5)
write.csv(matp5Results,"MA-Potato-f5-allreps-TP.csv")

matp9Results <- MATest(tp,mafeatures,9)
write.csv(matp9Results,"MA-Potato-f9-allreps-TP.csv")

matp10Results <- MATest(tp,mafeatures,10)
write.csv(matp10Results,"MA-Potato-f10-allreps-TP.csv")





#read all annotions
tp <- data.table(read.csv('PotatoAnnotations-TP.csv'))
tp$Classlabel<- rep('TruePositive', nrow(tp))

#Open all detected annotations
############################################################
tp1 <- data.table(read.csv("TT-Potato-f4-allreps-TP.csv"))
tp1$Classlabel<- rep('TruePositive', nrow(tp1))

tp2 <- data.table(read.csv("TT-Potato-f5-allreps-TP.csv"))
tp2$Classlabel<- rep('TruePositive', nrow(tp2))

tn1 <- data.table(read.csv("TT-Potato-f4-allreps-TN.csv"))
tn1$Classlabel<- rep('TrueNegative', nrow(tn1))

tn2 <- data.table(read.csv("TT-Potato-f5-allreps-TN.csv"))
tn2$Classlabel<- rep('TrueNegative', nrow(tn2))
############################################################

finaltp <- rbind(tp1,tp2)
test <- QuanAnova(finaltp,4,5)
final <- test[ClassLabel != 'unknown']
final <- test[Estimate != 0]
View(final)
logratio <- log(final$Estimate,2)
hist(logratio,breaks = 250)
View(logratio)

#Perform ROC analysis Treetop
############################################################

finaltp <- rbind(tp1,tp2)
finaltn <- rbind(tn1,tn2)
finalAnalysis <- rbind(finaltp,finaltn)
res <- QuanAnova(finalAnalysis,4,5)
final <- res[ClassLabel != 'unknown']
RocResults <- MyRocR(final)
RocResults[['AUC']]
plot(RocResults[[1]]@x.values[[1]], RocResults[[1]]@y.values[[1]], 
     type='l', xlab=RocResults[[1]]@x.name, ylab=RocResults[[1]]@y.name)



source('Functions.R')
#Open all detected annotations
############################################################
data <- fread("MA-Potato-f4-allreps-TP.csv")
udata <- data[Dilution != 'NA']
tp1 <- udata#GetUniquePeptides(udata)
tp1$Classlabel<- rep('TruePositive', nrow(tp1))

data <- fread("MA-Potato-f5-allreps-TP.csv")
udata <- data[Dilution != 'NA']
tp2 <- udata#GetUniquePeptides(udata)
tp2$Classlabel<- rep('TruePositive', nrow(tp2))

data <- fread("MA-Potato-f4-allreps-TN.csv")
udata <- data[Dilution != 'NA']
tn1 <- udata#GetUniquePeptides(udata)
tn1$Classlabel<- rep('TrueNegative', nrow(tn1))

data <- fread("MA-Potato-f5-allreps-TN.csv")
udata <- data[Dilution != 'NA']
tn2 <- udata#GetUniquePeptides(udata)
tn2$Classlabel<- rep('TrueNegative', nrow(tn2))
############################################################
finaltp <- rbind(tp1,tp2)
finaltn <- rbind(tn1,tn2)
finalAnalysis <- rbind(finaltp,finaltn)
res <- QuanAnova(finalAnalysis,4,5)
final <- res[ClassLabel != 'unknown']
RocResults <- MyRocR(final)
RocResults[['AUC']]
plot(RocResults[[1]]@x.values[[1]], RocResults[[1]]@y.values[[1]], 
     type='l', xlab=RocResults[[1]]@x.name, ylab=RocResults[[1]]@y.name)



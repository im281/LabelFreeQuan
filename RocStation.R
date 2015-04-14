source('Xtopia.R')
source('MassAnalyer.R')
source('Functions.R')
source('RocAnalysis.R')


#Read all annotations
#############################################################
#read all annotions
tp <- data.table(read.csv('PotatoAnnotations-TP.csv'))
tp$Classlabel<- rep('TruePositive', nrow(tp))

tn <- data.table(read.csv('PotatoAnnotations-TN.csv'))
tn$Classlabel<- rep('TrueNegative', nrow(tn))

#run Treetop tests (Feature lists from unknown feature detection against annotation list)
#############################################################
files <- c("130124_dilA_8_01_tt.csv","130124_dilA_8_02_tt.csv",
           "130124_dilA_8_03_tt.csv","130124_dilA_8_04_tt.csv")

files <- c("130124_dilA_6_01_xt.xls.csv","130124_dilA_6_02_xt.xls.csv",
           "130124_dilA_6_03_xt.xls.csv","130124_dilA_6_04_xt.xls.csv")

tp6results <- Test(tp,files,6)

write.csv(tp6results,"TT-Potato-f6-allreps-TP1.csv")

#MA tests

############################################################
mafeatures <- data.table(read.csv('MA-potato.csv'))

mafeatures <- fread('MA-potato.csv')
# Detect True Positives
#############################################################
matp1Results <- MATest(tp,mafeatures,1)
write.csv(matp1Results,"MA-Potato-f1-allreps-TP.csv")

matp2Results <- MATest(tp,mafeatures,2)
write.csv(matp2Results,"MA-Potato-f2-allreps-TP.csv")

matp4Results <- MATest(tp,mafeatures,4)
write.csv(matp4Results,"MA-Potato-f4-allreps-TP.csv")

matp5Results <- MATest(tp,mafeatures,5)
write.csv(matp5Results,"MA-Potato-f5-allreps-TP.csv")

matp6Results <- MATest(tp,mafeatures,6)
write.csv(matp6Results,"MA-Potato-f6-allreps-TP.csv")

matp7Results <- MATest(tp,mafeatures,7)
write.csv(matp7Results,"MA-Potato-f7-allreps-TP.csv")

matp8Results <- MATest(tp,mafeatures,8)
write.csv(matp8Results,"MA-Potato-f8-allreps-TP.csv")

matp9Results <- MATest(tp,mafeatures,9)
write.csv(matp9Results,"MA-Potato-f9-allreps-TP.csv")

matp10Results <- MATest(tp,mafeatures,10)
write.csv(matp10Results,"MA-Potato-f10-allreps-TP.csv")

#Detect True Negatives
##############################################################
matn6Results <- MATest(tn,mafeatures,6)
write.csv(matn6Results,"MA-Potato-f6-allreps-TN.csv")

matn7Results <- MATest(tn,mafeatures,7)
write.csv(matn7Results,"MA-Potato-f7-allreps-TN.csv")

matn8Results <- MATest(tn,mafeatures,8)
write.csv(matn8Results,"MA-Potato-f8-allreps-TN.csv")

matn5Results <- MATest(tn,mafeatures,5)
write.csv(matn5Results,"MA-Potato-f5-allreps-TN.csv")

matn9Results <- MATest(tn,mafeatures,9)
write.csv(matn9Results,"MA-Potato-f9-allreps-TN.csv")

matn10Results <- MATest(tn,mafeatures,10)
write.csv(matn10Results,"MA-Potato-f10-allreps-TN.csv")






#Open all detected annotations
############################################################
tp1 <- data.table(read.csv("TT-Potato-f9-allreps-TP1.csv"))
tp1$Classlabel<- rep('TruePositive', nrow(tp1))

tp2 <- data.table(read.csv("TT-Potato-f10-allreps-TP1.csv"))
tp2$Classlabel<- rep('TruePositive', nrow(tp2))

tn1 <- data.table(read.csv("TT-Potato-f9-allreps-TN1.csv"))
tn1$Classlabel<- rep('TrueNegative', nrow(tn1))

tn2 <- data.table(read.csv("TT-Potato-f10-allreps-TN1.csv"))
tn2$Classlabel<- rep('TrueNegative', nrow(tn2))
############################################################


#Perform ROC analysis Treetop
############################################################

finaltp <- rbind(tp1,tp2)
finaltn <- rbind(tn1,tn2)
finalAnalysis <- rbind(finaltp,finaltn)
res <- QuanAnova(finalAnalysis,9,10)
final <- res[ClassLabel != 'unknown']
RocResults <- MyRocR(final)
RocResults[['AUC']]
plot(RocResults[[1]]@x.values[[1]], RocResults[[1]]@y.values[[1]], 
     type='l', xlab=RocResults[[1]]@x.name, ylab=RocResults[[1]]@y.name)
############################################################



#Open all detected annotations
############################################################
data <- fread("MA-Potato-f7-allreps-TP.csv")
tp1 <- data
tp1$Classlabel<- rep('TruePositive', nrow(tp1))

data <- fread("MA-Potato-f8-allreps-TP.csv")
tp2 <- data
tp2$Classlabel<- rep('TruePositive', nrow(tp2))

data <- fread("MA-Potato-f7-allreps-TN.csv")
tn1 <- data
tn1$Classlabel<- rep('TrueNegative', nrow(tn1))

data <- fread("MA-Potato-f8-allreps-TN.csv")
tn2 <- data
tn2$Classlabel<- rep('TrueNegative', nrow(tn2))
#Run Roc Mass Analyzer
############################################################
finaltp <- rbind(tp1,tp2)
finaltn <- rbind(tn1,tn2)
finalAnalysis <- rbind(finaltp,finaltn)
res <- QuanAnova(finalAnalysis,7,8)
write.csv(res,"MA-ttest-7-8.csv")
final <- res[ClassLabel != 'unknown']
RocResults <- MyRocR(final)
RocResults[['AUC']]
plot(RocResults[[1]]@x.values[[1]], RocResults[[1]]@y.values[[1]], 
     type='l', xlab=RocResults[[1]]@x.name, ylab=RocResults[[1]]@y.name)

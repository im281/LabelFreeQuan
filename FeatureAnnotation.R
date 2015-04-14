source('Xtopia.R')
source('MassAnalyer.R')
source('Functions.R')

#read all annotions
tp <- data.table(read.csv('PotatoAnnotations-TP.csv'))
tp$Classlabel<- rep('TruePositive', nrow(tp))

tn <- data.table(read.csv('PotatoAnnotations-TN.csv'))
tn$Classlabel<- rep('TrueNegative', nrow(tn))

tpList <- tp[Dilution == "Dilution5"]
tpList <- tpList[Replicate == "Rep3"]
tnList <- tn[Dilution == "Dilution5"]
tnList <- tnList[Replicate == "Rep3"]


#read the output for the specific csv files each detector outputs
#mafeatures <- data.table(read.csv('MA-potato.csv'))
ttfeatures <- data.table(read.csv("130124_dilA_5_03_tt.csv"))
ttfeatures$Dilution <- rep('unkown', nrow(ttfeatures))
ttfeatures$Replicate <- rep('unknown', nrow(ttfeatures))

#ttfeatures$Dilution <- 'Dilution1'
#ttfeatures$Replicate<- 'Rep1'


#match the outputs per detector (untargeted modes) against the annotation list
TTtpresults <- GetXtopiaMatches(tpList,ttfeatures)

#get the best matches
TTtpresults <- GetUniquePeptides(TTtpresults)

View(TTtpresults)

#match the outputs per detector (untargeted modes) against the annotation list
TTtnresults <- Xtopia(tnList,ttfeatures)
#get the best matches
TTtnresults <- GetUniquePeptides(TTtnresults)

View(TTtnresults)

MAresults <- MA(tpList,mafeatures)
#get the best matches
MAresults <- GetUniquePeptides(MAresults)

View(MAresults)


files <- c("130124_dilA_2_01_tt.csv","130124_dilA_2_02_tt.csv",
           "130124_dilA_2_03_tt.csv","130124_dilA_2_04_tt.csv")

tnresults <- Test(tp,files,2)




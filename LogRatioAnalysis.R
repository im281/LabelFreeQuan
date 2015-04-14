source('Functions.R')
source('MassAnalyzer.R')
table <- GetPsmDataTableFromSqlite('MaxQuant-Fraction14.pdResult')
#Minora-QE-LogRatios.pdResult
#get maximum peak area per condition
dt <- table[,list(MaxArea=max(PeakArea)),
         by=list(Sequence, SpectrumFileName)]

table <- table[PeakArea > 0]
table <- table[PercolatorPEP < 0.00001]


dt <- table[grep('P0A', table$ParentProteinAccessions)]


test <- GLA(dt)

ma <-fread('MA-Fraction_14.csv')
MALogRatioTest(dt,ma)

peakList <- data.frame(ma)[,c(1,7,8,9,10,11,12,13,14,15)]


dt <- GetPsmDataTableFromSqlite('Minora-QE-LogRatios.pdResult')
h <-  dt[grep('OS=Homo sapiens', ParentProteinDescriptions)]
e <-  dt[grep('OS=Escherichia coli', ParentProteinDescriptions)]

t <- LogRatioAnalysis(e,'OS=Escherichia coli')
rt <- ComputeLogRatios(t)


#get maximum peak area per condition
test <- pep[,list(MaxArea=max(Area)),
         by=list(Sequence,SpectrumFileName)]

#dt1.2[grep('_2_', dt1.2$SpectrumFileName)]$Condition <- 'Dilution02

test <- test[grep('10ng',test$SpectrumFileName)]


max <- pep[which(pep$Area == max(pep$Area))]

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

#######################################################################################
dt <- GetPsmDataTableFromSqlite('LogRatio-TT-PPD-10-15ng.pdResult')
h <-  dt[grep('OS=Homo sapiens', ParentProteinDescriptions)]
e <-  dt[grep('OS=Escherichia coli', ParentProteinDescriptions)]

h <- h[PercolatorPEP < 0.001]
h <- h[PeakArea > 0]

e <- e[PercolatorPEP < 0.001]
e <- e[PeakArea > 0]

hum <- LogRatioAnalysisTT(h,'OS=Homo sapiens')
ecol <- LogRatioAnalysisTT(e,'OS=Escherichia coli')
humrt <- ComputeLogRatios(hum,10,15)
ecolrt <- ComputeLogRatios(ecol,10,15)


humrt$Ratio <- log(humrt$Ratio,2)
humrt <- humrt[Ratio != '-Inf']

ecolrt$Ratio <- log(ecolrt$Ratio,2)
ecolrt <- ecolrt[Ratio != '-Inf']

humrt <- na.omit(humrt)
ecolrt <- na.omit(ecolrt)


plot(density(humrt$Ratio),col="red",xlim=c(-3,3),main = "Xtopia-Target-PPD Log 2 Ratio 1.5:1") 
lines(density(ecolrt$Ratio),col = "blue")

sd(humrt$Ratio)
mean(humrt$Ratio)

sd(ecolrt$Ratio)
mean(ecolrt$Ratio)

#########################################################################################




hist(humrt$Ratio2,col="red",xlim=c(-2,2),main = "Log 2 Rato 2:1",breaks = 250) 
lines(ecolrt$Ratio2,col = "blue",breaks = 250)

library('lessR')
# generate 100 random normal data values
y <- rnorm(100)

# normal curve and general density curves superimposed over histogram
# all defaults
color.density(y)

# suppress the histogram, leaving only the density curves
# specify x-axis label per the xlab option for the plot function
color.density(y, col.hist="transparent", xlab="My Var")

# specify (non-transparent) colors for the curves,
# to make transparent, need alpha option for the rgb function
color.density(y, col.fill.nrm="darkgreen", col.fill.gen="plum")

# display only the general estimated density
#  so do not display the estimated normal curve
# specify the bandwidth for the general density curve,
#  use the standard bw option for the density function
color.density(y, type="general", bw=.6)

# display only the general estimated density and a corresponding
#  interval of unit width around x.pt
color.density(y, type="general", x.pt=2)






View(rt)

density(rt$Ratio2)
hist(rt$Ratio2,breaks = 250)

#get maximum peak area per condition
test <- pep[,list(MaxArea=max(Area)),
         by=list(Sequence,SpectrumFileName)]

#dt1.2[grep('_2_', dt1.2$SpectrumFileName)]$Condition <- 'Dilution02

test <- test[grep('10ng',test$SpectrumFileName)]


max <- pep[which(pep$Area == max(pep$Area))]





library("RSQLite")
drv <- dbDriver("SQLite")
con <- dbConnect(drv,"130124_dilA_1_01.pdResult")
data = dbGetQuery(con,"SELECT * from TargetPsms")
head(data)
str(data)

selectS <- select *,max(Xcorr_SequestHT) from TargetPsms where ParentProteinAccessions LIKE "sol%"
group by Sequence

# Read SQlite data base for Chawade et al., 2014 file 1 and 2 (highest concentration)
#Extract all unique peptides with maximum Xcorr score for feature annotation
library("RSQLite")
drv <- dbDriver("SQLite")
con <- dbConnect(drv,"130124_dilA_1_01.pdResult")
data = dbGetQuery(con, "select *,max(Xcorr_SequestHT) from TargetPsms where PercolatorqValue < 0.001 and  ParentProteinAccessions LIKE 'A0%' group by Sequence")
head(data)
str(data)

View(data)

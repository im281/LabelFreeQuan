require(reshape2)
require(data.table)
require(RSQLite)
require(ROCR)
require(pROC)

lappend <- function(lst, obj, name=NULL) {
  lst[[length(lst)+1]] <- obj
  if (!is.null(name)) {
    names(lst)[length(lst)] <- name
  }
  return(lst)
}

#' @title Prepare Data For Analysis
#' @author Joe Brown 
#' @description This function just reads in the 
#' data from the csv file that Iman provided,
#' and prepares the format for the data.table
#' that will be used in the t-test.
#' @param x A character string representing the
#' peptide file path that Iman sent.
jnb_PrepareDataForAnalysis <- function(x) {
  options(warn=-1)
  
  # Read in the data and parse the file
  peps <- read.csv(x, header=T)
  
  # Handle the first set of peptides
  peps1 <- peps[2:nrow(peps),1:5]
  
  #I call unlist because I was getting the row back as a list, 
  #and you can’t use a list to assign column names, you have to 
  #use a vector instead. The unlist, just converts the a list into a vector.
  colnames(peps1) <- unlist(peps[1,1:5])
    
  rownames(peps1) <- peps[2:nrow(peps),1]
  
  
  #The melt function comes from the reshape2 package, and it is a wonderful tool. 
  #It takes a data.frame and melts it down according to 3 things:
  #1.The identifiers that you specify. There’s an option where you don’t 
  #  have to specify id.vars, and I can’t remember what it takes then, rownames or 
  #  first column – something like that.  
  #2.       The column names, puts them in a column called ‘variable’.
  #3.       The data values, puts them in a column called ‘value’.  
  mData1 <- melt(peps1, id.vars='PeptideSequence')
  mData1 <- data.frame(
    FileID=rep('File11', nrow(mData1)),
    mData1)
  
  # Handle the second set of peptides
  peps2Rownames <- unique(peps[2:nrow(peps),7])
  peps2Rownames <- peps2Rownames[peps2Rownames != '']
  peps2 <- data.frame(
    PeptideSequence=peps2Rownames,
    peps[2:(length(peps2Rownames)+1),8:ncol(peps)])
  colnames(peps2)[2:ncol(peps2)] <- as.character(unlist(peps[1,8:ncol(peps)]))
  rownames(peps2) <- peps2Rownames
  mData2 <- melt(peps2, id.vars='PeptideSequence')
  mData2 <- data.frame(
    FileID=rep('File12', nrow(mData2)),
    mData2)
  
  # Combine peps1 and peps2
  pepData <- data.table(rbind(mData1, mData2))
  setkey(pepData, 'PeptideSequence') # Key the Peptide sequences for fast lookup
  
  return(pepData)
}

#' @title Perform the Test
#' @author Joe Brown
#' @description This method iterates through
#' the peptide sequences and selects the data
#' for each peptide individually. The data is 
#' only considered if the peptide was detected
#' in both conditions (FileID 11 and 12). The
#' peptide abundance values are then log transformed.
#' Missing values are imputed with a 0. A t-test
#' is performed on them, and finally the p-values
#' are adjusted using Benjamini-Hochberg.
#' @param x A data.table containing the data
#' to perform the test on. The PeptideSequence
#' column should already be keyed for fast 
#' lookup.
jnb_PerformTheTest <- function(x) {
  pep <- unique(x$PeptideSequence)
  lst <- list()
  
  for (i in 1:length(pep)) {
    # get the peptide
    peptideData <- x[J(pep[i])] 
    
    # Check that the peptide is present in both files    
    if (length(unique(peptideData$FileID))>1) {
      # impute missing values with 0
      peptideData[peptideData$value=='', 'value'] <- 1
      
      # Log transform value
      peptideData[,'value'] <- log(as.numeric(peptideData$value), 10)
      
      # run the t-test
      res <- t.test(value ~ FileID, peptideData)
      lst <- lappend(lst,
        data.frame(
          PeptideSequence=pep[i],
          Estimate=(res$estimate[2]-res$estimate[1]),
          Pvalue=res$p.value))
    }
  }
  res <- do.call(rbind, lst)
  res <- data.frame(
    res,
    AdjPvalue=p.adjust(as.numeric(r $Pvalue)))
  return(data.table(res))
}


#' @description
#' @param x Path to the database
#' @return A data.table representing the Psms table in the database.
GetPsmDataTableFromSqlite <- function(x) {
  
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, x)
  data = data.table(dbGetQuery(con,"SELECT * from TargetPsms"))
  setkey(data, 'SpectrumFileName')
  return(data)
}


#' @description
#' @param x Path to the database
#' @return A data.table representing the CHromatographicComponents table in the database.
GetChroTableFromSqlite <- function(x) {
  
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, x)
  data = data.table(dbGetQuery(con,"SELECT * from ChromatographicComponents"))
  setkey(data, 'FileId')
  return(data)
}

GetAnnotationsFromSqlite <- function(x) {
 
#     finalTable <- data.table(
#     FileID = character(),
#     PeptideSequence = character(),
#     ApexRetentionTime = numeric(),
#     Area = numeric(),
#     SN = numeric(),
#     MaxSN = numeric())
  
  
  # open table and get chromatograhic components
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, x)
  #data = data.table(dbGetQuery(con,"SELECT * from ChromatographicComponents"))
  data = data.table(dbGetQuery(con,"SELECT * from DilutionTable"))
  
  # Get unique list of peptide sequences and associated protein id
  uniquePep <- unique(data$TheoreticalSequence)

  pepList <- list()
  for(i in 1:length(uniquePep)){
           
#     sel <- paste("select FileId,TheoreticalSequence,MonoisotopicMass,TheoreticalTargetMass, ChargeState, PpmError,IsotopeModelScore,RetentionTimeAtApex,PeakArea,  
#                  SignaltoNoise, max(ChromatographicComponents.[PeakArea]),Dilution,Replicate ",
#       "as MaxSN from ChromatographicComponents where ",
#       "ChromatographicComponents.[IsotopeModelScore] > 0 and TheoreticalSequence ",
#       "LIKE \'", uniquePep[i], "\' group by ChromatographicComponents.[FileId]", sep="")
#    
    
    #here we select peptides from each dilution and each replicate where we return the maximum peak area
    # with the additional constraint that they must fit the theoretical isotope distribution
    sel <- paste("select FileId,TheoreticalSequence,MonoisotopicMass,TheoreticalTargetMass, ChargeState, PpmError,IsotopeModelScore,RetentionTimeAtApex,PeakArea,  
                 SignaltoNoise, max(DilutionTable.[PeakArea])",
                "as MaxPeakArea,Dilution,Replicate  from DilutionTable where ",
                "DilutionTable.[IsotopeModelScore] > 0 and TheoreticalSequence ",
                "LIKE \'", uniquePep[i], "\' group by DilutionTable.FileId", sep="")
   
    data <- data.table(dbGetQuery(con, sel))
    
    pepList <- lappend(pepList, data)
    
  }
  finalTable <- do.call(rbind, pepList)
 
  setkey(finalTable, 'FileId')
  
  return(finalTable)
}



GetAllPsmsTablesFromDatabase <- function(x) {
  dt <- data.table()
  for (i in 1:length(x)) {
    dt <- rbind(dt, GetPsmDataTableFromSqlite(x[i]))
  }
  return(dt)
}

GLA <- function(x){
  
  finaltbl <- data.table(
  Sequence = character(),
  massovercharge = numeric(),
  charge = numeric(),
  RT = numeric(),
  Area = numeric())
  
  uniquePeps <- unique(x$Sequence)
  for(i in 1:length(uniquePeps)){   
    temp <- x[Sequence == uniquePeps[i]] 
    
    maxDt <- temp[,MaxArea:=max(PeakArea), by=Sequence]
    
    temp <- maxDt[PeakArea == MaxArea]
    
    finaltbl <- rbind(finaltbl,data.table(
      Sequence = temp$Sequence,
      massovercharge = temp$MassOverCharge,
      charge = temp$Charge,
      RT = temp$ApexRetentionTime,
      Area = temp$PeakArea))
  }
  
  return (finaltbl)
}



#' @title Anova for Label-Free Quan
#' @author Joe Brown and Iman Mohtashemi
#' @description This function just reads in the 
#' data from the csv file that Iman provided,
#' and prepares the format for the data.table
#' that will be used in the t-test.
#' @param x A character string representing the
#' peptide file path that Iman sent.
MyCoolTtest <- function(x, peptideTypes) {
  
  minArea <- min(x$PeakArea, na.rm=TRUE)
  fileNames <- unique(x$SpectrumFileName)
  
  res <- list()
  # Create results table for each peptideType, store in a list
  for (i in 1:length(peptideTypes)) {
    tmpDataTable <- data.table(
      Peptide = character(),
      ComparisonCondition = character(),
      Estimate = numeric(),
      Pvalue = numeric())
    res <- lappend(res, tmpDataTable, peptideTypes[i])
    
    # Subset x for the specifice peptideType
    pepTypeDt <- x[like(ParentProteinAccessions, peptideTypes[i])]    
    setkey(pepTypeDt, 'Sequence')
    
    # Get unique list of peptide sequences and associated protein id
    uniquePep <- unique(pepTypeDt$Sequence)
    
    # Subset x for each peptide
    for (p in 1:length(uniquePep)) {
      dt <- pepTypeDt[J(uniquePep[p])]
      
      # Find any missing file names in subsetted data.table
      missingFileNames <- setdiff(fileNames, unique(dt$SpectrumFileName))  
      
      # Create a temporary data.table to store the missing file names
      missingDt <- pepTypeDt[which(pepTypeDt$SpectrumFileName%in%missingFileNames)]
      missingDt <- unique(data.table(
        Condition=missingDt$Condition,
        SpectrumFileName=missingDt$SpectrumFileName))
      
      # Interpolate / Gap-fill the PeakAreas in missing file name table
      missingDt$MaxArea = rep(minArea, nrow(missingDt))
      
      
      dt[is.na(PeakArea)]$PeakArea <- minArea
      dt <- dt[,list(MaxArea=max(log(PeakArea, 10))),
                   by=list(Condition, SpectrumFileName)]
      
      dt <- rbind(dt, missingDt)
      
      # Check that peptide is found in more than one file,
      if (length(unique(dt$Condition)) > 1) {
        # Perform ANOVA
        lmResult <- lm(MaxArea ~ Condition, dt)
        anovaResult <- anova(lmResult)
        
        conditionLevels <- unlist(lmResult$xlevels)
        
        for (l in 2:length(conditionLevels)) {
          est <- summary(lmResult)$coefficients[,'Estimate'][l]
          p.val <- summary(lmResult)$coefficients[,'Pr(>|t|)'][l]
          
          # Add result row to result data.table in list
          res[[i]] <- rbind(res[[i]],
                            data.table(
                              Peptide=uniquePep[p],
                              ComparisonCondition=conditionLevels[l],
                              Estimate=est,
                              Pvalue=p.val))
        }
      }      
    }
  } 
  
  return(res)
}

MinoraTtest <- function(x, peptideTypes) {
  
  minArea <- min(x$Area, na.rm=TRUE)
  fileNames <- unique(x$SpectrumFileName)
  
  res <- list()
  # Create results table for each peptideType, store in a list
  for (i in 1:length(peptideTypes)) {
    tmpDataTable <- data.table(
      Peptide = character(),
      ComparisonCondition = character(),
      Estimate = numeric(),
      Pvalue = numeric())
    res <- lappend(res, tmpDataTable, peptideTypes[i])
    
    # Subset x for the specifice peptideType
    pepTypeDt <- x[like(ParentProteinAccessions, peptideTypes[i])]    
    setkey(pepTypeDt, 'Sequence')
    
    # Get unique list of peptide sequences and associated protein id
    uniquePep <- unique(pepTypeDt$Sequence)
    
    # Subset x for each peptide
    for (p in 1:length(uniquePep)) {
      dt <- pepTypeDt[J(uniquePep[p])]
      
      # Find any missing file names in subsetted data.table
      missingFileNames <- setdiff(fileNames, unique(dt$SpectrumFileName))  
      
      # Create a temporary data.table to store the missing file names
      missingDt <- pepTypeDt[which(pepTypeDt$SpectrumFileName%in%missingFileNames)]
      missingDt <- unique(data.table(
        Condition=missingDt$Condition,
        SpectrumFileName=missingDt$SpectrumFileName))
      
      # Interpolate / Gap-fill the PeakAreas in missing file name table
      missingDt$MaxArea = rep(minArea, nrow(missingDt))
      
      
      dt[is.na(Area)]$Area <- minArea
      dt <- dt[,list(MaxArea=max(log(Area, 10))),
               by=list(Condition, SpectrumFileName)]
      
      dt <- rbind(dt, missingDt)
      
      # Check that peptide is found in more than one file,
      if (length(unique(dt$Condition)) > 1) {
        # Perform ANOVA
        lmResult <- lm(MaxArea ~ Condition, dt)
        anovaResult <- anova(lmResult)
        
        conditionLevels <- unlist(lmResult$xlevels)
        
        for (l in 2:length(conditionLevels)) {
          est <- summary(lmResult)$coefficients[,'Estimate'][l]
          p.val <- summary(lmResult)$coefficients[,'Pr(>|t|)'][l]
          
          # Add result row to result data.table in list
          res[[i]] <- rbind(res[[i]],
                            data.table(
                              Peptide=uniquePep[p],
                              ComparisonCondition=conditionLevels[l],
                              Estimate=est,
                              Pvalue=p.val))
        }
      }      
    }
  } 
  
  return(res)
}

FunctionsQuanAnova <- function(x, peptideTypes) {
  
  minArea <- min(x$Area, na.rm=TRUE)
  #fileNames <- unique(x$SpectrumFileName)
  
  res <- list()
  # Create results table for each peptideType, store in a list
  for (i in 1:length(peptideTypes)) {
    tmpDataTable <- data.table(
      Peptide = character(),
      ComparisonCondition = character(),
      Estimate = numeric(),
      Pvalue = numeric())
    res <- lappend(res, tmpDataTable, peptideTypes[i])
    
    # Subset x for the specifice peptideType
    pepTypeDt <- x    
    #setkey(pepTypeDt, 'PeptideSequence')
    
    # Get unique list of peptide sequences and associated protein id
    uniquePep <- unique(pepTypeDt$PeptideSequence)
    
    # Subset x for each peptide
    for (p in 1:length(uniquePep)) {
      dt <- pepTypeDt[J(uniquePep[p])]
      
      # Find any missing file names in subsetted data.table
      #missingFileNames <- setdiff(fileNames, unique(dt$SpectrumFileName))  
      
      # Create a temporary data.table to store the missing file names
      #missingDt <- pepTypeDt[which(pepTypeDt$SpectrumFileName%in%missingFileNames)]
      #missingDt <- unique(data.table(
      #  Condition=missingDt$Condition,
      #  SpectrumFileName=missingDt$SpectrumFileName))
      
      # Interpolate / Gap-fill the PeakAreas in missing file name table
      missingDt$MaxArea = rep(minArea, nrow(missingDt))
      
      
      dt[is.na(PeakArea)]$PeakArea <- minArea
      dt <- dt[,list(MaxArea=max(log(PeakArea, 10))),
               by=list(Condition, SpectrumFileName)]
      
      dt <- rbind(dt, missingDt)
      
      # Check that peptide is found in more than one file,
      if (length(unique(dt$Condition)) > 1) {
        # Perform ANOVA
        lmResult <- lm(MaxArea ~ Condition, dt)
        anovaResult <- anova(lmResult)
        
        conditionLevels <- unlist(lmResult$xlevels)
        
        for (l in 2:length(conditionLevels)) {
          est <- summary(lmResult)$coefficients[,'Estimate'][l]
          p.val <- summary(lmResult)$coefficients[,'Pr(>|t|)'][l]
          
          # Add result row to result data.table in list
          res[[i]] <- rbind(res[[i]],
                            data.table(
                              Peptide=uniquePep[p],
                              ComparisonCondition=conditionLevels[l],
                              Estimate=est,
                              Pvalue=p.val))
        }
      }      
    }
  } 
  
  return(res)
}


MyRoc <-function(x){
  
  a <- data.frame(x[1])
  b <- data.frame(x[2])
  
  a$ClassLabel <-rep('unknown',nrow(a))
  b$ClassLabel <-rep('unknown',nrow(b))
  
  
  for(i in 1:nrow(a)){  
    a[i,]$ClassLabel <- 'TruePositive'  
  }
  for(i in 1:nrow(b)){  
    b[i,]$ClassLabel <- 'TrueNegative'  
  }
  
  setnames(b, colnames(b), colnames(a))
  tbl <- do.call(rbind, list(a, b))
    
  # Basic example
  roc(tbl$ClassLabel, tbl$sol.Pvalue,
      levels=c("TruePositive", "TrueNegative"))
  # As levels a$ClassLabel == c("Good", "Poor"),
  # this is equivalent to:
  roc(tbl$ClassLabel, tbl$sol.Pvalue)
  # In some cases, ignoring levels could lead to unexpected results
  # Equivalent syntaxes:
  roc(ClassLabel ~ sol.Pvalue, a)
  roc(tbl$ClassLabel ~ tbl$sol.Pvalue)
  with(tbl, roc(ClassLabel, sol.Pvalue))
  with(tbl, roc(ClassLabel ~ sol.Pvalue))
  
  # With a formula:
  roc(ClassLabel ~ sol.Pvalue, data=tbl)
  
  # With controls/cases
  roc(controls=tbl$sol.Pvalue[tbl$ClassLabel=="TruePositive"], cases=tbl$sol.Pvalue[tbl$ClassLabel=="TrueNegative"])
  
  # Inverted the levels: "Poor" are now controls and "Good" cases:
  roc(tbl$ClassLabel, tbl$sol.Pvalue,
      levels=c("TruePositive", "TrueNegative"))
  
  # The result was exactly the same because of direction="auto".
  # The following will give an AUC < 0.5:
  roc(tbl$ClassLabel, tbl$sol.Pvalue,
      levels=c("TruePositive", "TrueNegative"), direction="<")
  
  # If we prefer counting in percent:
  roc(tbl$ClassLabel, tbl$sol.Pvalue, percent=TRUE)
  
  # Test the different algorithms:
  roc(tbl$ClassLabel, tbl$sol.Pvalue, algorithm = 1)
  roc(tbl$ClassLabel, tbl$sol.Pvalue, algorithm = 2)
  roc(tbl$ClassLabel, tbl$sol.Pvalue, algorithm = 3)
  if (require(microbenchmark)) {
    roc(tbl$ClassLabel, tbl$sol.Pvalue, algorithm = 0)
  }
  
  # Plot and CI (see plot.roc and ci for more options):
  roc(tbl$ClassLabel, tbl$sol.Pvalue,
      percent=TRUE, plot=TRUE, ci=TRUE)
  
  # Smoothed ROC curve
  roc(tbl$ClassLabel, tbl$sol.Pvalue, smooth=TRUE)
  # this is not identical to
  smooth(roc(tbl$ClassLabel, tbl$sol.Pvalue))
  # because in the latter case, the returned object contains no AUC
  
}

PrepDataFromRoc <- function(x) {
  a <- data.frame(x[1])
  b <- data.frame(x[2])
  
  a$ClassLabel <-rep('unknown',nrow(a))
  b$ClassLabel <-rep('unknown',nrow(b))
  
  
  for(i in 1:nrow(a)){  
    a[i,]$ClassLabel <- 'TruePositive'  
  }
  for(i in 1:nrow(b)){  
    b[i,]$ClassLabel <- 'TrueNegative'  
  }
  
  setnames(b, colnames(b), colnames(a))
  tbl <- do.call(rbind, list(a, b))
  return(tbl)
}

MyRocR <- function(x) {
  
  # Remove NAs
  tmp <- x[!is.na(x$Pvalue),]
  tmp$Pvalue <- -log(tmp$Pvalue, 10)
  
  pred <- data.frame(tmp$Pvalue)
  labels <- data.frame(tmp$ClassLabel)
      
  pred <- prediction(pred, labels)
  perf <- performance(pred, 'tpr', 'fpr')
  
  as.numeric(performance(pred,"auc")@y.values)
  
   return(list(
     Perf=perf,
     AUC=as.numeric(performance(pred,"auc")@y.values)))
  
}


Log2Ratio <- function(x, peptideTypes) {
  
  minArea <- min(x$PeakArea, na.rm=TRUE)
  fileNames <- unique(x$SpectrumFileName)
  conditions <- unique(x$Condition)
  
  res <- list()
  # Create results table for each peptideType, store in a list
  for (i in 1:length(peptideTypes)) {
    tmpDataTable <- data.table(
      Peptide = character(),
      ComparisonCondition = character(),
      Estimate = numeric(),
      Pvalue = numeric())
    res <- lappend(res, tmpDataTable, peptideTypes[i])
    
    db <- paste("OS=",peptideTypes[i],sep="")
    
    #Need to get regex for human and ecoli peptide list
    #another way to subset with regex
    #pepTypeDt <- x[grep(db, x$ParentProteinDescriptions)]
    
    # Subset x for the specifice peptideType
    pepTypeDt <- x[like(ParentProteinDescriptions, peptideTypes[i])]    
    
    setkey(pepTypeDt, 'Sequence')
    
    # Get unique list of peptide sequences and associated protein id
    uniquePep <- unique(pepTypeDt$Sequence)
    
    # Subset x for each peptide
    for (p in 1:length(uniquePep)) {
      dt <- pepTypeDt[J(uniquePep[p])]
      
      #another way to select rows
      #dt <- subset(pepTypeDt,Sequence == uniquePep[p])
      
      # Find any missing file names in subsetted data.table
      # missingFileNames <- setdiff(fileNames, unique(dt$SpectrumFileName)) 
      
      missingFileNames <- setdiff(conditions, unique(dt$Condition)) 
      
      if(length(missingFileNames) == 0){
                 
      
      #Create a temporary data.table to store the missing file names
      #missingDt <- pepTypeDt[which(pepTypeDt$SpectrumFileName%in%missingFileNames)]
      #missingDt <- unique(data.table(
      #  SpectrumFileName=missingDt$SpectrumFileName))
      
      #Interpolate / Gap-fill the PeakAreas in missing file name table
      #missingDt$MaxArea = rep(minArea, nrow(missingDt))
      
      
      #dt[is.na(PeakArea)]$PeakArea <- minArea
      #dt <- dt[,list(MaxArea=max(log(PeakArea, 10))),
      #        by=list(Condition, SpectrumFileName)]
      
      #dt[is.na(PeakArea)]$PeakArea <- minArea
      #dt <- dt[,list(MaxArea=max(PeakArea)),
      #         by=list(Condition, SpectrumFileName)]
      
      #dt <- rbind(dt, missingDt)
   
      
      # Check that peptide is found in more than one file,
      #if (length(unique(dt$Condition)) > 1) {
        
        # subset the data to only ge the +2 charge staes
        #dt<- subset(dt,dt$Charge == 2)
        
        dt <- dt[,list(MaxArea=max(PeakArea)),
                 by=list(Charge)]
        
        #get maximum peak area per condition
        dt <- dt[,list(MaxArea=max(PeakArea)),
                 by=list(Condition, SpectrumFileName)]
        
        #get the maximum peak area per condition
        dt <- dt[,list(MaxArea=max(MaxArea)),
                    by=list(Condition)]
               
        if(nrow(dt) > 0){
        
        
        conditionLevels <- dt$Condition
        
        #take log2 ratio of replicates (replicate 1 is 'control')
        #for(k in 1: nrow(dt)){
                   
          currentRatio <- log2((dt[2]$MaxArea/dt[1]$MaxArea))
          
          # Add result row to result data.table in list
          res[[i]] <- rbind(res[[i]],
                            data.table(
                              Peptide=uniquePep[p],
                              ComparisonCondition= d<- conditionLevels[1],
                              Estimate=currentRatio,                              Pvalue=1))
                                                                
        #}
       
      }      
            
     }

    }
  }   
  return(res)
}


AssignDilutionReplicates <- function(x){
  
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, x)
  mydata = data.table(dbGetQuery(con,"SELECT * from ChromatographicComponents"))
  
  idList <- unique(mydata$FileId)
  
  mydata$Dilution <- rep('Unknown', nrow(mydata))
  mydata$Replicate <- rep('Unknown', nrow(mydata))
   
  numberOfDilutions <- 10
  numberOfReplicates <- 4
  index <- 1
  fileID <- c()
  for(d in 1:numberOfDilutions){
    dilution <- paste('Dilution', d, sep='')
    fileIDs <- c()
    for (r in 1:numberOfReplicates) {
      replicate <- paste('Rep', r, sep='')
      fileID <- idList[index]
      mydata[FileId==fileID]$Replicate <- replicate
      fileIDs <- c(fileIDs, fileID)
      index <- index + 1
    }
    mydata[FileId%in%fileIDs]$Dilution <- dilution
  }
  
  return (mydata) 
}

WriteSqlite <-function(databaseFile,name,dataTable){
  
  # open table and get chromatograhic components
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, databaseFile)
  d <- dbReadTable(con, "ChromatographicComponents")
  
  #save another table with additional columns
  dbWriteTable(con,"DilutionTable", dataTable)
  
  
}

SnTtest <- function(x){
  
  res <- list()
  tmpDataTable <- data.table(
  sn = numeric())
  
  conditions = 10
  for(i in 1:conditions){
    c <- paste('Dilution',i,sep = '')
    s <- subset(x,x$Dilution == c)
    res[[i]] <- data.table(s$SignaltoNoise)
    
  }
  return(res)
  
}

PruneTable <-function(x){
  i
}

# BoxPlot <- function(x){
#   for(i in 1: length(x)){
#     data.frame(unlist(r[x])
#   }
#   
# }
box <- function(r){
  bp <- data.frame(unlist(r[1]),unlist(r[2]),unlist(r[3]),unlist(r[4]),unlist(r[5]),unlist(r[6]),unlist(r[7]),unlist(r[9]),unlist(r[10]))
  boxplot(bp)
}


MassAnalyzer <- function(dt,a){
  
  res <-list()
  
  resTable <- data.table(
    PeptideSequence = character(),
    TheoreticalMass = numeric(),
    DetectedMass = numeric(),
    PPMError = numeric(),
    ExpectedRT = numeric(),
    DetectedRT = numeric(),
    RtError = numeric(),
    Area = numeric())  

    
  dt <- data.frame(dt)[,c(1,7,8,10)]   
  dt$Time <- as.numeric(as.character(dt$Time))
  dt$m.z <- as.numeric(as.character(dt$m.z))  
  dt$X130124_dilA_1_01.raw <- as.numeric(as.character(dt$X130124_dilA_1_01.raw))  
  
  for(i in 1:nrow(a)){
    
    cA <- a[i,]
    cA$TheoreticalSequence <- as.character(cA$TheoreticalSequence)
    cA$TheoreticalTargetMass <- as.numeric(cA$TheoreticalTargetMass)
    cA$RetentionTimeAtApex <- as.numeric(cA$RetentionTimeAtApex)
    
    

    for(j in 1: nrow(dt)){
      
      tbl <- MyBS(dt,cA$TheoreticalTargetMass)
      
      
      if(j == 73){
        test2 = 0
      }
      unknown <- dt[j,]
      t1 <- cA$TheoreticalTargetMass
      t2 <- unknown$m.z
      t2Time <- unknown$Time
      rtError <- t2Time - cA$RetentionTimeAtApex
      area <- unknown$X130124_dilA_1_01.raw
      
      ppm <- abs(ComputeMassError(t1,t2))
      
      if(abs(ComputeMassError(t1,t2)) < 10
         & cA$RetentionTimeAtApex <= t2Time + 3.5
         &  cA$RetentionTimeAtApex >= t2Time - 3.5){
        
        t <-data.table(
          PeptideSequence = cA$TheoreticalSequence,
          TheoreticalMass = cA$TheoreticalTargetMass,
          DetectedMass = t2,
          PPM = ppm,
          ExpectedRT = cA$RetentionTimeAtApex,
          DetectedRT = t2Time,
          RtError = rtError,
          Area = area)        
        
        if(nrow(resTable) == 0){
          resTable = t
        }
        else{
          resTable <- rbind(resTable,t)
        }
       
      }
      
    }
  }  
  return(resTable)
}



ComputeMassError <- function(theoreticalMass,experimentalMass){
  
  return (((experimentalMass - theoreticalMass) / theoreticalMass) * 1e6);

}

Findtolerance <- function(tm,ppm){
  return (abs((ppm*tm/1e6) + tm)-tm)
}


MyBS <-function (dataTable,value){
  
  #sort array 
  dt <- dataTable[order(dataTable$m.z),]
  
  lowBound <- 1
  highBound <- nrow(dt)
  
  mid = 0
  while (lowBound <= highBound)
  {
    mid = floor((lowBound + highBound) / 2)
        
    if(abs(ComputeMassError(value,dt[mid,]$m.z)) < 10){
      return (dt[mid,])
    }      
    
    else if (dt[mid,]$m.z< value)#the element we search is located to the right from the mid point
    {
      lowBound = mid + 1      
    }
    else if (dt[mid,]$m.z > value)#the element we search is located to the left from the mid point
    {
      highBound = mid - 1      
    }
    #at this point low and high bound are equal and we have found the element or
    #arr[mid] is just equal to the value => we have found the searched element
#     else
#     {
#       return (mid)
#     }
  }  
return (-1)#value not found

}

MA <- function(targets,peakList){
  
  resTable <- data.table(
    PeptideSequence = character(),
    TheoreticalMass = numeric(),
    DetectedMass = numeric(),
    PPMError = numeric(),
    ExpectedRT = numeric(),
    DetectedRT = numeric(),
    RtError = numeric(),
    Area = numeric())  
  
  
  hitcount <- 0
  peakList$Time <- as.numeric(as.character(peakList$Time))
  peakList$m.z <- as.numeric(as.character(peakList$m.z))  
  peakList$X130124_dilA_1_01.raw <- as.numeric(as.character(peakList$X130124_dilA_1_01.raw))  
  
  
  for(i in 1:nrow(targets)){
    t <- Findtolerance(targets[i]$TheoreticalTargetMass,10)
    temp <-peakList[m.z < targets[i]$TheoreticalTargetMass + t &
                      m.z > targets[i]$TheoreticalTargetMass - t &
                      Time < targets[i]$RetentionTimeAtApex + 3.5 &
                      Time > targets[i]$RetentionTimeAtApex - 3.5]
    if(nrow(temp) > 0){
            
      for(j in 1:nrow(temp)){     
        
        hitcount <- hitcount + 1
        ppm <- abs(ComputeMassError(targets[i]$TheoreticalTargetMass,temp[j]$m.z))
        rtError <- temp[j]$Time - targets[i]$RetentionTimeAtApex
        
        t <-data.table(
          PeptideSequence = targets[i]$TheoreticalSequence,
          TheoreticalMass = targets[i]$TheoreticalTargetMass,
          DetectedMass = temp[j]$m.z,
          PPM = ppm,
          ExpectedRT = targets[i]$RetentionTimeAtApex,
          DetectedRT = temp[j]$Time,
          RtError = rtError,
          Area = temp[j]$X130124_dilA_1_01.raw)        
        
        if(nrow(resTable) == 0){
          resTable = t
        }
        else{
          resTable <- rbind(resTable,t)
        }
      }
    }
    
  }
  return (resTable)
}

GetUniquePeptides <-function(x){
  
  resTable <- data.table(
    PeptideSequence = character(),
    TheoreticalMass = numeric(),
    DetectedMass = numeric(),
    PPMError = numeric(),
    ExpectedRT = numeric(),
    DetectedRT = numeric(),
    RtError = numeric(),
    Area = numeric(),
    Dilution = character(),
    Replicate = character())  
  
  uniquePeptides <-unique(x$PeptideSequence)
  
  for(i in 1:length(uniquePeptides)){
    
    pep <- x[PeptideSequence == uniquePeptides[i]]
    
    max <- pep[which(pep$Area == max(pep$Area))]
    
    if(nrow(resTable) == 0){
      resTable = max
    }
    else{
      resTable <- rbind(resTable,max)
    }
    
  }
  return (resTable)
}



LogRatioAnalysisTT <- function(x,species){
  
  x <-  x[grep(species, ParentProteinDescriptions)]
  x <- x[MatchConfidence == 3]

  uniquePeptides <-unique(x$Sequence)
  
  resTable <- data.table(
    PeptideSequence = character(),
    DetectedMZ = numeric(),
    Concentration = numeric(),
    AveValue = numeric()) 
  
  for(i in 1:length(uniquePeptides)){
    
    #get all psms for this peptide sequence
    pep <- x[Sequence == uniquePeptides[i]]
    
    pep <- pep[PeakArea != 'NA']
    
    detectedMass <- mean(pep$MassOverCharge)
    
    #get maximum peak area per condition
    pep <- pep[,list(MaxArea=max(PeakArea)),
                by=list(Sequence,SpectrumFileName)]
    
    #get averages for each condition
    #get max for now since there is no gap filling
    dt1 <- pep[grep('_3ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt1$MaxArea))
    
    if(is.na(aveValue) == FALSE){
    t1 <- data.table(
      PeptideSequence = uniquePeptides[i],
      DetectedMZ = detectedMass,
      Concentration = 3,
      AveValue = aveValue)
    
    resTable <- rbind(resTable,t1)
    }
    
    dt2 <- pep[grep('_5ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt2$MaxArea))
    
    if(is.na(aveValue) == FALSE){
    t2 <- data.table(
      PeptideSequence = uniquePeptides[i],
      DetectedMZ = detectedMass,
      Concentration = 5,
      AveValue = aveValue)    
       
    resTable <- rbind(resTable,t2)
    }
    
    
    dt3 <- pep[grep('_7-5ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt3$MaxArea))
    
    if(is.na(aveValue) == FALSE){
    t3 <- data.table(
      PeptideSequence = uniquePeptides[i],
      DetectedMZ = detectedMass,
      Concentration = 7.5,
      AveValue = aveValue)
    
    resTable <- rbind(resTable,t3)
    }
    
    dt4 <- pep[grep('_10ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt4$MaxArea))
    
    if(is.na(aveValue) == FALSE){
    t4 <- data.table(
      PeptideSequence = uniquePeptides[i],
      DetectedMZ = detectedMass,
      Concentration = 10,
      AveValue = aveValue) 
    
    resTable <- rbind(resTable,t4)
    }
    
    dt5 <- pep[grep('_15ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt5$MaxArea))
    
    if(is.na(aveValue) == FALSE){
    t5 <- data.table(
      PeptideSequence = uniquePeptides[i],
      DetectedMZ = detectedMass,
      Concentration = 15,
      AveValue = aveValue) 
    
    resTable <- rbind(resTable,t5)
    }
        
  } 
  
  return (resTable)
  
}

LogRatioAnalysisMinora <- function(x,species){
  
  x <-  x[grep(species, ParentProteinDescriptions)]
  x <- x[MatchConfidence == 3]
  
  uniquePeptides <-unique(x$Sequence)
  
  resTable <- data.table(
    PeptideSequence = character(),
    DetectedMZ = numeric(),
    Concentration = numeric(),
    AveValue = numeric()) 
  
  for(i in 1:length(uniquePeptides)){
    
    #get all psms for this peptide sequence
    pep <- x[Sequence == uniquePeptides[i]]
    
    pep <- pep[Area != 'NA']
    
    detectedMass <- mean(pep$MassOverCharge)
    
    #get maximum peak area per condition
    pep <- pep[,list(MaxArea=max(Area)),
               by=list(Sequence,SpectrumFileName)]
    
    #get averages for each condition
    #get max for now since there is no gap filling
    dt1 <- pep[grep('_3ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt1$MaxArea))
    
    if(is.na(aveValue) == FALSE){
      t1 <- data.table(
        PeptideSequence = uniquePeptides[i],
        DetectedMZ = detectedMass,
        Concentration = 3,
        AveValue = aveValue)
      
      resTable <- rbind(resTable,t1)
    }
    
    dt2 <- pep[grep('_5ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt2$MaxArea))
    
    if(is.na(aveValue) == FALSE){
      t2 <- data.table(
        PeptideSequence = uniquePeptides[i],
        DetectedMZ = detectedMass,
        Concentration = 5,
        AveValue = aveValue)    
      
      resTable <- rbind(resTable,t2)
    }
    
    
    dt3 <- pep[grep('_7-5ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt3$MaxArea))
    
    if(is.na(aveValue) == FALSE){
      t3 <- data.table(
        PeptideSequence = uniquePeptides[i],
        DetectedMZ = detectedMass,
        Concentration = 7.5,
        AveValue = aveValue)
      
      resTable <- rbind(resTable,t3)
    }
    
    dt4 <- pep[grep('_10ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt4$MaxArea))
    
    if(is.na(aveValue) == FALSE){
      t4 <- data.table(
        PeptideSequence = uniquePeptides[i],
        DetectedMZ = detectedMass,
        Concentration = 10,
        AveValue = aveValue) 
      
      resTable <- rbind(resTable,t4)
    }
    
    dt5 <- pep[grep('_15ng_',pep$SpectrumFileName)]
    aveValue <-(max(dt5$MaxArea))
    
    if(is.na(aveValue) == FALSE){
      t5 <- data.table(
        PeptideSequence = uniquePeptides[i],
        DetectedMZ = detectedMass,
        Concentration = 15,
        AveValue = aveValue) 
      
      resTable <- rbind(resTable,t5)
    }
    
  } 
  
  return (resTable)
  
}

ComputeLogRatios <- function(x,d1,d2){
  
  #final table
  resTable <- data.table(
    PeptideSequence = character(),
    Ratio = numeric())
  
  
  #log transform Areas
  #x$AveValue <- log(x$AveValue,10)
  
  #get all unique peptides sequences
  uniquePeptides <-unique(x$PeptideSequence)
  
  for(i in 1:length(uniquePeptides)){
    
      pep <- x[PeptideSequence == uniquePeptides[i]]
      
      p1 <- pep[Concentration == d1]
      p2 <- pep[Concentration == d2]
      
      #if(nrow(p5) > 0 & nrow(p10) > 0 & nrow(p15) > 0){
      if(nrow(p1) > 0 & nrow(p2) > 0){
        resTable <- rbind(resTable,data.table(
          PeptideSequence = uniquePeptides[i],
          Ratio = (p2$AveValue/p1$AveValue)))
      }        
  }
  return (resTable)  
}



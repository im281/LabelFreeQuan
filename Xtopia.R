source('Functions.R')

GetXtopiaMatches <- function(targets,peakList){
    
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
  
  
  hitcount <- 0

   
  for(i in 1:nrow(targets)){
    t <- Findtolerance(targets[i]$TheoreticalTargetMass,10)
    temp <-peakList[Det.M.Z < targets[i]$TheoreticalTargetMass + t &
                      Det.M.Z > targets[i]$TheoreticalTargetMass - t &
                      RetentionTime < targets[i]$RetentionTimeAtApex + 3.5 &
                      RetentionTime > targets[i]$RetentionTimeAtApex - 3.5]
    if(nrow(temp) > 0){
      
      for(j in 1:nrow(temp)){     
        
        hitcount <- hitcount + 1
        ppm <- abs(ComputeMassError(targets[i]$TheoreticalTargetMass,temp[j]$Det.M.Z))
        rtError <- temp[j]$RetentionTime - targets[i]$RetentionTimeAtApex
        
        t <-data.table(
          PeptideSequence = targets[i]$TheoreticalSequence,
          TheoreticalMass = targets[i]$TheoreticalTargetMass,
          DetectedMass = temp[j]$Det.M.Z,
          PPM = ppm,
          ExpectedRT = targets[i]$RetentionTimeAtApex,
          DetectedRT = temp[j]$RetentionTime,
          RtError = rtError,
          Area = temp[j]$Area,
          Dilution = temp[j]$Dilution,
          Replicate = temp[j]$Replicate)        
        
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


Test <-function(targets,x,d){
  
  
  finalTable <- data.table(
    PeptideSequence = "NA",
    TheoreticalMass = 0,
    DetectedMass = 0,
    PPM = 0,
    ExpectedRT = 0,
    DetectedRT = 0,
    RtError = 0,
    Area = 0,
    Dilution = "NA",
    Replicate = "NA")   
   
  
  for(i in 1:length(x)){
    
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
    
    peakList <- data.table(read.csv(x[i]))
    
    dil <-paste('Dilution',d,sep ="")
    rep <-paste('Rep',i,sep="")
    tp <- targets[Dilution == dil]
    tp <- tp[Replicate == rep]  
    peakList$Dilution <- rep(dil, nrow(peakList))
    peakList$Replicate <- rep(rep, nrow(peakList))
    
    dil <-""
    rep <-""
        
    hitcount <- 0
       
    for(i in 1:nrow(tp)){
      tol <- Findtolerance(tp[i]$TheoreticalTargetMass,10)
      temp <-peakList[Det.M.Z < tp[i]$TheoreticalTargetMass + tol &
                        Det.M.Z > tp[i]$TheoreticalTargetMass - tol &
                        RetentionTime < tp[i]$RetentionTimeAtApex + 3.5 &
                        RetentionTime > tp[i]$RetentionTimeAtApex - 3.5]
      
      if(nrow(temp) == 0){
        t <-data.table(
          PeptideSequence = tp[i]$TheoreticalSequence,
          TheoreticalMass = tp[i]$TheoreticalTargetMass,
          DetectedMass = 0,
          PPM = 0,
          ExpectedRT = tp[i]$RetentionTimeAtApex,
          DetectedRT = 0,
          RtError = 0,
          Area = 0,
          Dilution = tp[i]$Dilution,
          Replicate = tp[i]$Replicate)        
        
        if(nrow(resTable) == 0){
          resTable = t
        }
        else{
          resTable <- rbind(resTable,t)
        }
      }

      if(nrow(temp) > 0){
        
        for(j in 1:nrow(temp)){     
          
          hitcount <- hitcount + 1
          ppm <- abs(ComputeMassError(tp[i]$TheoreticalTargetMass,temp[j]$Det.M.Z))
          rtError <- temp[j]$RetentionTime - tp[i]$RetentionTimeAtApex
          
          t <-data.table(
            PeptideSequence = tp[i]$TheoreticalSequence,
            TheoreticalMass = tp[i]$TheoreticalTargetMass,
            DetectedMass = temp[j]$Det.M.Z,
            PPM = ppm,
            ExpectedRT = tp[i]$RetentionTimeAtApex,
            DetectedRT = temp[j]$RetentionTime,
            RtError = rtError,
            Area = temp[j]$Area,
            Dilution = tp[i]$Dilution,
            Replicate = tp[i]$Replicate)        
          
          if(nrow(resTable) == 0){
            resTable = t
          }
          else{
            resTable <- rbind(resTable,t)
          }
        }
      }      
    } 
    finalTable <- rbind(finalTable,GetUniquePeptides(resTable))
  }
  return (finalTable)
}
source('Functions.R')

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


MATest <- function(targets,x,d){
  
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
  
  col <- 1
  if(d == 1){
    col <- 10
  }  
  if(d == 2){
    col <- 14
  }
  if(d == 3){
    col <- 20
  } 
  if(d == 4){
    col <- 28
  }
  if(d == 5){
    col <- 35
  }
  if(d == 6){
    col <- 39
  }
  if(d == 7){
    col <- 43
  }
  if(d == 8){
    col <- 47
  }
  if(d == 9){
    col <- 51
  }
  if(d == 10){
    col <- 55
  }
   
  for(i in 1:4){    
    
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
             
    peakList <- data.frame(x)[,c(1,7,8,col)]
    peakList <- data.table(peakList)  
    colnames(peakList)[4] <- 'Area'
    peakList <- peakList[Area > 0]
    peakList$Time <- as.numeric(as.character(peakList$Time))
    peakList$m.z <- as.numeric(as.character(peakList$m.z))   
    peakList$Area <- as.numeric(as.character(peakList$Area)) 
    
        
    dil <-paste('Dilution',d,sep ="")
    rep <-paste('Rep',i,sep="")
    tp <- targets[Dilution == dil]
    tp <- tp[Replicate == rep]   
    peakList$Dilution <- rep(dil, nrow(peakList))
    peakList$Replicate <- rep(rep, nrow(peakList))
    
    dil <-""
    rep <-""
    col <- col + 1
    
    for(i in 1:nrow(tp)){ 
      tol <- Findtolerance(tp[i]$TheoreticalTargetMass,10)
      temp <-peakList[m.z < tp[i]$TheoreticalTargetMass + tol &
                        m.z > tp[i]$TheoreticalTargetMass - tol &
                        Time < tp[i]$RetentionTimeAtApex + 3.5 &
                        Time > tp[i]$RetentionTimeAtApex - 3.5]
      
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
          
          ppm <- abs(ComputeMassError(tp[i]$TheoreticalTargetMass,temp[j]$m.z))
          rtError <- temp[j]$Time - tp[i]$RetentionTimeAtApex
          
          t <-data.table(
            PeptideSequence = tp[i]$TheoreticalSequence,
            TheoreticalMass = tp[i]$TheoreticalTargetMass,
            DetectedMass = temp[j]$m.z,
            PPM = ppm,
            ExpectedRT = tp[i]$RetentionTimeAtApex,
            DetectedRT = temp[j]$Time,
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
  return(finalTable)
}

MALogRatioTest <- function(targets,x){
  
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
   
  peakList <- data.frame(ma)[,c(1,7,8,9,10,11,12,13,14,15)] 
  peakList <- data.table(peakList) 
  peakList$Time <- as.numeric(as.character(peakList$Time))
  peakList$m.z <- as.numeric(as.character(peakList$m.z)) 
  peakList$MS.Area <- as.numeric(as.character(peakList$MS.Area)) 
  peakList$MS.Area.1 <- as.numeric(as.character(peakList$MS.Area.1 )) 
  peakList$MS.Area.2 <- as.numeric(as.character(peakList$MS.Area.2 )) 
  peakList$MS.Area.3 <- as.numeric(as.character(peakList$MS.Area.3 )) 
  peakList$MS.Area.4 <- as.numeric(as.character(peakList$MS.Area.4 )) 
  peakList$MS.Area.5 <- as.numeric(as.character(peakList$MS.Area.5 )) 
  
  
  for(i in 1:nrow(targets)){
    
    tol <- Findtolerance(targets[i]$MassOverCharge,10)
    temp <-peakList[m.z < targets[i]$MassOverCharge + tol &
                      m.z > targets[i]$MassOverCharge - tol &
                      Time < targets[i]$ApexRetentionTime + 2 &
                      Time > targets[i]$ApexRetentionTime - 2]
  }
}



  
  
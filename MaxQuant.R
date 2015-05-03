
#Modify data table
# mq$Dilution <- rep('Unknown', nrow(mq))
# mq$Replicate <- rep('Unknown', nrow(mq))
# 
# mq[grep('_1_', mq$Rawfile)]$Dilution <- 'Dilution1'
# mq[grep('_2_', mq$Rawfile)]$Dilution <- 'Dilution2'
# mq[grep('_3_', mq$Rawfile)]$Dilution <- 'Dilution3'
# mq[grep('_4_', mq$Rawfile)]$Dilution <- 'Dilution4'
# mq[grep('_5_', mq$Rawfile)]$Dilution <- 'Dilution5'
# mq[grep('_6_', mq$Rawfile)]$Dilution <- 'Dilution6'
# mq[grep('_7_', mq$Rawfile)]$Dilution <- 'Dilution7'
# mq[grep('_8_', mq$Rawfile)]$Dilution <- 'Dilution8'
# mq[grep('_9_', mq$Rawfile)]$Dilution <- 'Dilution9'
# mq[grep('_10_', mq$Rawfile)]$Dilution <- 'Dilution10'
# mq[grep('_11_', mq$Rawfile)]$Dilution <- 'Dilution11'
# mq[grep('_12_', mq$Rawfile)]$Dilution <- 'Dilution12'
# 
# mq[grep('_01', mq$Rawfile)]$Replicate <- 'Rep1'
# mq[grep('_02', mq$Rawfile)]$Replicate <- 'Rep2'
# mq[grep('_03', mq$Rawfile)]$Replicate <- 'Rep3'
# mq[grep('_04', mq$Rawfile)]$Replicate <- 'Rep4'
# mq[grep('_05', mq$Rawfile)]$Replicate <- 'Rep5'
# mq[grep('_06', mq$Rawfile)]$Replicate <- 'Rep6'
# mq[grep('_07', mq$Rawfile)]$Replicate <- 'Rep7'




MaxQuantSearch <-function(targets,x,d){
  
  
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
    
    peakList <- x
    
    dil <-paste('Dilution',d,sep ="")
    rep <-paste('Rep',i,sep="")
    tp <- targets[Dilution == dil]
    tp <- tp[Replicate == rep]  
    peakList <- peakList[Dilution == dil]
    peakList <- peakList[Replicate == rep]    

    rep <- paste('Rep',i + 1,sep ="")
    
    for(i in 1:nrow(tp)){
      tol <- Findtolerance(tp[i]$TheoreticalTargetMass,10)
      temp <-peakList[Uncalibratedmz < tp[i]$TheoreticalTargetMass + tol &
                        Uncalibratedmz > tp[i]$TheoreticalTargetMass - tol &
                        RT < tp[i]$RetentionTimeAtApex + 3.5 &
                        RT > tp[i]$RetentionTimeAtApex - 3.5]
      
      if(nrow(temp) == 0){
        t  <-data.table(
          PeptideSequence = tp[i]$TheoreticalSequence,
          TheoreticalMass = tp[i]$TheoreticalTargetMass,
          DetectedMass = temp[j]$Uncalibratedmz,
          PPM = ppm,
          ExpectedRT = tp[i]$RetentionTimeAtApex,
          DetectedRT = temp[j]$RT,
          RtError = rtError,
          Area = temp[j]$Intensity,
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
               
          ppm <- abs(ComputeMassError(tp[i]$TheoreticalTargetMass,temp[j]$Uncalibratedmz))
          rtError <- temp[j]$RT - tp[i]$RetentionTimeAtApex
          
          t <-data.table(
            PeptideSequence = tp[i]$TheoreticalSequence,
            TheoreticalMass = tp[i]$TheoreticalTargetMass,
            DetectedMass = temp[j]$Uncalibratedmz,
            PPM = ppm,
            ExpectedRT = tp[i]$RetentionTimeAtApex,
            DetectedRT = temp[j]$RT,
            RtError = rtError,
            Area = temp[j]$Intensity,
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
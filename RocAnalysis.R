QuanAnova <- function(x,dil1,dil2) {
  
  dilution1 <- paste('Dilution',dil1,sep ='')
  dilution2 <- paste('Dilution',dil2,sep ='')
  
  minArea <- min(x$Area, na.rm=TRUE)
  #fileNames <- unique(x$SpectrumFileName)
  
  res <- list()
  
  tmpDataTable <- data.table(
    Peptide = character(),
    ComparisonCondition = character(),
    Estimate = numeric(),
    Pvalue = numeric(),  
    ClassLabel = character())
  
  
  # Get unique list of peptide sequences and associated protein id
  uniquePep <- as.character(unique(x$PeptideSequence)) 
  
  # Subset x for each peptide
  for (p in 1:length(uniquePep)) {
    dt <- x[PeptideSequence == uniquePep[p]]   
    
    #ignore annotations that are not annotated in all replicates
    if(nrow(dt) < 8){
      # Add result row to result data.table in list
      tmpDataTable <- rbind(tmpDataTable,
                            data.table(
                              Peptide=uniquePep[p],
                              ComparisonCondition = 'Dilution',
                              Estimate= 0,
                              Pvalue=1,
                              ClassLabel = 'unknown'))
    }       
    else{
      
      dt$Area <- log(dt$Area,10)
      
      temp <- dt[Area =='-Inf']
      if(nrow(temp) == 0){
        
        if(nrow(dt) > 1){
          #Perform t test
          a1 <- dt[Dilution == dilution1]
          a2 <- dt[Dilution == dilution2]
          r <- t.test(a1$Area,a2$Area)
          
          #get ratio          
          ave1 <- mean(dt[Dilution == dilution1]$Area)
          ave2 <- mean(dt[Dilution == dilution2]$Area)
          
          # Add result row to result data.table in list
          tmpDataTable <- rbind(tmpDataTable,
                                data.table(
                                  Peptide=uniquePep[p],
                                  ComparisonCondition= 'Dilution',
                                  Estimate=(10^ave1)/(10^ave2),
                                  Pvalue= (as.numeric(r[3])),
                                  ClassLabel = dt[1]$Classlabel))
                                
        }
      }
      # if a row has a Inf it means that it was undetected and therefore we
      else{
        # Add result row to result data.table in list
        tmpDataTable <- rbind(tmpDataTable,
                              data.table(
                                Peptide=uniquePep[p],
                                ComparisonCondition= 'Dilution',
                                Estimate= 0,
                                Pvalue=1,
                                ClassLabel = dt[1]$Classlabel))
      }
    }
  }
  
  return(tmpDataTable)
}


RunRoc <- function(detector){
  par(mfrow=c(3,3))
  aucList <-list()
  d = ''
  tp =''
  for(i in 1:9){
    
    if(i == 2 |i == 3){
     #do nothing
    }
    else{
      if(detector == 'MassAnalyzer'){
        d = "MA"
        tp1 <- paste(d,"-Potato-f",i,"-allreps-TP.csv",sep = "")
        tp2 <- paste(d,"-Potato-f",i+1,"-allreps-TP.csv",sep = "")
        tn1 <- paste(d,"-Potato-f",i,"-allreps-TN.csv",sep = "")
        tn2 <- paste(d,"-Potato-f",i+1,"-allreps-TN.csv",sep = "")
      }
      if(detector == 'Xtopia'){
        d = "TT"
        tp1 <- paste(d,"-Potato-f",i,"-allreps-TP2.csv",sep = "")
        tp2 <- paste(d,"-Potato-f",i+1,"-allreps-TP2.csv",sep = "")
        tn1 <- paste(d,"-Potato-f",i,"-allreps-TN1.csv",sep = "")
        tn2 <- paste(d,"-Potato-f",i+1,"-allreps-TN1.csv",sep = "")
      }
      if(detector == 'MaxQuant'){
        d = "MaxQuant"
        tp1 <- paste(d,"-Potato-f",i,"-allreps-TP-Nomatchbetween.csv",sep = "")
        tp2 <- paste(d,"-Potato-f",i+1,"-allreps-TP-Nomatchbetween.csv",sep = "")
        tn1 <- paste(d,"-Potato-f",i,"-allreps-TN-Nomatchbetween.csv",sep = "")
        tn2 <- paste(d,"-Potato-f",i+1,"-allreps-TN-Nomatchbetween.csv",sep = "")
      }
      if(detector == 'Minora'){
        d = "Minora"
      }
 
    #Open all detected annotations
    ############################################################
    tp1 <- data.table(read.csv(tp1))
    tp1$Classlabel<- rep('TruePositive', nrow(tp1))
    
    tp2 <- data.table(read.csv(tp2))
    tp2$Classlabel<- rep('TruePositive', nrow(tp2))
    
    tn1 <- data.table(read.csv(tn1))
    tn1$Classlabel<- rep('TrueNegative', nrow(tn1))
    
    tn2 <- data.table(read.csv(tn2))
    tn2$Classlabel<- rep('TrueNegative', nrow(tn2))
    ############################################################
    
    #Perform ROC analysis Treetop
    ############################################################
    
    finaltp <- rbind(tp1,tp2)
    finaltn <- rbind(tn1,tn2)
    finalAnalysis <- rbind(finaltp,finaltn)
    res <- QuanAnova(finalAnalysis,i,i + 1)
    final <- res[ClassLabel != 'unknown']
    RocResults <- MyRocR(final)
    RocResults[['AUC']]
    title <- paste("d",i," ","vs"," ","d",i + 1," ","AUC="," ",round(RocResults[['AUC']],2),sep = "")
    plot(RocResults[[1]]@x.values[[1]], RocResults[[1]]@y.values[[1]], 
         type='l', xlab=RocResults[[1]]@x.name, ylab=RocResults[[1]]@y.name,main = title)
    aucList[[i]] <-  RocResults[['AUC']]
    #aucList[[i+1]] <- finalAnalysis
    #aucList[[i+2]] <- res
    ############################################################
    }
  }
  return (aucList)
}



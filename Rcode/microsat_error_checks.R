##############################################################################
##                   Error Checks for Microsat Data in CSV                  ##
##############################################################################


##checks repeat motifs, expected size ranges, missing alleles, 
##out of order genotypes, and identifies rare alleles


#setWD
setwd("")


##import genotypes in CSV file
dataFile<-""
genoDat<-read.csv(dataFile,header=T)


#(if not whole numbers there is an error)

##load functions
mode <- function(v) {
  uniqv <- na.omit(unique(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
errorCheck<-function(data,column1,column2,repeatMotif,expectedMin,expectedMax,rareCount=NA,print=F){
  unique1<-unique(unlist(data[column1]),na.rm=T)
  unique2<-unique(unlist(data[column2]),na.rm=T)

  uniqueGeno<-unique(c(unlist(unique1),unlist(unique2)))
  uniqueGeno<-uniqueGeno[!is.na(uniqueGeno)]
  
  if(is.integer(uniqueGeno)==F){
    print("ERROR: may be text in datafile")
    stop()
  }
  
  motifError<-uniqueGeno[which(((uniqueGeno-mode(unlist(data[column1])))/repeatMotif)-floor((uniqueGeno-mode(unlist(data[column1])))/repeatMotif)!=0)]
  motifError<-sort(as.numeric(c(levels(droplevels(data[sapply(data[column1], function(x) x %in% motifError),]$Sample)),levels(droplevels(data[sapply(data[column2], function(x) x %in%motifError),]$Sample)))))
  
  outsideRange<-uniqueGeno[c(which(uniqueGeno < expectedMin),which(uniqueGeno > expectedMax))]
  outsideRange<-sort(as.numeric(c(levels(droplevels(data[sapply(data[column1], function(x) x %in% outsideRange),]$Sample)),levels(droplevels(data[sapply(data[column2], function(x) x %in%outsideRange),]$Sample)))))
  
  counts<-as.data.frame(table(c(as.numeric(unlist(data[column1])),as.numeric(unlist(data[column2])))))
  counts<-counts[order(counts$Freq),]
  if(!is.na(rareCount)){
    counts<-subset(counts,counts$Freq<=rareCount)
  }
  rareGenoVals<-as.numeric(levels(droplevels(counts$Var1)))
  rareGeno<-sort(as.numeric(c(levels(droplevels(data[sapply(data[column1], function(x) x %in% rareGenoVals),]$Sample)),levels(droplevels(data[sapply(data[column2], function(x) x %in%rareGenoVals),]$Sample)))))
  
  alleleOrder<-as.numeric(c(levels(droplevels(data[sapply(data[column1], function(x) x > data[column2]),]$Sample))))
  
  missingAllele<-c(levels(droplevels(data[sapply(data[column1], function(x) is.na(x) != is.na(data[column2])),]$Sample)))
  
  
  if(print==T){
    print("Genotypes:")
    print(sort(as.integer(uniqueGeno)))
    if(length(motifError)>0){
      print("Potential Motif Error, Sample #:")
      print(motifError)
    }
    else{
      print("No Motif Errors Found")
    }
    if(length(outsideRange)>0){
      print("Genotype Outside Expected Range, Sample #:")
      print(outsideRange)
    }
    else{
      print("No Genotypes Outside Expected Range")
    }
    if(length(rareGeno)>0){
      print("Rare Genotypes Found, Sample #:")
      print(outsideRange)
    }
    else{
      print("No Rare Genotypes")
    }
    if(length(alleleOrder)>0){
      print("Allele 2 > Allele 1, Sample #:")
      print(alleleOrder)
    }
    else{
      print("No Allele Order Errors Found")
    }
    if(length(missingAllele)>0){
      print("Missing Single Alelle, Sample #:")
      print(missingAllele)
    }
    else{
      print("No Missing Alleles (ploidy issues) Found")
    }
  }
  returnList<-list(motifError,outsideRange,rareGeno,alleleOrder,missingAllele)
  names(returnList)<-c("motifError","outsideRange","rareGeno","alleleOrder","missingAllele")
  return(returnList)
}

#run error check for all loci-- could also set up in loop for larger number of loci
C113<-errorCheck(data=genoDat,column1 = "Sfo.C113a",column2 = "Sfo.C113b",repeatMotif = 3,expectedMin = 125,expectedMax = 170,rareCount=2)
D75<-errorCheck(data=genoDat,column1 = "Sfo.D75a",column2 = "Sfo.D75b",repeatMotif = 4,expectedMin = 165,expectedMax = 250,rareCount=2)
C88<-errorCheck(data=genoDat,column1 = "Sfo.C88a",column2 = "Sfo.C88b",repeatMotif = 3,expectedMin = 170,expectedMax = 205,rareCount=2)
D100<-errorCheck(data=genoDat,column1 = "Sfo.D100a",column2 = "Sfo.D100b",repeatMotif = 4,expectedMin = 200,expectedMax = 276,rareCount=2)
C24<-errorCheck(data=genoDat,column1 = "Sfo.C24a",column2 = "Sfo.C24b",repeatMotif = 3,expectedMin = 110,expectedMax = 190,rareCount=2)
C129<-errorCheck(data=genoDat,column1 = "Sfo.C129a",column2 = "Sfo.C129b",repeatMotif = 3,expectedMin = 215,expectedMax = 270,rareCount=2)
C115<-errorCheck(data=genoDat,column1 = "Sfo.C115a",column2 = "Sfo.C115b",repeatMotif = 2,expectedMin = 225,expectedMax = 376,rareCount=2)
D237<-errorCheck(data=genoDat,column1 = "Sfo.D237a",column2 = "Sfo.D237b",repeatMotif = 4,expectedMin = 270,expectedMax = 494,rareCount=2)

#join together in list
errors<-list(C113,D75,C88,D100,C24,C129,C115,D237)
#create list of loci names
loci<-c("C113","D75","C88","D100","C24","C129","C115","D237")

#create blank lists for storing error data
errorLocus<-list()
errorType<-list()
errorSample<-list()

#loop through lists and save samples with potential errors
counter=1
for(i in 1:length(errors)){
  ##repeat motifs
  numMotifErrors<-length(errors[[i]]$motifError)
  if(numMotifErrors>0){
    errorLocus[counter:(counter+numMotifErrors-1)]<-rep(loci[i],numMotifErrors)
    errorType[counter:(counter+numMotifErrors-1)]<-rep("motif",numMotifErrors)
    errorSample[counter:(counter+numMotifErrors-1)]<-errors[[i]]$motifError
    counter=counter+numMotifErrors
  }
  ##expected size range
  numRangeErrors<-length(errors[[i]]$outsideRange)
  if(numRangeErrors>0){
    errorLocus[counter:(counter+numRangeErrors-1)]<-rep(loci[i],numRangeErrors)
    errorType[counter:(counter+numRangeErrors-1)]<-rep("range",numRangeErrors)
    errorSample[counter:(counter+numRangeErrors-1)]<-errors[[i]]$outsideRange
    counter=counter+numRangeErrors
  }
  ##rare allele (may not be actual error, but potentially worth double checking)
  numRareErrors<-length(errors[[i]]$rareGeno)
  if(numRareErrors>0){
    errorLocus[counter:(counter+numRareErrors-1)]<-rep(loci[i],numRareErrors)
    errorType[counter:(counter+numRareErrors-1)]<-rep("rare",numRareErrors)
    errorSample[counter:(counter+numRareErrors-1)]<-errors[[i]]$rareGeno
    counter=counter+numRareErrors
  }
  ##allele order- column 1 should be smaller/lower
  numAlleleOrderErrors<-length(errors[[i]]$alleleOrder)
  if(numAlleleOrderErrors>0){
    errorLocus[counter:(counter+numAlleleOrderErrors-1)]<-rep(loci[i],numAlleleOrderErrors)
    errorType[counter:(counter+numAlleleOrderErrors-1)]<-rep("alleleOrder",numAlleleOrderErrors)
    errorSample[counter:(counter+numAlleleOrderErrors-1)]<-errors[[i]]$alleleOrder
    counter=counter+numAlleleOrderErrors
  }
  ##missing allele call?
  missingAlleleErrors<-length(errors[[i]]$missingAllele)
  if(missingAlleleErrors>0){
    errorLocus[counter:(counter+missingAlleleErrors-1)]<-rep(loci[i],missingAlleleErrors)
    errorType[counter:(counter+missingAlleleErrors-1)]<-rep("missingAlelle",missingAlleleErrors)
    errorSample[counter:(counter+missingAlleleErrors-1)]<-errors[[i]]$missingAllele
    counter=counter+missingAlleleErrors
  }
}

#create dataframe and export errors
errorOutDat<-data.frame(unlist(errorLocus),unlist(errorType),unlist(errorSample))
write.csv(errorOutDat,"errorData.csv")

#identify samples in original genotype file, export data
genoDat$errorCheck<-sapply(genoDat$Sample, function(x) if(x %in% errorOutDat$unlist.errorSample.){return(1)}else{return(0)})
write.csv(genoDat,"data_with_potentialErrors.csv")


##subset to error check data, output to new file
errorSamp<-sort(as.numeric(levels(droplevels(genoDat[sapply(genoDat$Purpose, function(x) grepl("error", x)),]$Sample))))
errorDat<-subset(genoDat,genoDat$Sample %in% errorSamp)
write.csv(errorDat,"errorCheckData.csv")
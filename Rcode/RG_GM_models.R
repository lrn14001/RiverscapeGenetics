setwd("C:/Users/lrn14001/Desktop/BKT/fall2017Analyses/LG_focal/FINAL_ANALYSIS/")

##load packages
library("gstudio")
library("adegenet")
library("hierfstat")
library("lme4")
library("scales")
library("doBy")

##load custom functions
readFROMcolumn<-function(df,locus.columns){
  #copied from gstudio package
  end.col.meta.data <- min(locus.columns) - 1
  locus.columns <- seq(min(locus.columns), (max(locus.columns) - 1), by = 2)
  ret <- df[, 1:end.col.meta.data]
  for (locCol in locus.columns) {
    alleles <- df[, locCol:(locCol + 1)]
    tmp <- paste0(alleles[,1],":",alleles[,2])
    tmp[which(tmp==":")]<-NA
    locus_name <- names(df)[locCol]
    ret[[locus_name]] <- tmp
  }
  return(ret)
}
mat2listFULL <- function(x,name) {
  list<-as.list(as.data.frame(x))
  IDs<-names(list)
  pair1<-unlist(lapply(IDs,function(x) rep(x,length(IDs))))
  pair2<-rep(IDs,length(IDs))
  val<-unlist(list)
  outList<-data.frame(pair1,pair2,val)
  names(outList)<-c("site1","site2",name)
  return(outList)
}
singlePair2Doublepair<-function(x,site1,site2,vals){
  doubles<-x[c(site2,site1,vals[1:length(vals)])]
  names(doubles)<-c("site1","site2",vals[1:length(vals)])
  dat2bind<-doubles[c("site1","site2",vals[1:length(vals)])]
  outDat<-rbind(x,dat2bind)
  return(outDat)
}
logTrans<-function(dataframe,vars2trans){
  for(v in vars2trans){
    col2Trans<-which(names(dataframe)==v)
    dataframe[,col2Trans]<-log(dataframe[,col2Trans])
    dataframe[,col2Trans][which(!is.finite(dataframe[,col2Trans]))] <- 0
  }
  return(dataframe)
}
find.r2 <- function(lme.object){
  #modified from Dileo et al. 2014 see cornus data files
  yhat <- unlist(predict( lme.object ))
  #modifed for varying genetic variable
  obs <- unlist(lme.object@resp$`.->y`)
  f <- lm( obs ~ yhat )
  summary(f)$adj
}
SHEDS_mean<-function(data,variable,distance){
  dataTemp<-data
  dataTemp$SumVar<-(dataTemp[variable]*dataTemp[distance])
  resis_paths<-aggregate(dataTemp$SumVar, by=list(Category=dataTemp$Name), FUN=sum)
  dist_paths<-aggregate(dataTemp[distance], by=list(Category=dataTemp$Name), FUN=sum)
  paths<-merge(resis_paths, dist_paths, by="Category")
  paths[paste("avg_",variable,sep="")]<-paths[2]/paths[3]
  #pairID<-c("site1","site2")
  #siteIDs<-data[c(2:4)]
  #paths[pairID] <- lapply(pairID, function(x) siteIDs[[x]][match(paths$Name,siteIDs$Name)])
  #names(paths)<-c("Pair","SumVar","dist",paste("Avg_",variable,sep=""),"site1","site2")
  names(paths)[1]<-"Name"
  return(paths)
}


#####################################################################
#                        read in data
#####################################################################
#read in model data
SRDAT<-read.csv("SR_modelDAT.csv")
FRDAT<-read.csv("FR_modelDAT.csv")

datasets=list(FRDAT,SRDAT)
datasetNames=c("WF","EF") #labeled as "WF", western focal area and "EF", eastern focal area in datasets



#####################################################################
#                         set model parameters candidate 
#####################################################################

#specify variables to log transfrom
vars2trans<-c("FST","b.avg_impervious","b.avg_slope_pcnt","b.avg_surfcoarse_upstrm",
              "a.alloffnet","a.allonnet","a.AreaSqKM","a.may_prcp_mm","a.impervious","a.surfCoarse",
              "a.dsAREA")


#specify variables in each model
natural<-c("prcp","slope","surfcoarse","Area","AREA")
anthro<-c("imperv","canopy")
land<-c("imperv","canopy","slope")
water<-c("prcp","net","ann_tmax","surfcoarse","Area","AREA")

#model names
modelNames<-c("distance","dams","dambin","betweenSites","atSite","global",
              "natural.global","natural.between","natural.site",
              "anthro.global","anthro.between","anthro.site",
              "land.global","land.between","land.site",
              "water.global","water.between","water.site")


#specify distance max threshold
distMax<-30000


##remove sites (following post hoc residual evaluation)
removeSites<-c(103,104,94,110,109,98,83)
removeSites<-NA


#####################################################################
#                         run candidate models
#####################################################################
#create blank lists to store model outputs
dataSet<-variable<-level<-C_modRsq<-C_modAIC<-C_modBIC<-C_modLL<-formulas<-list()

#start counter
counter=1

#loop through both datasets
for(d in 1:length(datasets)){
  #create temporary DF
  tempData<-datasets[[d]]
  #remove sites, if applicable
  tempData<-tempData[which(!tempData$site1 %in% removeSites),]
  tempData<-tempData[which(!tempData$site2 %in% removeSites),]
  #subset to distance threshold
  tempData<-subset(tempData,tempData$distance<distMax)
  #log transform variables
  dataTrans<-logTrans(tempData,names(tempData)[grepl(paste(vars2trans,collapse="|"),x = names(tempData))])
  dataTrans$distance<-scale(dataTrans$distance)
  #make doubles of all datapoints
  dataTrans<-singlePair2Doublepair(dataTrans,"site1","site2",names(dataTrans)[c(2,4:length(dataTrans))])
  
  #get variable names based on models specified before loop (allows for differences in spatial level)
  varNames<-names(dataTrans[c(10:length(dataTrans))])
  natVars<-varNames[grepl(paste(natural,collapse="|"),x = varNames)]
  anthroVars<-varNames[grepl(paste(anthro,collapse="|"),x = varNames)]
  landVars<-varNames[grepl(paste(land,collapse="|"),x = varNames)]
  waterVars<-varNames[grepl(paste(water,collapse="|"),x = varNames)]
  #create models
  models<-c("distance","dams","damsbin",
            paste(names(dataTrans)[grepl("^b.",x = names(dataTrans))],collapse="+"),
            paste(names(dataTrans)[grepl("^a.",x = names(dataTrans))],collapse="+"),
            paste(varNames,collapse="+"),
            paste(natVars,collapse="+"),
            paste(natVars[grepl("^b.",x = natVars)],collapse="+"),
            paste(natVars[grepl("^a.",x = natVars)],collapse="+"),
            paste(anthroVars,collapse="+"),
            paste(anthroVars[grepl("^b.",x = anthroVars)],collapse="+"),
            paste(anthroVars[grepl("^a.",x = anthroVars)],collapse="+"),
            paste(landVars,collapse="+"),
            paste(landVars[grepl("^b.",x = landVars)],collapse="+"),
            paste(landVars[grepl("^a.",x = landVars)],collapse="+"),
            paste(waterVars,collapse="+"),
            paste(waterVars[grepl("^b.",x = waterVars)],collapse="+"),
            paste(waterVars[grepl("^a.",x = waterVars)],collapse="+"))
  #analyze model
  for(mod in 1:length(modelNames)){
      C_form<-as.formula(paste("FST~",paste(models[mod],"+(1|site1)",collapse="+")))
      C_mod<-lmer(C_form,dataTrans,REML=T)
      #store values
      dataSet[counter]<-datasetNames[d]
      variable[counter]<-modelNames[mod]
      C_modRsq[counter]<-find.r2(C_mod)
      C_modAIC[counter]<-AIC(C_mod)
      C_modBIC[counter]<-BIC(C_mod)
      C_modLL[counter]<-logLik(C_mod)
      formulas[counter]<-as.character(C_form)[3]
      
      counter=counter+1
  }
}

#store results in DF
modelRESULTS<-data.frame(unlist(dataSet),unlist(variable),
                         unlist(C_modRsq),unlist(C_modAIC),unlist(C_modBIC),unlist(C_modLL),unlist(formulas))
colnames(modelRESULTS)<-c("dataset","model","C_Rsq","C_AIC","C_BIC","C_LL","formula")



#####################################################################
#                         model selection weights 
#####################################################################
#create generic ID
modelRESULTS$modID<-1:nrow(modelRESULTS)
#specify selection criteria
modSelecValues<-c("C_AIC","C_BIC")
#create blank table
weightTable<-as.data.frame(matrix(nrow = nrow(modelRESULTS),ncol=length(modSelecValues)+1))
names(weightTable)<-c("C_AIC_wt","C_BIC_wt","modID")
#start counter
counter=1
#calculate weigths
for(i in datasetNames){
  modelRESULTSsub<-subset(modelRESULTS,modelRESULTS$dataset==i)
  for(l in 1:length(modSelecValues)){
    vals<-modelRESULTSsub[modSelecValues[l]]
    valsNoNA<-na.omit(vals)
    min<-min(na.exclude(vals))
    delta<-vals-min
    RL<-exp(-0.5*(delta))
    sumRL<-sum(na.exclude(RL))
    weight<-RL/sumRL
    ranks<-rank(-weight)
    ranks[which(is.na(vals))]<-NA
    weightTable[counter:(counter+nrow(modelRESULTSsub)-1),l]<-unlist(round(weight,digits=3))
  }
  weightTable[counter:(counter+nrow(modelRESULTSsub)-1),ncol(weightTable)]<-modelRESULTSsub$modID
  counter=counter+nrow(modelRESULTSsub)
}

##join with output data
outDATA<-merge(modelRESULTS,weightTable,by="modID")
#subset data
finalTable<-outDATA[which(outDATA$C_BIC_wt>0),c(2,3,8,6,10,4)]
#reorder by BIC weight
finalTable<-finalTable[order(finalTable$C_BIC_wt,decreasing = T),]
#save output
write.csv(finalTable,paste0("LG_model_output_",Sys.Date(),".csv"),row.names = F)


#####################################################################
#                         residual diagnostics 
#####################################################################
#FR model
tempData<-datasets[[1]]
tempData<-tempData[which(!tempData$site1 %in% removeSites),]
tempData<-tempData[which(!tempData$site2 %in% removeSites),]
dataTrans<-tempData
dataTrans<-logTrans(tempData,names(tempData)[grepl(paste(vars2trans,collapse="|"),x = names(tempData))])
#make doubles of all datapoints
dataTrans<-singlePair2Doublepair(dataTrans,"site1","site2",names(dataTrans)[c(2,4:length(dataTrans))])
top_form<-as.formula(FST~b.avg_impervious_upstrm + b.avg_tree_canopy_rip50+(1 | site1))
top_modFR<-lmer(top_form,dataTrans,REML=T)

#model diagnostic plots
library(ggplot2)
plotDatFR<-as.data.frame(cbind(resid(top_modFR),fitted(top_modFR)));names(plotDatFR)<-c("Residuals","Fitted")
ggplot(plotDatFR,aes(x=Fitted,y=Residuals))+geom_point()+theme_classic()+geom_hline(yintercept=0,linetype="dotted")+
  theme(axis.title=element_text(size=24),axis.text=element_text(size=14))
ggplot(plotDatFR,aes(x=Residuals),color="white")+geom_histogram()+theme_classic()+ylab("Count")

#SR model
tempData<-datasets[[2]]
dataTrans<-logTrans(tempData,names(tempData)[grepl(paste(vars2trans,collapse="|"),x = names(tempData))])
dataTrans$distance<-scale(dataTrans$distance)
#make doubles of all datapoints
dataTrans<-singlePair2Doublepair(dataTrans,"site1","site2",names(dataTrans)[c(2,4:length(dataTrans))])
top_form<-as.formula(FST~distance + (1 | site1))
top_modSR<-lmer(top_form,dataTrans,REML=T)

#model diagnostic plots
plotDatSR<-as.data.frame(cbind(resid(top_modSR),fitted(top_modSR)));names(plotDatSR)<-c("Residuals","Fitted")
ggplot(plotDatSR,aes(x=Fitted,y=Residuals))+geom_point()+theme_classic()+geom_hline(yintercept=0,linetype="dotted")+
  theme(axis.title=element_text(size=24),axis.text=element_text(size=14))
ggplot(plotDatSR,aes(x=Residuals),color="white")+geom_histogram()+theme_classic()+ylab("Count")




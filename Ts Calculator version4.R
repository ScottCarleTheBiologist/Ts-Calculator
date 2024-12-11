##Ts Calculator - a script to compile spike wetting test data, filter it, 
#and calculate and export Ts values (sprouting score threshold intercepts)
#for use in data analysis and mapping attempts.

####first, we need to define a few things. Every line within the input area must be run, because all these variables need to be defined.
#designate the thresholds that you want to use.
ThresholdVector<-c(1.1,1.5,2.5,3.5,5.5,7.5)
#next, specify your working directory
setwd("E:/")
#then, read in your spike wetting test file.
#This should be formatted as specified in the documentation at: https://github.com/ScottCarleTheBiologist/Ts-Calculator/tree/main
Filename<-"example data Ts Calculator testing"
incomingfile<-read.csv(paste(Filename,".csv",sep=""))
###the name column is the column that averages will be grouped by
#here, it is hardcoded as the first column in your incoming csv file.
#if you want to use a different column, entry number, for instance, then change the 'namecolumn' variable to the number of the column that has the data you desire for making averages from
namecolumn<-1
###the datacolumns variable should include: Day1, Day2, Day3, Day4, Day5, Day6, Day7
###these are hardcoded to be columns 2 through 8 in your incoming csv file.
#if you want to change these to 7 columns at a different position, you can change those numbers below
#however, if you want to change it to not be 7 columns, then talk with Scott (scott.carle@wsu.edu), he's not guaranteeing that flexibility in this version, but it can be done
datacolumns<-c(2:8)
###here we set the max and min scores that it's possible to have
#these are the floor and ceiling of the analysis
#with our scoring method, we are setting these to 1 and 10, but experimentally, we have considered setting an upper limit at 5.
###these are being made into variables to allow experimenting with your population if you desire.
minscore<-1
maxscore<-10
##if you have errors in your data set, and want to omit specific rows, put those rows in the following line.
RowsToOmit<-c()

################Input area ended.
#####################Everything below requires no manual input
#It should automatically produce the output files in your working directory
#############step 1
slopevector<-vector()
offsetvector<-vector()
#Core loop, several parts to this loop
i<-1
for(i in c(1:nrow(incomingfile)))
{
  ##for the first part of this loop we're just taking 1 row of data at a time, and reformatting it
  PlantObservationsSubset<-as.data.frame(cbind(c(1:length(incomingfile[i,datacolumns])),t(incomingfile[i,datacolumns])))
  colnames(PlantObservationsSubset)<-c("Day","Score")
  NAmaker<-vector()
  #The following allows for omissions from the entry area, if there are 1 or more rows from the input data that the user desires to omit
  #These omissions can help with data anomalies that could otherwise break the script or produce flawed results
  if(i %in% RowsToOmit){
    SlopeValue<-NA
    OffsetValue<-NA
  }
  if(!(i %in% RowsToOmit)){
    #Here inserting logic to translate preliminary 1's and subsequent 10's into NA's for the linear fit
    minrow<-as.numeric(max(PlantObservationsSubset[PlantObservationsSubset$Score==minscore,1]))
    maxrow<-as.numeric(min(PlantObservationsSubset[PlantObservationsSubset$Score>=maxscore,1]))
    if(minrow==-Inf){minrow<-1}
    if(minrow==as.numeric(nrow(PlantObservationsSubset))){minrow<-1}
    if(maxrow==Inf){maxrow<-as.numeric(nrow(PlantObservationsSubset))}
    if(maxrow==1){maxrow<-1}
    if(minrow>1){NAmaker<-c(1:(minrow-1))}
    if(maxrow<as.numeric(nrow(PlantObservationsSubset))){NAmaker<-c(NAmaker,maxrow+1:as.numeric(nrow(PlantObservationsSubset)))}
    #replace things
    PlantObservationsSubset[NAmaker,2]=NA
    #now, we calculate a linear regression for each head
    attach(PlantObservationsSubset)
    linearmodel1<-lm(Score~Day,data=PlantObservationsSubset)
    SlopeValue<-as.numeric(round(linearmodel1$coefficients[2],digits=4))
    OffsetValue<-as.numeric(round(linearmodel1$coefficients[1],digits=4))
    detach(PlantObservationsSubset)
  }
  slopevector<-c(slopevector,SlopeValue)
  offsetvector<-c(offsetvector,OffsetValue)
}
###loop ended, now staple together
LinearResults<-cbind(incomingfile[,namecolumn],slopevector,offsetvector)
colnames(LinearResults)[1]<-names(incomingfile)[namecolumn]
###Now, all your samples, have linear regression formulae calculated for them, and in these Linear Results.

##There is a philosophical conundrum about the next step.
#One can either average the slopes and offsets now, which makes use of data that lie on the floor and ceiling the entire time
#or, one can do what I did in the paper, and omit from the intercept calculations those samples with a slope of 0.
#There is a case to be made for either method.
#here we can export a file that contains all the linear models, if we unhash the next line
#write.csv(LinearResults,file=paste(Filename,"Individual Sample Linear Models.csv"))
####Step 2
##Here we are going to calculate the intercepts of the linear regressions with the various 
InterceptMatrix<-LinearResults[,1]
j<-1
for(j in c(1:length(ThresholdVector)))
{
  InterceptVector<-((ThresholdVector[j] - as.numeric(offsetvector))/as.numeric(slopevector))
  InterceptMatrix<-cbind(InterceptMatrix,InterceptVector)
}
InterceptNamesWIP<-rep(c("Ts"),times=length(ThresholdVector))
InterceptNames<-paste(InterceptNamesWIP,ThresholdVector,sep="")
colnames(InterceptMatrix)<-c(names(incomingfile)[namecolumn],InterceptNames)
#export the individual head Ts values
write.csv(InterceptMatrix,file=paste(Filename,"Intercepts Per Sample.csv"),row.names=FALSE)
###########################
#################Step 3
for(i in c(2:ncol(InterceptMatrix))){InterceptMatrix[,i]<-as.numeric(sub("Inf",NA,InterceptMatrix[,i]))}
#Average Calculations for each unique category in the name column
traitcols<-c(2:ncol(InterceptMatrix))
LoopColNames<-colnames(InterceptMatrix)[traitcols]
genotypelist<-unique(InterceptMatrix[,1])
compilation<-c(names(incomingfile)[namecolumn],LoopColNames)
i<-1
for(i in c(1:length(genotypelist))){
  currentdata<-InterceptMatrix[InterceptMatrix[,1]==genotypelist[i],]
  staple<-vector()
  for(j in c(1:length(traitcols))){staple<-c(staple,mean(as.numeric(currentdata[,traitcols[j]]),na.rm=TRUE))}
  compilation<-rbind(compilation,c(genotypelist[i],staple))
}
#now we format the final output and round the output numbers to 3 decimal places
colnames(compilation)<-compilation[1,]
compilation<-compilation[-1,]
for(i in c(1:nrow(compilation))){
  for(j in c(2:ncol(compilation))){
    compilation[i,j]<-round(as.numeric(compilation[i,j]),3)
  }
}
#finally, export the end result, the Ts means of each genotype/entry/etc. in the name column
write.csv(compilation,file=paste(Filename,"Ts Means.csv",sep=""),row.names=FALSE,col.names=FALSE)
###########################
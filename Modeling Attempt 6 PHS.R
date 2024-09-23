

setwd("F:/Camille GWAS Project/")
###pull in observations
PlantObservations<-read.csv("Gold Standard 03.07.23 Raw TCAP Spike Wetting Test scores fixed.csv")
##
##need to change two things with this formatting, should pull in columns beforehand and rework this line
PlantObservationsSubset<-as.data.frame(cbind(c(1:length(PlantObservations[53,8:12])),c(3,4,5,6,7),t(PlantObservations[53,8:12])))
colnames(PlantObservationsSubset)<-c("Datarow","Day","Score")

####This is our model ... in progress
attach(PlantObservationsSubset)
#PHSModel1<-function(Day, alpha, displacement) ((5.5)+(9/pi)*atan(alpha*(Day+displacement)))
#model1 <- nls(Score~PHSModel1(alpha,displacement), data=PlantObservationsSubset, start=list(alpha=4,displacement=-6))

modelworking<-nls(Score~((5.5)+(9/pi)*atan(alpha*(Day+displacement))),data=PlantObservationsSubset, start=list(alpha=4,displacement=-6))
nls(Score~((5.5)+(9/pi)*atan(alpha*(Day+displacement))),data=PlantObservationsSubset, start=list(alpha=4,displacement=-6))
plot(Score~Day,xlim=c(3,7),ylim=c(1,10))
lines(Day,predict(modelworking))
par(new=T)
summary(modelworking)
modelworking
######
Alpha<-as.numeric(modelworking$m$getAllPars()[1])
Disp<-as.numeric(modelworking$m$getAllPars()[2])
Arctanfit<-as.numeric(modelworking$m$getAllPars()[2])
as.numeric(modelworking$m$deviance())
residualvector<-as.vector(modelworking$m$resid())
residualvector
modelworking$m
tempcurve<-function(Day)((5.5)+(9/pi)*atan(Alpha*(Day+Disp)))
curve(tempcurve,3,7,ylim=c(1,10))
##this if for manually testing our assumptions with the arctan model
((5.5)+((4.5/(pi/2)))*atan(4*(3-6)))
####

##This module is for truncating the data to the part that is increasing, removing additional 1's or 10's at the ends
minrow<-max(PlantObservationsSubset[PlantObservationsSubset$Score==1,1])
maxrow<-min(PlantObservationsSubset[PlantObservationsSubset$Score==10,1])
if(minrow==-Inf){minrow<-1}
if(maxrow==Inf){maxrow<-as.numeric(nrow(PlantObservationsSubset))}
PlantObservationsSubsetTruncate<-PlantObservationsSubset[minrow:maxrow,]
detach(PlantObservationsSubset)
detach(PlantObservationsSubsetTruncate)
attach(PlantObservationsSubsetTruncate)
######################################

###module for testing linear model
linearmodel1<-nls(Score~((alpha*Day)+displacement),data=PlantObservationsSubsetTruncate, start=list(alpha=2,displacement=-6))
plot(Score~Day)
lines(Day,predict(linearmodel1))
####


min<-min(PlantObservationsSubset[PlantObservationsSubset$Score==10,1])
traitvector<-vector()
a<-1
#pull in observations 1 plant at a time
for(a in c(1:nrow(PlantObservations)))
{


  
traitvector<-c(traitvector,X)


ThresholdVector<-c()
b<-1
#for(b in c(1:length(ThresholdVector))){
}
PlantObservations<-cbind(PlantObservations,traitvector)



#}


##Formula to model by: attempt 1
#((5.5)+4.5*atan(alpha(days+displacement)))

#function where we vary 2 variables, right? alpha and displacement? 
#Alpha is the slope, displacement is the shifting of x-axis to line up 
#with rise by days? Maybe? We'll see.

#PHSModel1<-function(days, alpha, displacement) ((5.5)+4.5*atan(alpha*days+(displacement)))
#model1 <- nls(fluorI ~ PHSModel1(days,alpha,displacement), data=PlantObservations, start=list(myA=10,myT=5))

#code reference = "https://martinlab.chem.umass.edu/r-fitting-data/"
#eDecay <- function(t, ampl, tau) (ampl*exp(-t/tau))
#nls(fluorI ~ eDecay(t,myA,myT), data=ExpData, start=list(myA=10,myT=5))
#model1 <- nls(fluorI ~ eDecay(t,myA,myT), data=ExpData, start=list(myA=10,myT=5))

#plot(ExpData$t, residuals(model1), main="Residuals - Single Exponential fit", xlab="time (ns)", ylab="Residuals (predicted - observed)")
#plot(ExpData$t,ExpData$fluorI,xlab="time (ns)", ylab="fluorescence",main="Fluorescence Decay Assay - Double Exponential Fit")
#lines(ExpData$t,predict(model2))



####compare two data types for each observation, save all data: formula constants, fit statistic, residual string

####approx or approxfun
###equivalent


thresholds<-c(1.1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,9.9)

####okay, intercept attempt 1:

#for arctan function
arctanresiduals<-matrix(nrow=5,ncol=2)
###oh jeez, I'm going to rbind this in a loop? 5200 times? well.... it's going to start slowing down at the end.... but it's not that bad, right?
#it will be fine, yeah, totally fine..... <_< no problems here... >_> look.... it's a nicer output.... it's not going to take any longer than converting the wide data to long data would.... well... let's not do this for 100,000 samples.


####example code that I'm borrowing from https://stackoverflow.com/questions/53360961/how-to-find-x-intercept-for-intersection-of-two-curves-in-r
# Create an dataframe:
x <- seq(0, 10, 0.01)
y1 <- sin(x)
y2 <- cos(x)
df <- data.frame(x = x, y1 = y1, y2 = y2)
# Let's plot the curves
# I like to use the colorblind-friendly palette from
# Wong, Bang. 2011. "Points of view: Color blindness." Nature Methods 8:441.
wong_palette <- c("#e69f00", "#56b4e9")
# Plot the first curve
plot(x, y1, type = "l", xlab = "V1", ylab = "V2", col = wong_palette[1])
# Add the second
lines(x, y2, col = wong_palette[2])
# What are the intersections?
equivalent <- function(x, y, tol = 0.005) abs(x - y) < tol
xmin <- 3
xmax <- 8
intersection_indices <- which(equivalent(y1, y2) & x >= xmin & x <= xmax)
x[intersection_indices]
#> [1] 3.93 7.07
points(x[intersection_indices], y1[intersection_indices])

#########modifying
###okay, this works well enough for single points (though there are mild concerns about the tol value in the equivalent function being too stringent in the middle and too biased at the ends), 
#need to streamline and scale up
x <- seq(-10, 20, 0.01)
alpha<-Alpha
displacement<-Disp
y1 <- ((5.5)+(9/pi)*atan(alpha*(x+displacement)))
y2 <- 6.3
df <- data.frame(x = x, y1 = y1, y2 = y2)
# Let's plot the curves
# I like to use the colorblind-friendly palette from
# Wong, Bang. 2011. "Points of view: Color blindness." Nature Methods 8:441.
wong_palette <- c("#e69f00", "#56b4e9")
# Plot the first curve
plot(x, y1, type = "l", xlab = "V1", ylab = "V2", col = wong_palette[1])
# Add the second
#lines(x, y2, col = wong_palette[2])
abline(h=y2)
# What are the intersections?
#swc adding these contingencies to improve stability and accuracy<
if(y2<=9|y2>=2){z<-.1}
if(y2>9|y2<2){z<-.02}
###>
equivalent <- function(x, y, tol = z) abs(x - y) < tol
intersection_indices <- equivalent(y1, y2)
x[intersection_indices]
intersection_index<-mean(x[intersection_indices])
points(intersection_index, ((5.5)+(9/pi)*atan(alpha*(intersection_index+displacement))))
intersection_index

###Look, I know I probably should have just done the algebra/trigonometry/calculus to solve for X at a given Y with that arctan function instead of using this approximation, but this is (ironically) easier for me at the moment
##I really could look this up... but I'm not in the mood
#the approximations are useful in being general purpose and not requiring annoying algebra
#that sounds like a cop-out... ug, I sure know how to make myself feel guilty even while I get answers that are honestly close enough, and the other answers ought not have significant digits better than I'm getting anyway....
#whatever, I'm not doing algebra on this tonight. You humans can't make me!
#frak... I'm going to do algebra aren't I? yeah... here goes attempt 1:
#step 0 list formula;
# score=((5.5)+(9/pi)*atan(alpha*(Day+displacement)))
#step 1 subtract 5.5 from both sides
#score - 5.5 = (9/pi)*atan(alpha*(day+displacement))
#step 2 divide both sides by 9/pi
#(pi/9)*(score - 5.5) = atan(alpha*(day+displacement))
#step 3 take the tangent of both sides? yes, I think that's correct
#tan((pi/9)*(score - 5.5)) = alpha*(day+displacement)
#step 4 divide both sides by alpha:
#((tan((pi/9)*(score - 5.5))))/alpha) = day + displacement
#step 5 subtract displacement from both sides:
#((tan((pi/9)*(score - 5.5))))/alpha)-displacement = day
####ummm, did that do it? I think so?
#let's test the results
testingday<-(((tan((pi/9)*(y2-5.5))))/alpha)-displacement
testingday
round(testingday,digits=3)
###huzzah!!! I am math genius! well, no, but I am at least moderately competent.
#Maybe that perfect score on my math SAT wasn't a fluke after all
###well, now we don't have to do that ugly approximation function.
#Algebra-man, what a great super-hero
#the linear method is even easier to work with.

####okay, here goes my attempt at coding the script from the pieces
###what do we need?
#1: we need the right loop structure: 
#2: We need the right inputs
#3: we need to save the right outputs
#4: we need to have our organization correct
####to accomplish these things, I propose doing this in modules:
##Module 1: Here we run all arctans and save all arctans

#############Here we go, module 1 coded
setwd("F:/Camille GWAS Project/")
PlantObservations<-read.csv("Gold Standard 03.07.23 Raw TCAP Spike Wetting Test scores fixed.csv")
##
alphavector<-vector()
dispvector<-vector()
fitvector<-vector()
residualmatrix<-matrix(nrow=0,ncol=5)
colnames(residualmatrix)<-c("Day3resid","Day4resid","Day5resid","Day6resid","Day7resid")
#Core loop
i<-1
RowsToOmit<-c(159,517,1393,2779,2981,3662,4198,4311,4614,4908,4995,5256,5282)
for(i in c(1:nrow(PlantObservations)))
{
PlantObservationsSubset<-as.data.frame(cbind(c(1:length(PlantObservations[i,8:12])),c(3,4,5,6,7),t(PlantObservations[i,8:12])))
colnames(PlantObservationsSubset)<-c("Datarow","Day","Score")
##encountered row 159 with all 2's.... and I'm making an overriding list of omitted rows
if(i %in% RowsToOmit){
Alpha<-NA
Disp<-NA
Arctanfit<-NA
residualvector<-c(NA,NA,NA,NA,NA)
}
if(!(i %in% RowsToOmit)){
##Conditional statement to test if anything rises above 1
if(sum(PlantObservationsSubset$Score)<=5){
Alpha<-"All ones"
Disp<-"All ones"
Arctanfit<-"All ones"
residualvector<-c("All ones","All ones","All ones","All ones","All ones")
}
if(sum(PlantObservationsSubset$Score)>5 & sum(PlantObservationsSubset[1:4,3])<=4){
  Alpha<-"lastday"
  Disp<-"lastday"
  Arctanfit<-"lastday"
  residualvector<-c("lastday","lastday","lastday","lastday","lastday")
}
if(sum(PlantObservationsSubset$Score)>5 & sum(PlantObservationsSubset[1:4,3])>4){
attach(PlantObservationsSubset)
modelworking<-nls(Score~((5.5)+(9/pi)*atan(alpha*(Day+displacement))),data=PlantObservationsSubset, start=list(alpha=4,displacement=-6))
Alpha<-as.numeric(modelworking$m$getAllPars()[1])
Disp<-as.numeric(modelworking$m$getAllPars()[2])
Arctanfit<-as.numeric(modelworking$m$deviance())
residualvector<-as.vector(modelworking$m$resid())
detach(PlantObservationsSubset)
}}
#saving to vectors and matrices
alphavector<-c(alphavector,Alpha)
dispvector<-c(dispvector,Disp)
fitvector<-c(fitvector,Arctanfit)
residualmatrix<-rbind(residualmatrix,residualvector)
}
###loop ended, now staple together
ArctanResults<-cbind(PlantObservations,fitvector,residualmatrix,alphavector,dispvector)
write.csv(ArctanResults,file="ArctanModelResults_TCAP_Gold_Fixed_march23.csv")

##
i
alphavector
i
#notes on omitted rows
#159 was all 2's
###nothing apparently wrong with row 517, but took too many iterations
#1393 is all 2's
#2779 seems fine, but took too many iterations
#same for 2981
#same for 3662, I understand, not enough informative points
#4198 was ironically too susceptible for this model with our range
#4311 seems fine, but took too many iterations
#4614 seems fine, but took too many iterations
#4908 seems fine, but took too many iterations
#4995
#5256, all 2's
#5282, not sure, step factor

###yay, it worked
#now to make plots

##########Module 2
PlottingResiduals<-read.csv("ArctanModelResults_residuals_for_plotting_fixed.csv")
boxplot(PlottingResiduals$ResidualOfScore~PlottingResiduals$Day,main="Residuals Of Arctangent model over TCAP",ylab="Residual of Score",xlab="Day",col="green")

PlottingCurves<-read.csv("ArctanModelResults_alpha_disp.csv")
i<-1
for(i in c(1:nrow(PlottingCurves))){
Alpha<-PlottingCurves[i,1]
Disp<-PlottingCurves[i,2]
tempcurve<-function(Day)((5.5)+(9/pi)*atan(Alpha*(Day+Disp)))
curve(tempcurve,3,7,ylim=c(1,10))
par(new=T)
}
#it's a messy plot, but I like it.

####################################################################################
#####################################################################
##########Module 3
############Time to do stuff for the linear version of all this
#############Here we go, module 1 coded
setwd("D:/Camille GWAS Project/")
PlantObservations<-read.csv("Gold Standard 03.07.23 Raw TCAP Spike Wetting Test scores fixed.csv")
##
slopevector<-vector()
offsetvector<-vector()
fitvector<-vector()
residualmatrix<-matrix(nrow=0,ncol=5)
colnames(residualmatrix)<-c("Day3resid","Day4resid","Day5resid","Day6resid","Day7resid")
#Core loop
i<-1
RowsToOmit<-c()
for(i in c(1:nrow(PlantObservations)))
{
  PlantObservationsSubset<-as.data.frame(cbind(c(1:length(PlantObservations[i,8:12])),c(3,4,5,6,7),t(PlantObservations[i,8:12])))
  colnames(PlantObservationsSubset)<-c("Datarow","Day","Score")
  NAmaker<-vector()
  ##encountered row 159 with all 2's.... and I'm making an overriding list of omitted rows
  if(i %in% RowsToOmit){
    SlopeValue<-NA
    OffsetValue<-NA
    residualvector<-c(NA,NA,NA,NA,NA)
  }
  if(!(i %in% RowsToOmit)){
    #Here inserting logic to translate preliminary 1's and subsequent 10's into NA's for the linear fit
    minrow<-max(PlantObservationsSubset[PlantObservationsSubset$Score==1,1])
    maxrow<-min(PlantObservationsSubset[PlantObservationsSubset$Score==10,1])
    if(minrow==-Inf){minrow<-1}
    if(minrow==as.numeric(nrow(PlantObservationsSubset))){minrow<-1}
    if(maxrow==Inf){maxrow<-as.numeric(nrow(PlantObservationsSubset))}
    if(maxrow==1){maxrow<-1}
    if(minrow>1){NAmaker<-c(1:(minrow-1))}
    if(maxrow<as.numeric(nrow(PlantObservationsSubset))){NAmaker<-c(NAmaker,maxrow+1:as.numeric(nrow(PlantObservationsSubset)))}
    #replace things
    PlantObservationsSubset[NAmaker,3]=NA
    ####huzzah, I feel clever. The several rows above are a little verbose, but 'should' take care of all contingencies.
    ##removed the other logic that was necessary for the arctan function. A line only requires 2 points afterall.
      attach(PlantObservationsSubset)
      linearmodel1<-lm(Score~Day,data=PlantObservationsSubset)
      SlopeValue<-as.numeric(round(linearmodel1$coefficients[2],digits=4))
      OffsetValue<-as.numeric(round(linearmodel1$coefficients[1],digits=4))
      residualvector<-as.vector(round(linearmodel1$residuals,digits=4))
      detach(PlantObservationsSubset)
    }
  ##saving to vectors and matrices
  #removed fitvector because it was being a pain, and was giving warning messages about unreliable values with perfect fits
  #I don't think I need it...I'll add it back in if I'm wrong
  
  #I got a bit of help from the internet, if I do decide to add it back in, this code will help
  #summary(linearmodel1)$adj.r.squared
  slopevector<-c(slopevector,SlopeValue)
  offsetvector<-c(offsetvector,OffsetValue)
  residualmatrix<-rbind(residualmatrix,residualvector)
}
###loop ended, now staple together
LinearResults<-cbind(PlantObservations,residualmatrix,slopevector,offsetvector)
write.csv(LinearResults,file="LinearModelResults_TCAP_Gold_Fixed_march2023.csv")
i

#####################
###########Errr, got some warnings, but I think it worked? Wait.... maybe not...
#checking that my slope and intercept variables are correct
#They were wrong, now they're fixed
linearmodel1$model
summary(linearmodel1)
linearmodel1$coefficients
#let's take a look at the individual points and lines, to see how well it's working.

pdf(file="Examining fits of linear models for TCAP.pdf")
par(mfrow=c(5,5))
for(i in c(1:nrow(LinearResults)))
{
  plot(c(3:7),LinearResults[i,8:12],xlim=c(3,7),ylim=c(1,10))
  abline(a=LinearResults[i,28],b=LinearResults[i,27])
}
dev.off()

###calculated slope and intercept properly, yay
#however, I think the residuals are not correct, there should be lots of NA's, and there aren't.


#################Alright, here's where things stand.
#We have successful mathematical models for all the samples using both the arctangent and linear regression models
#the residuals we calculated are okay for the arctangent models, but don't work for the linear regression model
#However, the method I used previously for plotting residuals by day with boxplots is, upon reflection, actually hot garbage
#You heard me Module 2, your boxplots are not the right way to look at this.
#I didn't account for the fact that all the curves were at different stages. That's not good.
#so, how do we fix it? By consulting the Sacred Chickens? No!
#We will fix this with beautiful, brilliant, abstract glory
#we will do algebra on the vectors, and solve for the X intercept at Y=5.5 (the vertical midpoint)
#we will then displace all our equations by this value.
#At this point, the logic splits, where we want 1 set that is scaled for slope, and one that isn't
#we backtrack a little bit, we take the residuals on each day (these need to be properly recalculated for the linear models)
#then we offset the X value associated with those residuals by the X value we obtained for the vertical midpoint.
#then, we plot them all... all of them... O_O all of them.

#we do want to look at genotype effects within this. Maybe plot out all the different curves for each genotype, corrected about the midpoint, on one plot, with used n=x displayed in the corner
#

#okay, time to just start trying this, there are a few aspects I still need to tinker with for the scaled part, and where I draw the range cutoffs.
#I should consider using the time residuals instead of the score residuals for assessing the linear model
setwd("D:/Camille GWAS Project/")
Linear<-read.csv("LinearModelResults_TCAP_Gold_Fixed_march2023.csv")
Arctan<-read.csv("ArctanModelResults_TCAP_Gold_Fixed_march23.csv")
ArctanSmaller<-Arctan[,c(17,29,30)]
LinearSmaller<-Linear[,c(17,28,29)]

LinearSlopes<-LinearSmaller[LinearSmaller$slopevector!=0,]
MidpointInterceptLinearSlopes<-((5.5 - LinearSlopes$offsetvector)/LinearSlopes$slopevector)
LinearPartlySolved<-cbind(LinearSlopes,MidpointInterceptLinearSlopes)
#y=mx+b
#y-b=mx
#(y-b)/m=x
#y=5.5
#x is in days

##test this
randoms<-sample(c(1:4000),16)

i<-1
par(mfrow=c(1,1))
for(i in c(1:length(randoms)))
{
plot(LinearPartlySolved[randoms[i],4],5.5,ylim=c(1,10),xlim=c(1,15),col="red",xlab="day",ylab="score",main="Testing Intercepts")
abline(h=5.5,col="blue")
abline(a=c(LinearPartlySolved[randoms[i],3]),c(b=LinearPartlySolved[randoms[i],2]))
par(new=TRUE)
#print(paste("slope",i,LinearPartlySolved[randoms[i],2]))
#print(paste("offset",i,LinearPartlySolved[randoms[i],3]))
#print(paste("intersect",i,randoms[i]))
}
par(new=FALSE)
####well, yay, I did the test, and isn't that a pretty graph.
#note to future self: It will be different every time

###Okay, for the linear stuff I have the time offsets that I need to adjust for
#So, Y=mx+b, but instead, Y=m(x-d)+b
#This is going to get weird....
#these are the variables I need to use
Slope<-LinearPartlySolved$slopevector
Offset<-LinearPartlySolved$offsetvector
Alignment<-LinearPartlySolved$MidpointInterceptLinearSlopes
#I want to plot out all the lines
#I also want to calculate all the residuals with respect to time, but not quite yet
#let's start by plotting everything

#Here I render the plot, spoiler alert, it works.
png(filename = "All TCAP PHS Linear Models Aligned.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-5,5),main="All TCAP PHS Linear Models Aligned",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
for(i in c(1:nrow(LinearPartlySolved))){
tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
curve(tempcurve,-5,5,ylim=c(1,10),xlab="",ylab="",main="")
par(new=TRUE)
}
par(new=FALSE)
dev.off()
################I also want a version without the alignment
png(filename = "All TCAP PHS Linear Models Before Alignment.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models",xlab="Time (days)",ylab="Score")
par(new=TRUE)
for(i in c(1:nrow(LinearPartlySolved))){
  tempcurve<-function(Day)((Slope[i]*(Day))+Offset[i])
  curve(tempcurve,1,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
dev.off()
###### Okay, those plots are rendered and look good
#I could do all of this for the arctangent models immediately, but... I want to deal with the residuals first
#Now I need to do the residual calculations where the residuals are in time.
###so, in order to do this I kind of need to do the same thing I did before, only backwards
#I need a loop, where I manually calculate the temporal residuals on each observation that's non-zero
#and the filtering I did so far prevents me from doing this on the PartlySolved dataset
#residual calculations starting now... my brain is exhausted already, and it's Friday night, but who cares, R is my social life, and math is all the company I need. Eventually, I'll probably die hungry and alone, but I'll be able to look back and say "you know, those R plots I made, those were really cool"

#what are the conditions under which I don't want residuals calculated?
#1) For values of 1 and 10, I want the same filtering I used earlier
#2) For slopes of 0, I don't want to calculate temporal residuals.
#Those are the main two criteria
#however, I also want to sample the dataset for those individuals that have at least 2 non-1 and non-10 points.

#Original Points
png(filename = "All TCAP PHS Just Datapoints.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(3,7),main="All TCAP PHS Pure Datapoints",xlab="Time (days)",ylab="Score")
par(new=TRUE)
for(i in c(1:nrow(Linear))){
  plot(c(3,4,5,6,7),Linear[i,9:13],ylim=c(1,10),xlim=c(3,7),main="",xlab="",ylab="",type="o")
  par(new=TRUE)
}
par(new=FALSE)
dev.off()

###converted our original datasets to something nicer to deal with
Linear2<-Linear[,-c(1,18:27)]
Arctan2<-Arctan[,-c(18:23)]
###Now, Let's focus on the linear one, we need to get the intercepts at 0 and 10 for those that have a non-zero slope
#let's get the code working
PosSlopeLinear2<-Linear2[Linear2$slopevector!=0,]
#while we're at it, let's also make the version with at least 2 non-tail datapoints
temp1<-PosSlopeLinear2[,8:12]!=(1|10)
temp1[temp1==TRUE]<-1
Countingvector<-rowSums(temp1)
PosSlopeMultiObsLinear2<-PosSlopeLinear2[Countingvector>=2,]


EndpointInterceptPosLinearSlopes2<-((10 - PosSlopeLinear2$offsetvector)/PosSlopeLinear2$slopevector)
MidpointInterceptPosLinearSlopes2<-((5.5 - PosSlopeLinear2$offsetvector)/PosSlopeLinear2$slopevector)
StartpointInterceptPosLinearSlopes2<-((1 - PosSlopeLinear2$offsetvector)/PosSlopeLinear2$slopevector)
PosSlopeLinear2prime<-cbind(PosSlopeLinear2,StartpointInterceptPosLinearSlopes2,MidpointInterceptPosLinearSlopes2,EndpointInterceptPosLinearSlopes2)
#same thing but for the multiple (3 to 5) meaningful observation subset
EndpointInterceptPosLinearSlopesMultiObs2<-((10 - PosSlopeMultiObsLinear2$offsetvector)/PosSlopeMultiObsLinear2$slopevector)
MidpointInterceptPosLinearSlopesMultiObs2<-((5.5 - PosSlopeMultiObsLinear2$offsetvector)/PosSlopeMultiObsLinear2$slopevector)
StartpointInterceptPosLinearSlopesMultiObs2<-((1 - PosSlopeMultiObsLinear2$offsetvector)/PosSlopeMultiObsLinear2$slopevector)
PosSlopeMultiObsLinear2prime<-cbind(PosSlopeMultiObsLinear2,StartpointInterceptPosLinearSlopesMultiObs2,MidpointInterceptPosLinearSlopesMultiObs2,EndpointInterceptPosLinearSlopesMultiObs2)
#############Now to calculate residuals in 2 directions, before the correction
#First, need to reconvert all those non-informative points to NA's
#oh jeez, I can sort of see how to do all this in the abstract in just one step, but I can't quite hold it all in my head, and am afraid I'd screw it up
#So, I'll do things step by step
#converting non informative points to NA's, borrowing earlier code and applying it differently

#I'm just going to use this same loop, and tweak the input and output.
#I'm actually going to save the output as a list, instead of as a larger matrix.
#The list will give me more versatility, with each observation having its own matrix.j
#hard to export lists, but internally, they're super nifty.
###Here, define what you're inputting
PlantObservations<-PosSlopeLinear2prime
#change the observation columns for future use, and bear in mind that this is set up for days 3:7
ObservationColumns<-c(8:12)
i<-1
NewList<-list()
for(i in c(1:nrow(PlantObservations)))
{
PlantObservationsSubset<-as.data.frame(cbind(c(1:length(PlantObservations[i,ObservationColumns])),c(3,4,5,6,7),t(PlantObservations[i,ObservationColumns])))
colnames(PlantObservationsSubset)<-c("Datarow","Day","Score")
NAmaker<-vector()
#identifying the bounds
minrow<-max(PlantObservationsSubset[PlantObservationsSubset$Score==1,1])
maxrow<-min(PlantObservationsSubset[PlantObservationsSubset$Score==10,1])
if(minrow==-Inf){minrow<-1}
if(minrow==as.numeric(nrow(PlantObservationsSubset))){minrow<-1}
if(maxrow==Inf){maxrow<-as.numeric(nrow(PlantObservationsSubset))}
if(maxrow==1){maxrow<-1}
if(minrow>1){NAmaker<-c(1:(minrow-1))}
if(maxrow<as.numeric(nrow(PlantObservationsSubset))){NAmaker<-c(NAmaker,maxrow+1:as.numeric(nrow(PlantObservationsSubset)))}
#replacing things
PlantObservationsSubset[NAmaker,3]=NA
#save
NewList[[i]]<-PlantObservationsSubset
}
PosSlopeLinear2primeList<-NewList
#brilliant, it works. 
#It gives warnings about infinity, but I have touched the whirling infinite, and such warnings are not needed for one with conditionals such as I
#going to quickly replicate this for the multiObs set

PlantObservations<-PosSlopeMultiObsLinear2prime
#change the observation columns for future use, and bear in mind that this is set up for days 3:7
ObservationColumns<-c(8:12)
i<-1
NewList<-list()
for(i in c(1:nrow(PlantObservations)))
{
  PlantObservationsSubset<-as.data.frame(cbind(c(1:length(PlantObservations[i,ObservationColumns])),c(3,4,5,6,7),t(PlantObservations[i,ObservationColumns])))
  colnames(PlantObservationsSubset)<-c("Datarow","Day","Score")
  NAmaker<-vector()
  #identifying the bounds
  minrow<-max(PlantObservationsSubset[PlantObservationsSubset$Score==1,1])
  maxrow<-min(PlantObservationsSubset[PlantObservationsSubset$Score==10,1])
  if(minrow==-Inf){minrow<-1}
  if(minrow==as.numeric(nrow(PlantObservationsSubset))){minrow<-1}
  if(maxrow==Inf){maxrow<-as.numeric(nrow(PlantObservationsSubset))}
  if(maxrow==1){maxrow<-1}
  if(minrow>1){NAmaker<-c(1:(minrow-1))}
  if(maxrow<as.numeric(nrow(PlantObservationsSubset))){NAmaker<-c(NAmaker,maxrow+1:as.numeric(nrow(PlantObservationsSubset)))}
  #replacing things
  PlantObservationsSubset[NAmaker,3]=NA
  #save
  NewList[[i]]<-PlantObservationsSubset
}
PosSlopeMultiObsLinear2primeList<-NewList


##now we need to input Y (score) values into the formula, and extract X values for predictions.
#the specific scores we need are the scores within each list there. Perhaps having the NA's will even help synchronize things.
#let's try it. So, we open up each list entry and we append a few columns to it.
#additional columns: Temporal Residuals, Score Residuals, Temporal Range Corrected, Temporal Residuals Range Corrected

Sourcedata<-PosSlopeLinear2prime
SourceList<-PosSlopeLinear2primeList

Sourcedata<-PosSlopeMultiObsLinear2prime
SourceList<-PosSlopeMultiObsLinear2primeList
i<-1
for(i in c(1:nrow(Sourcedata))){
#clear vectors for safety
TemporalResiduals<-vector()
ScoreResiduals<-vector()
TemporalRangeCorrected<-vector()
TemporalResidualsRangeCorrected<-vector()
###Variables used
SourceList[[i]]$Day
SourceList[[i]]$Score
Slope<-Sourcedata[i,17]
Intercept<-Sourcedata[i,18]
Startpoint<-Sourcedata[i,19]
Midpoint<-Sourcedata[i,20]
Endpoint<-Sourcedata[i,21]
#The Slope and intercept define the regression line, and the Start, Mid and Endpoints define the X values (days)
#1 , 5.5 and 10 are the Y values (score) for the Start, Mid and Endpoints
#Alright, now that we have our variables more clearly defined, let's calculate our vectors
#For temporal residuals, we calculate for the score values we have, SourceList[[i]]$Score, the X values according to the model. Then we subtract the models X values from ours.
#The temporal values are calculated with:
#((SourceList[[i]]$Score - Intercept)/Slope)
TemporalResiduals<-(SourceList[[i]]$Day -((SourceList[[i]]$Score - Intercept)/Slope))
#the above should work
#to calculate score residuals, we solve for what the model thinks the score should be on that day, and then we subtract the predicted score from the observed score
ScoreResiduals<-(SourceList[[i]]$Score - ((Slope*SourceList[[i]]$Day)+Intercept))
#should work
#For correcting the data, we're going to transpose it onto a relative 0 to 100 scale
#So, Endpoint - Startpoint = range
#subtract startpoint from all day calculations, so the start is always 0
#then divide the days by the range, then multiply the rangecorrected days by 100
#we also have to do this for all the residuals, so get it smooth....
TemporalRangeCorrected<- 100*((SourceList[[i]]$Day-Startpoint)/(Endpoint-Startpoint))
#the above should work, wait, minor problem, with NA's not crossing over, fixing it with the following line:
TemporalRangeCorrected<-(SourceList[[i]]$Score/SourceList[[i]]$Score)*TemporalRangeCorrected
###errr, I think, because the residuals are already based on the curve, that I just need to divide by the range
TemporalResidualsRangeCorrected<-100*((TemporalResiduals)/(Endpoint-Startpoint))
#all of the above appears to be good, but I will need further testing within plots
#however, to try to confirm that what I've done is actually correct, I also want to do the same transformation on the fit lines
#Need a new slope and a new intercept: if everything is done correctly, the slope should be .1 and the intercept should be 0 in all cases
#so, let's apply out math to the formula... somehow ... brain... why aren't you working.... aaaaaggghhhhh
#each observation has an intercept range for its line, we calculated that
#All we're doing is applying a horizontal offset, and multiplying everything by 100/range ... right?
#if that's true, then, we ... apply those terms to the Y side of the equation, and then shift them over to the X? wait... ug,
#okay, here's the final observation, let's use that
#Slope
#Intercept
#we're changing the x variable here, so the x term needs to be expanded
#plot(NA,NA,xlim=c(-50,150),ylim=c(1,100))
#lines((Slope*((100/(Endpoint-Startpoint))*(SourceList[[i]]$Day-Startpoint)))+Intercept)
#Nevermind, this.... needs me to rethink this. I'm doing it wrong.
#The logic that I used to correct the points should all work.

#well, let's go for it
SourceList[[i]]<-cbind(SourceList[[i]],TemporalResiduals,ScoreResiduals,TemporalRangeCorrected,TemporalResidualsRangeCorrected)
}
SourceList[[i]]

#Nifty trick: did you know you can save lists? The commands are as follows
write.csv(Sourcedata,file="PosSlopeLinear2PrimeMatrix.csv")
saveRDS(SourceList,file="PosSlopeLinear2PrimeListResults.RData")
write.csv(Sourcedata,file="PosSlopeMultiObsLinear2PrimeMatrix.csv")
saveRDS(SourceList,file="PosSlopeMultiObsLinear2PrimeListResults.RData")

#okay, small issue, for most of the stuff it's working great, but for the temporal range corrected it's not using the NA's for the ommitted points
#so, I need to fix that
#if I so something like: TemporalRangeCorrected<-(SourceList[[i]]$Score/SourceList[[i]]$Score)*TemporalRangeCorrected
#Good, fixed.

#Now, I need to save and load the 4 files:

PosSlopeLinear2primeList<-readRDS("PosSlopeLinear2PrimeListResults.RData")
PosSlopeMultiObsLinear2primeList<-readRDS("PosSlopeMultiObsLinear2PrimeListResults.RData")
PosSlopeLinear2prime<-read.csv("PosSlopeLinear2PrimeMatrix.RData")
PosSlopeMultiObsLinear2prime<-read.csv("PosSlopeMultiObsLinear2PrimeMatrix.RData")

#okay, so here's our data
PosSlopeLinear2prime
PosSlopeLinear2primeList
PosSlopeMultiObsLinear2prime
PosSlopeMultiObsLinear2primeList


################Time to make plots of all this crazy math
#alternate these two inputs to the plot I want
InputFile<-PosSlopeLinear2prime
InputList<-PosSlopeLinear2primeList
png(filename = "All TCAP PHS Linear Models Without Any Correction Prime.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models Points No Correction Prime",xlab="Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
for(i in c(1:nrow(InputFile))){
  plot(InputList[[i]]$Day,InputList[[i]]$Score,ylim=c(1,10),xlim=c(1,10),col="red",xlab="",ylab="",main="")
  abline(a=Offset[i],b=Slope[i])
  par(new=TRUE)
}
par(new=FALSE)
dev.off()

#for this one, change the plot function #Takes about 45 seconds to run
png(filename = "All TCAP PHS Linear Models Midpoint Sync Prime.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-10,10),main="All TCAP PHS Linear Models Points Midpoint Sync Prime",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$Day-Alignment[i]),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-10,10),col="red",xlab="",ylab="",main="")
  par(new=TRUE)
  tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
  curve(tempcurve,-10,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
dev.off()


#for this one, change the plot function #Takes about 45 seconds to run
png(filename = "All TCAP PHS Linear Models Range Sync Prime.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-50,150),main="All TCAP PHS Linear Models Points Range Sync Prime",xlab="Percent Completion",ylab="Score")
abline(a=1,b=9/100)
par(new=TRUE)
plot(c(0,100),c(1,10),ylim=c(1,10),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(1,0,0,.1))
  par(new=TRUE)
}
par(new=FALSE)
dev.off()

#######################
###########Duplicating the above but changing the inputs
#####################
InputFile<-PosSlopeMultiObsLinear2prime
InputList<-PosSlopeMultiObsLinear2primeList
png(filename = "All TCAP PHS Linear Models MultiObs Without Any Correction Prime.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models MultiObs Points No Correction Prime",xlab="Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
for(i in c(1:nrow(InputFile))){
  plot(InputList[[i]]$Day,InputList[[i]]$Score,ylim=c(1,10),xlim=c(1,10),col="green4",xlab="",ylab="",main="")
  abline(a=Offset[i],b=Slope[i])
  par(new=TRUE)
}
par(new=FALSE)
dev.off()

#for this one, change the plot function #Takes about 45 seconds to run
png(filename = "All TCAP PHS Linear Models MultiObs Midpoint Sync Prime.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-10,10),main="All TCAP PHS Linear Models MultiObs Points Midpoint Sync Prime",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$Day-Alignment[i]),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-10,10),col="green4",xlab="",ylab="",main="")
  par(new=TRUE)
  tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
  curve(tempcurve,-10,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
dev.off()


#for this one, change the plot function #Takes about 45 seconds to run
png(filename = "All TCAP PHS Linear Models MultiObs Range Sync Prime.png",width = 480, height = 480)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Points Range Sync Prime",xlab="Percent Completion",ylab="Score")
abline(a=1,b=9/100)
par(new=TRUE)
plot(c(0,100),c(1,10),ylim=c(1,10),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(.1,1,.1,.1))
  par(new=TRUE)
}
par(new=FALSE)
dev.off()

#############################################################
##################################
####################################################################################
#Now, to take all that, and stick it in 1 PDF, 6 plots per page, 1st page is overall set

GrandInput1<-PosSlopeLinear2prime
GrandInputList1<-PosSlopeLinear2primeList
GrandInput2<-PosSlopeMultiObsLinear2prime
GrandInputList2<-PosSlopeMultiObsLinear2primeList

#pdf(file="Grand Design Linear Models TCAP PHS All and Individuals.pdf")
png(filename = "Grand Design Litear Models TCAP PHS All.png",width = 1600, height = 800)
par(mfrow=c(2,4))
############First, do everything
InputList<-GrandInputList1
InputFile<-GrandInput1
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models Points No Correction Prime",xlab="Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
for(i in c(1:nrow(InputFile))){
  plot(InputList[[i]]$Day,InputList[[i]]$Score,ylim=c(1,10),xlim=c(1,10),col="red",xlab="",ylab="",main="")
  abline(a=Offset[i],b=Slope[i])
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-10,10),main="All TCAP PHS Linear Models Points Midpoint Sync Prime",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$Day-Alignment[i]),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-10,10),col="red",xlab="",ylab="",main="")
  par(new=TRUE)
  tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
  curve(tempcurve,-10,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-50,150),main="All TCAP PHS Linear Models Points Range Sync Prime",xlab="Percent Completion",ylab="Score")
abline(a=1,b=9/100)
par(new=TRUE)
plot(c(0,100),c(1,10),ylim=c(1,10),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(1,0,0,.1))
  par(new=TRUE)
}
par(new=FALSE)
#4th plot, with residuals:
i<-1
plot(NA,NA,ylim=c(-3,3),xlim=c(-50,150),main="All TCAP PHS Linear Models Score Residuals",xlab="Percent Completion",ylab="Score Residual")
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-3,3),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$ScoreResiduals,ylim=c(-3,3),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(1,0,0,.1))
  par(new=TRUE)
}
par(new=FALSE)
#########Doing everything for 2nd set of plots with more filtered observations
InputFile<-GrandInput2
InputList<-GrandInputList2
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models MultiObs Points No Correction Prime",xlab="Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
for(i in c(1:nrow(InputFile))){
  plot(InputList[[i]]$Day,InputList[[i]]$Score,ylim=c(1,10),xlim=c(1,10),col="green4",xlab="",ylab="",main="")
  abline(a=Offset[i],b=Slope[i])
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-10,10),main="All TCAP PHS Linear Models MultiObs Points Midpoint Sync Prime",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopesMultiObs2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$Day-Alignment[i]),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-10,10),col="green4",xlab="",ylab="",main="")
  par(new=TRUE)
  tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
  curve(tempcurve,-10,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
##########
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Points Range Sync Prime",xlab="Percent Completion",ylab="Score")
abline(a=1,b=9/100)
par(new=TRUE)
plot(c(0,100),c(1,10),ylim=c(1,10),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(.1,1,.1,.1))
  par(new=TRUE)
}
par(new=FALSE)
#4th plot, with residuals:
i<-1
plot(NA,NA,ylim=c(-3,3),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Score Residuals",xlab="Percent Completion",ylab="Score Residual")
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-3,3),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$ScoreResiduals,ylim=c(-3,3),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(.1,1,.1,.1))
  par(new=TRUE)
}
par(new=FALSE)
###################################
######End of totality graphs
dev.off()


###########################Interdependent but important module for assessing residuals
####Wait, want cross comparison between score residuals and time residuals, I think apples to apples
###Also, let's make a version that has a much lower opacity, 1% opacity, instead of 10%
###grab multiobs, it's the more conservative (less centrally weighted) one for this
png(filename = "Grand Design Residual Comparisons MultiObs.png",width = 1200, height = 800)
InputFile<-GrandInput2
InputList<-GrandInputList2
par(mfrow=c(2,3))
#original
i<-1
plot(NA,NA,ylim=c(-3,3),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Score Residuals",xlab="% Completion",ylab="Score Residual",cex.axis=2,cex.lab=1.8)
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-3,3),xlim=c(-50,150),col="blue",xlab="",ylab="",main="",axes=FALSE)
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$ScoreResiduals,ylim=c(-3,3),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(0,0,0,.1),axes=FALSE)
  par(new=TRUE)
}
par(new=FALSE)
#temporal residuals
i<-1
plot(NA,NA,ylim=c(-3,3),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Time Residuals",xlab="% Completion",ylab="Relative Days Residual",cex.axis=2,cex.lab=1.8)
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-3,3),xlim=c(-50,150),col="blue",xlab="",ylab="",main="",axes=FALSE)
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$TemporalResiduals,ylim=c(-3,3),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(0,0,0,.1),axes=FALSE)
  par(new=TRUE)
}
par(new=FALSE)
#Temporal Residuals Range Corrected
i<-1
plot(NA,NA,ylim=c(-25,25),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs %Completion Time Residuals",xlab="% Completion",ylab="% Completion Residual",cex.axis=2,cex.lab=1.8)
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-25,25),xlim=c(-50,150),col="blue",xlab="",ylab="",main="",axes=FALSE)
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$TemporalResidualsRangeCorrected,ylim=c(-25,25),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(0,0,0,.1),axes=FALSE)
  par(new=TRUE)
}
par(new=FALSE)
####################Okay, do all of that again, but change the opacity, but change the 
i<-1
plot(NA,NA,ylim=c(-3,3),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Score Residuals",xlab="% Completion",ylab="Score Residual 1%-opacity",cex.axis=2,cex.lab=1.8)
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-3,3),xlim=c(-50,150),col="blue",xlab="",ylab="",main="",axes=FALSE)
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$ScoreResiduals,ylim=c(-3,3),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(0,0,0,.01),axes=FALSE)
  par(new=TRUE)
}
par(new=FALSE)
#temporal residuals
i<-1
plot(NA,NA,ylim=c(-3,3),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Time Residuals",xlab="% Completion",ylab="Relative Days Residual 1%-opacity",cex.axis=2,cex.lab=1.8)
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-3,3),xlim=c(-50,150),col="blue",xlab="",ylab="",main="",axes=FALSE)
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$TemporalResiduals,ylim=c(-3,3),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(0,0,0,.01),axes=FALSE)
  par(new=TRUE)
}
par(new=FALSE)
#Temporal Residuals Range Corrected
i<-1
plot(NA,NA,ylim=c(-25,25),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs %Completion Time Residuals",xlab="% Completion",ylab="% Completion Residual 1%-opacity",cex.axis=2,cex.lab=1.8)
abline(a=0,b=0)
par(new=TRUE)
plot(c(0,100),c(0,0),ylim=c(-25,25),xlim=c(-50,150),col="blue",xlab="",ylab="",main="",axes=FALSE)
par(new=TRUE)
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$TemporalResidualsRangeCorrected,ylim=c(-25,25),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(0,0,0,.01),axes=FALSE)
  par(new=TRUE)
}
par(new=FALSE)
dev.off()

##############################




#Alright, that should take care of page 1, now for all the other pages
Genotypes<-unique(GrandInput1$Genotype)
#cycling through all genotypes for individual pages
j<-3
for(j in c(1:length(Genotypes))){
samples<-c(GrandInput1$Genotype==Genotypes[j])
samples<-samples*c(1:length(samples))
samples<-samples[samples!=0]
InputFile<-GrandInput1[samples,]
InputList<-as.list(GrandInputList1[[samples]])

InputList<-GrandInputList1[[2]]

i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models Points No Correction Prime",xlab="Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
for(i in c(1:nrow(InputFile))){
  plot(InputList[[i]]$Day,InputList[[i]]$Score,ylim=c(1,10),xlim=c(1,10),col="red",xlab="",ylab="",main="")
  abline(a=Offset[i],b=Slope[i])
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-10,10),main="All TCAP PHS Linear Models Points Midpoint Sync Prime",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$Day-Alignment[i]),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-10,10),col="red",xlab="",ylab="",main="")
  par(new=TRUE)
  tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
  curve(tempcurve,-10,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-50,150),main="All TCAP PHS Linear Models Points Range Sync Prime",xlab="Percent Completion",ylab="Score")
abline(a=1,b=9/100)
par(new=TRUE)
plot(c(0,100),c(1,10),ylim=c(1,10),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(1,0,0,.1))
  par(new=TRUE)
}
par(new=FALSE)
###########Duplicating the above but changing the inputs
InputFile<-GrandInput2[GrandInput1$Genotype==Genotypes[j],]
InputList<-GrandInputList2[[GrandInput1$Genotype==Genotypes[j]]]
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(1,10),main="All TCAP PHS Linear Models MultiObs Points No Correction Prime",xlab="Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
for(i in c(1:nrow(InputFile))){
  plot(InputList[[i]]$Day,InputList[[i]]$Score,ylim=c(1,10),xlim=c(1,10),col="green4",xlab="",ylab="",main="")
  abline(a=Offset[i],b=Slope[i])
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-10,10),main="All TCAP PHS Linear Models MultiObs Points Midpoint Sync Prime",xlab="Relative Time (days)",ylab="Score")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$Day-Alignment[i]),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-10,10),col="green4",xlab="",ylab="",main="")
  par(new=TRUE)
  tempcurve<-function(Day)((Slope[i]*(Day+Alignment[i]))+Offset[i])
  curve(tempcurve,-10,10,ylim=c(1,10),xlab="",ylab="",main="")
  par(new=TRUE)
}
par(new=FALSE)
i<-1
plot(NA,NA,ylim=c(1,10),xlim=c(-50,150),main="All TCAP PHS Linear Models MultiObs Points Range Sync Prime",xlab="Percent Completion",ylab="Score")
abline(a=1,b=9/100)
par(new=TRUE)
plot(c(0,100),c(1,10),ylim=c(1,10),xlim=c(-50,150),col="blue",xlab="",ylab="",main="")
par(new=TRUE)
Slope<-InputFile$slopevector
Offset<-InputFile$offsetvector
Alignment<-InputFile$MidpointInterceptPosLinearSlopes2
for(i in c(1:nrow(InputFile))){
  plot((InputList[[i]]$TemporalRangeCorrected),InputList[[i]]$Score,ylim=c(1,10),xlim=c(-50,150),xlab="",ylab="",main="",pch=19,col=rgb(.1,1,.1,.1))
  par(new=TRUE)
}
par(new=FALSE)
}
dev.off()




## Author: Bruce Maxwell
## 04/05/2018
## Description: linear and non-linear regression analysis on yield and protein and 
## use of the models to do comparative net return calculations for different
## side-dress nitrogen rate application strategies on the field given the models.
## Farm: Oeste Field: Oeste25   Year harvested 2018   Crop: Winter wheat.

##*************************************************************************************
## Mapping previous same crop yields in field and show nitrogen rate experiment layout
##*************************************************************************************
## Author: Nick Silverman
## 05/11/2017
##************************* 

## Import libraries
library(rgdal)
library(raster)

## Set constants
prev.yr <- "2016"  #Don't know if this is correct
field.name <- "Oeste25"

## Bring in data
#dpath <- '../../../data/'  #Set the path to the working directory where all files reside
d <- read.csv("oeste_v1.csv",header=T)
colnames(d)
hist(d$fert_rate)  #Easy way to see if there are  fert rates that are way too high to be real
which(d$fert_rate>500) #Show list of rows from data file where Nitrogen fert rates exceed 500 kg/ha
d <- d[!(d$fert_rate > 500), ]          # get rid of outlier N rates
hist(d$fert_rate) #Frequency distribution of fert rates after unrealistically high rates removed

## *************  Map data ****************************
## Clean up data for yield and nitrogen maps    #
Ypts <- data.frame(x=d$x, y=d$y, YldPrev=d$prev_yl, Yld=d$yield, N.Rate=d$fert_rate)   ## YldPrev=d$yl_avg2016
abs.max <- max(cbind(Ypts$YldPrev, Ypts$Yld))
abs.min <- min(cbind(Ypts$YldPrev, Ypts$Yld))
cuts <- seq(0, round(abs.max, -1), 10)

## Convert yield data to raster
Ymat <- as.matrix(Ypts)
sp <- SpatialPoints(coords=Ymat[, c('x', 'y')])
e <- extent(Ypts[, c('x', 'y')])
rast <- raster(ext=e, resolution=10)
rastYldPrev <- rasterize(sp, rast, Ypts$YldPrev, fun=mean)
rastN <- rasterize(sp, rast, Ypts$N.Rate, fun=mean)
rastYld <- rasterize(sp, rast, Ypts$Yld, fun=mean)

## Plot previous year yield
png(paste(field.name, "_prev_yield.png", sep=""), width=5, height=10, units='in', res=100)
plot(rastYldPrev, main=paste("Observed Yield ", prev.yr), col=rev(heat.colors(100)),
     legend.args=list(text='kg/ha', side=4, font=2, line=2.5, cex=1.0))
dev.off()

## Plot current year yield
png(paste(field.name, "_yield.png", sep=""), width=5, height=10, units='in', res=100)
plot(rastYld, main="Observed Yield 2017", col=rev(heat.colors(100)),
     legend.args=list(text='kg/ha', side=4, font=2, line=2.5, cex=1.0))
dev.off()

## Plot as-applied nitrogen rates
png(paste(field.name, "_nitrogen_app.png", sep=""), width=5, height=10, units='in', res=100)
plot(rastN, main="Nitrogen Rate Experiment", col=rainbow(8),
     legend.args=list(text='kg-N/acre', side=4, font=2, line=2.5, cex=1.0))
dev.off()

#Reduce the data by removing random points from middle fert_rates so they do not have so much leverage on curve
plot(d$fert_rate,d$yield) #Before data is trimmed to give balance
mod0 <- lm(d$yield~d$fert_rate)
summary(mod0)
abline(mod0,lwd=2,lty=2,col=7)
## Remove random data points in center of cloud
d1 <- subset(d,d$fert_rate<140)
d2 <- subset(d,d$fert_rate>280)
d3 <- subset(d,d$fert_rate>=140 & d$fert_rate<=280) 
d4 <- d3[sample(nrow(d3),dim(d2)[1] + dim(d1)[1]),]
d5 <- rbind(d1,d2,d4)
plot(d5$fert_rate,d5$yield)
d <- d5

write.csv(d, "Oeste25_reduced_data.csv")  #Save the reduced data so that it can be consistent for adding 
#                                            spatial parameter and model comparisons
#####################################################

##### spataial parameter calc ######
library(RANN)

d <- read.csv("Oeste25_reduced_data.csv")
d <- d[,-1]
Ypts <- data.frame(x=d$x, y=d$y, YldPrev=d$prev_yl, Yld=d$yield, N.Rate=d$fert_rate)   

spCalc <- function(x,y){
  w <- x
  y1 <- y[[1]]
  y2 <- y[[2]]
  y2[y2==0] <- 0.0000000001
  lf <- ncol(y1)-1
  a <- 2
  
  # set up data
  for(i in 2:ncol(y2)){
    w[paste('nD',i-1, sep="")] <- y2[,i]
  }
  for(i in 2:ncol(y1)){
    w[paste('n',i-1, sep="")] <- y1[,i]
  }
  w[,(ncol(w)+1):(ncol(w)+(lf))] <- 0
  colnames(w)[((ncol(w))-(lf-1)):ncol(w)] <- paste('nVal',1:(lf),sep="")
  if(any(grepl("yl", colnames(w))|grepl("Yl",colnames(w)))){
    for(j in ((ncol(w))-(lf-1)):ncol(w)){
      for(i in 1:nrow(w)){
        w[i,j] <- w[(w[i,j-(lf)]),('YldPrev')]
      }
    }
  }else{
    for(j in ((ncol(w))-(lf-1)):ncol(w)){
      for(i in 1:nrow(w)){
        w[i,j] <- w[(w[i,j-(lf)]),('Pro')]
      }
    }
  }
  
  # calculate sp
  w[,(ncol(w)+1):(ncol(w)+(lf))] <- 0
  colnames(w)[((ncol(w))-(lf-1)):ncol(w)] <- paste('nNum',1:(lf),sep="")
  if(any(grepl("yl", colnames(w))|grepl("Yl",colnames(w)))){
    for(j in ((ncol(w))-(lf-1)):ncol(w)){
      for(i in 1:nrow(w)){
        w[i,j] <- (w[i,j-(lf*3)]^-a)*(w[i,j-lf]-w[i,'YldPrev'])
      }
    }
  }else{
    for(j in ((ncol(w))-(lf-1)):ncol(w)){
      for(i in 1:nrow(w)){
        w[i,j] <- (w[i,j-(lf*3)]^-a)*(w[i,j-lf]-w[i,'Pro'])
      }
    }
  }
  w[,ncol(w)+1] <- 0
  colnames(w)[ncol(w)] <- "sp"
  if(lf==1){
    if(any(grepl("yl", colnames(w))|grepl("Yl",colnames(w)))){
      for(i in 1:nrow(w)){
        w[i,"sp"] <- (w[i,(ncol(w)-lf)]/(
          w[i,(ncol(w)-(lf*4))]^-a))
      }
    }else{
      for(i in 1:nrow(w)){
        w[i,"sp"] <- (w[i,(ncol(w)-lf)]/(
          w[i,(ncol(w)-(lf*4))]^-a))
      }
    }
  }else{
    if(any(grepl("yl", colnames(w))|grepl("Yl",colnames(w)))){
      for(i in 1:nrow(w)){
        w[i,"sp"] <- (apply(
          w[i,c((ncol(w)-(lf)):(ncol(w)-1))],1,sum)/apply(
            w[i,c(((ncol(w)-((lf*4)))):((ncol(w)-((lf*3)+1))))]^-a,1,sum))
      }
    }else{
      for(i in 1:nrow(w)){
        w[i,"sp"] <- (apply(
          w[i,c((ncol(w)-(lf)):(ncol(w)-1))],1,sum)/apply(
            w[i,c(((ncol(w)-((lf*4)))):((ncol(w)-((lf*3)+1))))]^-a,1,sum))
      }
    }
  }
  
  return(w)
}

pc <- proc.time()
dPts <- Ypts 
dPtsNbrs <- nn2(data = cbind(dPts$x,dPts$y), query = cbind(dPts$x,dPts$y), k = 11)  # gets nn ID and dist
dPts <- spCalc(dPts,dPtsNbrs) 
d$diff_nn_yld <- dPts$sp
(proc.time()-pc)/60

write.table(d, "d_sp.txt", sep="\t")
d <- read.table("d_sp.txt", sep="\t")

###### Predictive Models ############################

#Linear regression with reduced data set
plot(d$fert_rate,d$yield,ylab="Grain yield (kg/ha)",xlab="N fertilizer rate (kg/ha)")
mod0 <- lm(d$yield~d$fert_rate)
summary(mod0)
abline(mod0,lty=2)
apprxR2 <- 1 - (sum(residuals(mod0)^2)/sum((d$yield-mean(d$yield))^2))  #Calculated approximate R^2
apprxR2

d$fertsq <- d$fert_rate^2
mod1 <- lm(d$yield~d$fert_rate + d$fertsq + d$prev_yl + d$elev_m + d$ec30 + d$ec90 + d$diff_nn_yld,data=d)
summary(mod1)
AIC(mod1)
points(d$fert_rate,predict.lm(mod1),col=2,cex=.5)
cor(d$ec30,d$ec90) #could probably drop ec30 because highly correlated with ec90
length(which(d$fert_rate==0))
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fit non-linear predictive models
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Define non-linear model function
hyperbolic <- function(alpha, beta, gamma, N, tao){
  ### Hyperbolic function ###
  ## alpha => baseline effect (value when nitrogen=0)
  ## beta => max fertilizer effect (asymptote)
  ## gamma => fertilizer efficiency (how quickly the curve rises)
  ## N => nitrogen applied
  y <- alpha + ((beta - alpha) * N)/((1.0/gamma) + N) + tao
  return(y)
}
## -------------------------------------------------------------------
##************ Yield model *******************
#---------------------------------------------------------------------
# 
b1=0
gamma=0.001
t1=0
mod2 <- nls(yield ~ hyperbolic(a0 +
                                   a1*elev_m +
                                   a2*ec30 +
                                   a3*ec90, 
                               b0 +
                                   b1*prev_yl,
                                 gamma, fert_rate,
                                 t1*diff_nn_yld),
                                 
            data=d, nls.control(maxiter = 1000, minFactor=1e-10),
            start=list(a0=-1, a1=10, a2=-1, a3=-1, 
                       b0=7000)) #
summary(mod2)
AIC(mod2)
apprxR2 <- 1 - (sum(residuals(mod2)^2)/sum((d$yield-mean(d$yield))^2))
apprxR2

## Notes: 
##  AIC = 237860.8                                   
## -------------------------------------------------
## Visualize scatter using all variables
# First make a data frame with regression coefficients and there std. errors
# If some variables do not come into the model this data frame must isert the variable 
# in their appropriate position in list and give their values 0 so when the Rcpp
# functions use the variables the calculations will come out correct and all of the
# equations do not need to be changed for each field.
Ycm1 <- summary(mod2)$coefficients[1:5,1:2] #Yield equation coefficients
b1 <- c(0,0)
gamma <- c(0.001,0.0001)
Ycm <-rbind(Ycm1,b1,gamma)
Ycm

Yield.prev = mean(d$prev_yl) ## 
#Aspect.cos = mean(d$aspect_cos)
#Aspect.sin = mean(d$aspect_sin)
#TPI = mean(d$tpi30)
Elev = mean(d$elev_m)
#Slope = mean(d$slope_deg)
#NDVI.2016 = mean(d$ndvi_2016)
#NDVI.2017 = mean(d$ndvi_2017) 
#NDVI.2015 = mean(d$ndvi_2015)
#NDVI.2014 = mean(d$ndvi_2014)
#NDVI.2013 = mean(d$ndvi_2013)
#NN.Yld = mean(d$diff_nn_yld)
EC.30 = mean(d$ec30)
EC.90 = mean(d$ec90)

Dtmp <- data.frame(Nrate=seq(0,500))
Dtmp$alpha = Ycm[1,1] + Ycm[2,1]*Elev + Ycm[3,1]*EC.30 + Ycm[4,1]*EC.90
Dtmp$beta = Ycm[5,1] + Ycm[6,1]*Yield.prev
gamma=Ycm[7,1]
tao=0
Dtmp$yld <- hyperbolic(Dtmp$alpha, Dtmp$beta, gamma, Dtmp$Nrate,tao)
summary(Dtmp$yld)
#Show mean predicted yields on the observed yields
plot(d$fert_rate, d$yield, ylab="Yield (kg/ha)", xlab="Nitrogen (kg/ha)")  ## plotted with 2016 nitrogen inputs, d$yl_nn2016
lines(Dtmp$Nrate,Dtmp$yld,col=2,lwd=2,lty=1)
#Add predicted yield for each point using all vriables as points on top of observed values
d$alpha = Ycm[1,1] + Ycm[2,1]*d$elev_m + Ycm[3,1]*d$ec30 + Ycm[4,1]*d$ec90
d$beta = Ycm[5,1] + Ycm[6,1]*d$prev_yl
gamma = Ycm[7,1]
d$predYield <- hyperbolic(d$alpha, d$beta, gamma, d$fert_rate,tao)      ## 2017 nitrogen inputs
summary(d$predYield)  #which(d$predYield<=0)
d$predYield <- ifelse(d$predYield<0,0,d$predYield)
points(d$fert_rate,d$predYield, pch=19, cex=.5, col=2)     ## 2016 nitrogen inputs
legend(400,4000,c("observed","predicted"),cex=.8,pch=c(1,19),col=c(1,2))

##*************** Sigmoidal model **********************
## Define sigmoidal model function 
sigmoidal <- function(alpha, Beta, delta, gamma, N, tao){
  #Beta is max yield
  y <- (alpha + (Beta - alpha)/(1 + exp(delta-gamma*N))) + tao 
  return(y)
}
## -------------------------------------------------------------------
##************ Yield model *******************
#---------------------------------------------------------------------
# 
gamma=.02
mod3 <- nls(yield ~ sigmoidal(a0 +
                                a1*prev_yl +
                                a2*ec30 +
                                a3*ec90 +
                                a4*elev_m, 
                              Beta,
                              delta,
                              gamma, 
                              fert_rate,
                              t1*diff_nn_yld),
            data=d, nls.control(maxiter = 500, minFactor=1e-10),
            start=list(a0=150,a1=.1,a2=-7,a3=.8,a4=1.5,
                       Beta=60, delta=6,t1=3)) #
summary(mod3)
AIC(mod3)
apprxR2 <- 1 - (sum(residuals(mod3)^2)/sum((d$yield-mean(d$yield))^2))
apprxR2

## AIC = 237339.6   ~R^2 = 0.1267063
##*********************************************************************
Ycm1 <- summary(mod3)$coefficients[1:7,1:2] #Yield equation coefficients
gamma <- c(gamma,.0001)
t1 <- summary(mod3)$coefficients[8,1:2]
Ycm <- rbind(Ycm1,gamma,t1)
Ycm

Yield.2017 = mean(d$prev_yl) ## 
EC.30 = mean(d$ec30)
EC.90 = mean(d$ec90)
Elev = mean(d$elev_m)
NN.Yld = mean(d$diff_nn_yld)

Dtmp <- data.frame(Nrate=seq(0,1500))
Dtmp$alpha = Ycm[1,1] + Ycm[2,1]*Yield.2017 + Ycm[3,1]*EC.30 + Ycm[4,1]*EC.90 + Ycm[5,1]*Elev 
Beta = Ycm[6,1] 
delta = Ycm[7,1]
gamma=Ycm[8,1]
tao= Ycm[9,1]*NN.Yld 
Dtmp$yld <- sigmoidal(Dtmp$alpha, Beta, delta, gamma, Dtmp$Nrate, tao)
summary(Dtmp$yld)
#Show mean predicted yields on the observed yields
#png(paste(field.name, "_predYld.png", sep=""), width=5, height=5, units='in', res=100)
plot(d$fert_rate, d$yield, ylab="Yield (kg/ha)", xlab="Nitrogen (kg/ha)")  ##
lines(Dtmp$Nrate,Dtmp$yld,col=2,lwd=2,lty=1)
#Add predicted yield for each point using all vriables as points on top of observed values
d$alpha = Ycm[1,1] + Ycm[2,1]*d$prev_yl + Ycm[3,1]*d$ec30 + Ycm[4,1]*d$ec90 + Ycm[5,1]*d$elev_m
Beta = Ycm[6,1] 
delta = Ycm[7,1]
gamma=Ycm[8,1]
tao= Ycm[9,1]*d$diff_nn_yld 
d$Yield <- sigmoidal(d$alpha, Beta, delta, gamma, d$fert_rate, tao)      
#d$Yield <- ifelse(d$n_lbs_ac>100,30,hyperbolic(alpha,beta,gamma,d$n_lbs_ac)) #Put upper asymptote on
summary(d$Yield)  #which(d$Yield<=0)
d$Yield <- ifelse(d$Yield<0,0,d$Yield)
points(d$fert_rate,d$Yield, pch=19, cex=.5, col=2)     
legend(1000,4000,c("observed","predicted"),cex=.8,pch=c(1,19),col=c(1,2))

## Write yield equation parameter values to a data file
write.table(Ycm, "YieldParms_oeste25.txt", sep="\t")

#Export data file with predicted yield 
write.table(d, "Predicted_yield_Oeste25.txt", sep="\t")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Net Return Calculation @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#*******************************************************************************************************
#The predicted Yield and protein values from above remain associated with yield points (= most dense data)
# for net return claculation below. 

library(Rcpp)

#################### Run the code below to compile the C++ functions #############
#Create an Rcpp function that calculates the crop yield, protein and net return given economic parameters and
# returns the profit maximizing nirtogen fertilizer rate and the highest average mean net return ($/acre) for
# the site-specific VRA, the no N-fert application and the producers chosen uniform N-fert application rate.
cppFunction('
            NumericMatrix NRoptC(NumericMatrix dnr, NumericVector Ycc, double Bp, 
            double CN, double ssAC, double FC, double fcnr,
            int Nratecutoff, int rr){
            NumericMatrix NRopt(rr,8);
            double alphay;
            double betay;
            double deltay;
            double gammay;
            double taoy;
            double yld;
            double NR;
            double NRmax = 0;
            int Nrateopt;
            double Ynrmax;
            double NR0;
            double NRfcnr;
            for(int jj = 0; jj < rr; jj++){
            for(int Nrate = 0; Nrate <= Nratecutoff; Nrate++){ 
            alphay = Ycc(0) + Ycc(1)*dnr(jj,4) + Ycc(2)*dnr(jj,5) + Ycc(3)*dnr(jj,6) + Ycc(4)*dnr(jj,7);
            betay = Ycc(5); 
            deltay = Ycc(6);
            gammay = Ycc(7);
            taoy = Ycc(8)*dnr(jj,8);
            yld = (alphay + (betay - alphay)/(1 + exp(deltay-gammay*Nrate))) + taoy;
            NR = (yld*Bp) - CN*Nrate - FC; 
            if(Nrate==0){
            NRmax = NR;
            Nrateopt = 0;
            Ynrmax = yld;
            NR0 = NR;
            }else{
            if(NR-ssAC > NRmax){
            NRmax = NR-ssAC;
            Nrateopt = Nrate;
            Ynrmax = yld;
            }    
            }
            if(Nrate==fcnr){
            NRfcnr = NR; 
            }
            }
            NRopt(jj,0) = Bp;
            NRopt(jj,1) = dnr(jj,0);
            NRopt(jj,2) = dnr(jj,1);
            NRopt(jj,3) = Ynrmax;
            NRopt(jj,4) = Nrateopt;
            NRopt(jj,5) = NRmax;
            NRopt(jj,6) = NR0;
            NRopt(jj,7) = NRfcnr;
            }
            return NRopt;
            }
            ')
## End of site-specific optimization

#Create an Rcpp function that calculates the crop yield, protein and net return given economic parameters and
# returns the profit maximizing nirtogen fertilizer rate for full field uniform application. 
cppFunction('
            NumericMatrix NRoptffC(NumericMatrix dnr, NumericVector Ycc,
            double Bp, double CN, double ssAC, double FC, 
            int Nratecutoff, int rr){
            NumericMatrix NRff(Nratecutoff+1,2);
            double alphay;
            double betay;
            double deltay;
            double gammay;
            double taoy;
            double yld;
            double NR;
            double NRffld;
            double YldTot;
            for(int Nrate = 0; Nrate <= Nratecutoff; Nrate++){
            NRffld = 0;
            YldTot = 0;
            for(int jj = 0; jj < rr; jj++){
            alphay = Ycc(0) + Ycc(1)*dnr(jj,4) + Ycc(2)*dnr(jj,5) + Ycc(3)*dnr(jj,6) + Ycc(4)*dnr(jj,7);
            betay = Ycc(5); 
            deltay = Ycc(6);
            gammay = Ycc(7);
            taoy = Ycc(8)*dnr(jj,8);
            yld = (alphay + (betay - alphay)/(1 + exp(deltay-gammay*Nrate))) + taoy;
            NR = (yld*Bp) - CN*Nrate - FC - ssAC;  
            NRffld = NRffld + NR;
            YldTot = YldTot + yld;
            }
            NRff(Nrate,0) = Nrate;
            NRff(Nrate,1) = NRffld; 
            } 
            return NRff;
            } 
            ')
## End full-field optimization

## End function declaration *********************************************************************************

########################################
##***********************************************************************************************************

#Bring in yield point data with all independent variables and predicted values for each point
DD <- read.table("Predicted_yield_oeste25.txt", sep="\t")
names(DD)
dim(DD)

#Get regression model parameter files
Ycm <- read.table("YieldParms_oeste25.txt", sep="\t")
Ycm

## ******************************************************************************************
## User can vary these parameters below 
#--------------------------------------------------------------------------------------------

#These are variables that can be varied to get estimate of probability of net return
#Base Price (Bp) for winter wheat in $/bu
#Prices received for winter wheat and cost of N in last 16 years to create variation
Prc <- read.csv("MT_Organic_vs_Conv_wheat_N_prices.csv", header = T) #**********
#hist(Prc$HRWWconv)

field.name <- "oeste_25"
field.res <-  70*0.3048    # resolution in meters
fieldsize = 1642.27 # oeste_25 acres   = 664.68 hectares #************
## Costs of production
#CN = #Cost of nitrogen fertilizerin $/lb for side-dress application
ssAC = 2.00 #Site-specific application additional cost  #********** 
FC = 71.31 #All other Fixed/ownership costs. #********* 
fcnr = 180 #names(which.max(table(d$fert_rate)))    #The nitrogen rate that the farmer was going to apply uniformly to the field #**********
Nratecutoff = 500 # ## Set the upper asymptote for protein and yield if curves do not flatten appropriately
if(max(DD$fert_rate)<Nratecutoff){
  Nratecutoff=max(DD$fert_rate)
} 
Nratecutoff <- as.integer(as.numeric(Nratecutoff))
DDs <- subset(DD,DD$fert_rate>=Nratecutoff)
yasymptot = mean(DDs$Yield) #Yield maximum (asymptot). Artificially imposes a maximum yield value when fit curve does not flatten out.

#Premium or dockage to add to BP
#PD<- read.table("Billings_PremDock_2016.csv",header=T,sep=",") #*********
#colnames(PD)<-c("PCpro","Bill","pro","PremDock")
#PD$prosq = PD$pro^2
#fm = lm(PD$PremDock~PD$pro+PD$prosq) #Create linear estimate function for protein premium/dockage
#summary(fm)
#B0pd = as.vector(coef(fm)[1])
#B1pd = as.vector(coef(fm)[2])
#B2pd = as.vector(coef(fm)[3])
#PD$PredPrem=B0pd + B1pd*PD$pro + B2pd*PD$prosq
#lines(PD$pro,PD$PredPrem,lwd=2)

##--------------------------------------------------------------------------------------
# Set up data and parameter matrices that can go into Rcpp functions
DNR <- as.data.frame(DD) #Create a new data frame to get rid of character columns that mess with Rcpp
#Parameter values for yield prediction equation
Yc <- data.frame(B0=Ycm[1,1], By1=Ycm[2,1], B4=Ycm[3,1], B5=Ycm[4,1], B6=Ycm[5,1],Beta=Ycm[6,1],delta=Ycm[7,1],gamma=Ycm[8,1],tao=Ycm[9,1])
Ycc <- as.vector(as.numeric(Yc))

dnr <- matrix(0,dim(DNR)[1],9)
DNR <- as.matrix(DNR)
#head(DNR)
dnr[,1:2] <- as.numeric(DNR[,1:2]) #column numbers vary b/w data sets
dnr[,3] <- as.numeric(DNR[,4])
dnr[,4] <- as.numeric(DNR[,6])
dnr[,5] <- as.numeric(DNR[,8])
dnr[,6:7] <- as.numeric(DNR[,10:11])
dnr[,8] <- as.numeric(DNR[,7])
dnr[,9] <- as.numeric(DNR[,12])

#colnames(dnr) <- c("x","y","yield","fert_rate","prev_yl","ec30","ec90","elev","diff_nn_yld")
#dnr <- dnr[1:1000,]
rr = as.integer(dim(dnr)[1])

#@@@@@@@@@@@@@@@@@@@@@@
## metric to imperial (for initial analysis)
# 2.2 lbs/kg
# 2.47 ac/ha
# 60 lbs/bu
Prc[,2:3] <- (Prc[,2:3]*2.2)/60 # prices of wheat = $/kg
Prc[,4] <- Prc[,4]*2.2 # price N = $/kg

#@@@@@@@@@@@@@@@@@@@@@@
##*******************************************************************************
# Determine the variation in net return for each nitrogen application strategy
# by using the variation in price received for WW ($/bu) from 2000 to 2016
#********************************************************************************
pc <- proc.time() #This determines the time it takes to run this loop below
sPr = 100 #Number of simulations with different prices
Bp.var <- matrix(0,nrow=sPr,ncol=7)
colnames(Bp.var) <- c("BaseP","NR.ssopt","NR.0","NR.fs","ffopt.Nrate","NR.ffopt", "NR.act")
for(bp in 1:sPr){
  rp <- as.integer(runif(1,1,length(Prc$Year)))
  Prc$Year[rp]
  Bp=as.numeric(Prc[rp,'HRWWconv'])
  CN=as.numeric(Prc[rp,'Ncost'])
  
  #Send to C++ function to get net return on ss, 0, fs, and org nitrogen fertilizer application strategies
  NRopt <- NRoptC(dnr, Ycc, Bp, CN, ssAC, FC, fcnr, Nratecutoff, rr)
  colnames(NRopt) <- c("BaseP","x","y","yld.NRmax","N.rate.opt","NR.opt","NR.0","NR.fcnr")
  Bp.var[bp,'BaseP'] <- Bp
  Bp.var[bp,'NR.ssopt'] <- mean(NRopt[,'NR.opt'])
  Bp.var[bp,'NR.0'] <- mean(NRopt[,'NR.0'])
  Bp.var[bp,'NR.fs'] <- mean(NRopt[,'NR.fcnr'])
  
  #Send to C++ function to get full-field uniform nitrogen rate that maximizes net return
  NRff <- NRoptffC(dnr, Ycc, Bp, CN, ssAC, FC, Nratecutoff, rr)
  colnames(NRff) <- c("N.rate","NR.ffopt")
  NRffmax <- subset(NRff, NRff[,'NR.ffopt']==max(NRff[,'NR.ffopt']))
  Bp.var[bp,'ffopt.Nrate'] <- NRffmax[,'N.rate']
  Bp.var[bp,'NR.ffopt'] <- NRffmax[,'NR.ffopt']/rr #Full-field profit maximizing
  Bp.var[bp,'NR.act'] = sum((dnr[,3]*Bp) - CN*dnr[,4] - FC)/rr #Actual net return using the rate experiment
}
png(paste(field.name, "_avgNR_box.png", sep=""), width=7, height=7, units='in', res=100)
boxplot(Bp.var[,'NR.ssopt'],Bp.var[,'NR.0'],Bp.var[,'NR.fs'],Bp.var[,'NR.ffopt'],
        Bp.var[,'NR.act'], ylab="Average net return ($/ha)", 
        names=c("NR.ssopt","NR.0","NR.fs","NR.ffopt","NR.act"),col=3)
dev.off()
(proc.time()-pc)/60 #minutes for simulation to run

#calculate the propabability that site-specific application is more profitable than the other app methods
rt = 1000 #number of replicate tests
tot23=0
tot24=0
tot25=0
tot26=0
for(nn in 1:rt){
  rp <- as.integer(runif(1,1,length(Prc$HRWWconv)))
  tot23 <- ifelse(Bp.var[rp,'NR.ssopt']>Bp.var[rp,'NR.0'],1+tot23,0+tot23)
  tot24 <- ifelse(Bp.var[rp,'NR.ssopt']>Bp.var[rp,'NR.fs'],1+tot24,0+tot24)
  tot25 <- ifelse(Bp.var[rp,'NR.ssopt']>Bp.var[rp,'NR.ffopt'],1+tot25,0+tot25)
  tot26 <- ifelse(Bp.var[rp,'NR.ssopt']>Bp.var[rp,'NR.act'],1+tot26,0+tot26)
}
Table.p <- data.frame(tot23/rt, tot24/rt, tot25/rt, tot26/rt)
colnames(Table.p) <- c("NR.ssopt>NR.0","NR.ssopt>NR.fs","NR.ssopt>NR.ffopt","NR.ssopt>NR.act")
row.names(Table.p) <- c("Probability")
Table.p #Probability that site-specific optimum will be greater than other N application strategies
#write.table(Table.p,"TableP.csv")

#Present results as a bar chart with error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
#Average net return over the field for the different N application methods without organic
TF3 <- data.frame(NR.ssopt=mean(Bp.var[,'NR.ssopt']), NR.0=mean(Bp.var[,'NR.0']), 
                  NR.fs=mean(Bp.var[,'NR.fs']), NR.ffopt=mean(Bp.var[,'NR.ffopt'])) 
TF3
#If plotting metric values don't forget to relabel y-axis in barplot
#Put std. error bars on
TF3.means <- as.matrix(TF3)
TF3.sd <- as.matrix(data.frame(NR.ssopt=sd(Bp.var[,'NR.ssopt']), NR.0=sd(Bp.var[,'NR.0']), 
                               NR.fs=sd(Bp.var[,'NR.fs']), NR.ffopt=sd(Bp.var[,'NR.ffopt'])))

png(paste(field.name, "_avgNR_bar.png", sep=""), width=7, height=7, units='in', res=100)
barTF <- barplot(TF3.means, ylim=c(0,max(TF3)+20), col=3, axis.lty=1, ylab="Average net return ($/ha)",
                 names=c("NR.ssopt","NR.0","NR.fs","NR.ffopt"))
error.bar(barTF, TF3.means, 1.96*TF3.sd/10)
text(.7, TF3$NR.ssopt+2, as.integer(TF3$NR.ssopt))
text(1.9,TF3$NR.0+2, as.integer(TF3$NR.0))
text(3.1,TF3$NR.fs+2, as.integer(TF3$NR.fs))
text(4.3,TF3$NR.ffopt+2, as.integer(TF3$NR.ffopt))
dev.off()

#**************************************************************************************************************
## Total nitrogen fertilizer used with each strategy
## This requires running through the NRopt one more time to obtain
## output with 2016 price and cost
##*******************************************************************
rp <- length(Prc$HRWWconv)
Bp=as.numeric(Prc[rp,'HRWWconv'])
CN=as.numeric(Prc[rp,'Ncost'])
#Send to C++ function to get net return on ss, 0, fs, and org nitrogen fertilizer application strategies
NRopt <- NRoptC(dnr, Ycc, Bp, CN, ssAC, FC, fcnr, Nratecutoff, rr)
colnames(NRopt) <- c("BaseP","x","y","yld.NRmax","N.rate.opt","NR.opt","NR.0","NR.fcnr")

#Send to C++ function to get full-field uniform nitrogen rate that maximizes net return
NRff <- NRoptffC(dnr, Ycc, Bp, CN, ssAC, FC, Nratecutoff, rr)
colnames(NRff) <- c("N.rate","NR.ffopt")
NRffmax <- subset(NRff, NRff[,'NR.ffopt']==max(NRff[,'NR.ffopt']))  

summary(NRopt[,'N.rate.opt'])
TF4 <- data.frame(N.ssopt=(sum(NRopt[,'N.rate.opt'])*fieldsize)/dim(NRopt)[1], N.0=0,N.fs=fcnr*fieldsize,
                  N.ffopt=NRffmax[,'N.rate']*fieldsize) 
TF4
counts <- as.matrix(TF4)
ymax=max(TF4)+1
ymin=min(TF4)
#If plotting metric values don't forget to relabel yaxis
png(paste(field.name, "_NfertUsed_bar.png", sep=""), width=7, height=7, units='in', res=100)
barplot(counts,ylim=c(ymin,ymax),col=4,ylab="Nitrogen used on field (kg/field)",
        main="Compare nitrogen used with each application strategy",names.arg=c("N.ssopt","N.0","N.fs","N.ffopt"))
text(.7, TF4$N.ssopt-5000, as.integer(TF4$N.ssopt),col=0)
text(1.9,TF4$N.0-5000, as.integer(TF4$N.0),col=0)
text(3.1,TF4$N.fs-5000, as.integer(TF4$N.fs),col=0)
text(4.3,TF4$N.ffopt-5000, as.integer(TF4$N.ffopt),col=0)
#Percent decrease or increase in nitrogen use with site-specific profit maximizing over the full-field optimum rate
if(as.integer(TF4$N.ffopt)==0){
  text(1.3,TF4$N.ssopt+200, as.integer(((TF4$N.ssopt-TF4$N.ffopt)/1)))
  text(2.9,TF4$N.ssopt+200,"% more N used in N.ssopt than N.ffopt")
}else{
  if(as.integer((1-TF4$N.ssopt/(TF4$N.ffopt)) * 100)<0){    ## N.ffopt = 0 -> creates NA
    text(1.3,TF4$N.ssopt-1000, as.integer((1-TF4$N.ssopt/(TF4$N.ffopt)) * -100))
    text(2.9,TF4$N.ssopt-1000,"% more N used in N.ssopt than N.ffopt")
  }else{
    text(1.3,TF4$N.ssopt-100, as.integer((1-TF4$N.ssopt/TF4$N.ffopt) * 100))
    text(2.5,TF4$N.ssopt-100,"% less N used in N.ssopt than N.ffopt", cex=.8)
  }
}
dev.off()

#Plot the distribution of N rates across the field in the ssopt approach
png(paste(field.name, "_NoptRates_hist.png", sep=""), width=7, height=7, units='in', res=100)
hist(NRopt[,'N.rate.opt'],xlab="lbs/acre",main="Site-specific profit maximizing nitrogen rates",col=4)
points(mean(NRopt[,'N.rate.opt']),10,pch=19,col=2,cex=2) #Mean site-specific optimum N rate = red dot
dev.off()

#Plot the distribution and mean of full field uniform application profit maximizing N rate
mean(Bp.var[,'ffopt.Nrate'])
sd(Bp.var[,'ffopt.Nrate'])
#write.table(data.frame(mean = mean(Bp.var[,'ffopt.Nrate']), sd = sd(Bp.var[,'ffopt.Nrate'])), "ffopt.csv")
png(paste(field.name, "_ffNoptRates_hist.png", sep=""), width=7, height=7, units='in', res=100)
hist(Bp.var[,'ffopt.Nrate'],col=4, xlab="profit maximizing top-dress nitrogen rate (lbs/acre)",main="Variation in best full-field uniform N rate over years")
points(mean(Bp.var[,'ffopt.Nrate']),0,pch=19,col=2, cex=2) # Mean full field optimum N rate = yellow dot
dev.off()
#************************************************************************************************************

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Profit maximizing site-specific nitrogen fertilizer rate prescription map
## ...model used N application for 2017
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Clean up data for net return and nitrogen top-dress prescription maps
D <- as.data.frame(NRopt)
NRatepts <- data.frame(x=D$x,y=D$y,OptN.Rate=D$N.rate.opt,NRopt=D$NR.opt,NR.0=D$NR.0,NR.fcnr=D$NR.fcnr)
abs.max <- max(NRatepts$OptN.Rate)
abs.min <- min(NRatepts$OptN.Rate)
cuts <- seq(0, round(abs.max, -1), 10)

## Convert opt nitrogen rate data to raster
Ymat <- as.matrix(NRatepts)
sp <- SpatialPoints(coords=Ymat[, c('x', 'y')])
e <- extent(NRatepts[, c('x', 'y')])
rast <- raster(ext=e, resolution=10)
rastOptN <- rasterize(sp, rast, NRatepts$OptN.Rate, fun=mean)

## Plot optimum (profit maximizing) nitrogen rate for current year
png(paste(field.name, "_opt_N_prescrip.png", sep=""), width=5, height=10, units='in', res=100)
plot(rastOptN, main=paste("Profit maximizing nitrogen rates 2017"), col=rev(rainbow(8)),
     legend.args=list(text='kg/ha', side=4, font=2, line=2.5, cex=1.0))
dev.off()

##-----------------------------------------------------------------------------------------
## Net return map 2017 with profit maximizing site specific nitrogen application in 2016
D <- as.data.frame(NRopt)
NRatepts <- data.frame(x=D$x,y=D$y,OptN.Rate=D$N.rate.opt,NRopt=D$NR.opt,NR.0=D$NR.0,NR.fcnr=D$NR.fcnr)
abs.max <- max(NRatepts$NRopt)
abs.min <- min(NRatepts$NRopt)
cuts <- seq(0, round(abs.max, -1), 10)

## Convert opt nitrogen rate data to raster
Ymat <- as.matrix(NRatepts)
sp <- SpatialPoints(coords=Ymat[, c('x', 'y')])
e <- extent(NRatepts[, c('x', 'y')])
rast <- raster(ext=e, resolution=10)
rastOptNR <- rasterize(sp, rast, NRatepts$NRopt, fun=mean)

## Plot optimum (profit maximizing) nitrogen rate for current year
png(paste(field.name, "_ssopt_NR.png", sep=""), width=5, height=10, units='in', res=100)
plot(rastOptNR, main=paste("Net return with ss opt N rate 2017"), col=rev(topo.colors(100)),
     legend.args=list(text='$/ha', side=4, font=2, line=2.5, cex=1.0))
dev.off()

##-----------------------------------------------------------------------------------------
## Actual net return map
Bp <- as.numeric(Prc[93,'HRWWconv'])
D <- data.frame(x=dnr[,1],y=dnr[,2],Yld=dnr[,3],N.rate=dnr[,4],NR.act=(dnr[,3]*Bp) - Prc[93,'Ncost']*dnr[,4] - FC)
abs.max <- max(D$NR.act)
abs.min <- min(D$NR.act)
cuts <- seq(0, round(abs.max, -1), 10)

Ymat <- as.matrix(D)
sp <- SpatialPoints(coords=Ymat[, c('x', 'y')])
e <- extent(D[, c('x', 'y')])
rast <- raster(ext=e, resolution=10)
rastNR <- rasterize(sp, rast, D$NR.act, fun=mean)

## Plot net return for current year
png(paste(field.name, "_NR_actual.png", sep=""), width=5, height=10, units='in', res=100)
plot(rastNR, main=paste("Net return actual 2017"), col=rev(topo.colors(100)),
     legend.args=list(text='$/ha', side=4, font=2, line=2.5, cex=1.0))
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
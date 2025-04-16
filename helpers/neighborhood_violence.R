# OG author: Chris Dharma

library(haven)
library(dplyr)
library(fastDummies)
library(scales)
library(truncnorm)
library(tidyverse)

set.seed(123)

N = 1000

#Function to create a beta distribution continuous variable
betadist <- function(min,max,shape1,shape2,N) {
  rbetasample <- rbeta(N, shape1 = shape1, shape2 = shape2)
  return(scales::rescale(rbetasample,to=c(min,max)))
} 

# Create the number of communities
community_max=100
#truncate it to keep them as integers
community_code<-trunc(betadist(min=1,max=community_max,shape1=2,shape2=15,N=N))
hist(community_code)
#The frequencies within each community is skewed, some communities have a large number of people, some have smaller

#Create the community level variables using beta distribution 
Income<-trunc(betadist(min=12015,max=241300,shape1=2,shape2=10,N=community_max))
OwnerVac<-betadist(min=0,max=27.5,shape1=1,shape2=10,N=community_max)
RenterVac<-betadist(min=0,max=29.00,shape1=1,shape=5,N=community_max)
edu<-betadist(min=0,max=75.4,shape1=2,shape=5,N=community_max)
PubAssist<-betadist(min=0,max=20.40,shape1=2,shape=9,N=community_max)
Poor<-betadist(min=1.10,max=58,shape1=2,shape=9,N=community_max)
Unemployed<-betadist(min=2.30,max=34.60,shape1=2,shape=9,N=community_max)
Vehicles0<-betadist(min=0,max=75.2,shape1=1,shape=19,N=community_max)
Marital1<-betadist(min=7.9,max=93.7,shape1=3,shape=7,N=community_max) #%never married
Marital3<-betadist(min=0,max=7.6,shape1=2,shape=4,N=community_max) #%separated

##Predict the crime level (exposure) according to other variables at the community level, later will call this A
x.pred.exp<-scale(cbind(Income,OwnerVac,RenterVac,edu,PubAssist,Poor,Unemployed,Vehicles0,Marital1,Marital3))
x.matrix.A<-as.matrix(cbind(rep(1,community_max),x.pred.exp))
coef.A<-as.matrix(c(32.64,-5,4.1,2.5,-5.6,10.5,10.1,10.2,2.1,0.4,0.1))
error_com=rnorm(community_max,0,1)
pred.A <- (x.matrix.A %*% coef.A) + error_com
summary(pred.A)
hist(pred.A)
#Original data = 32.64 ( 2.03 - 223.2)
##Add by a factor of 21 to ensure that there were no negative values. Just wanted to follow the distribution
Ascore <- pred.A + 21
summary(Ascore)
hist(Ascore)

A <- Ascore

ARIMA.pred.12moavg_yr<-c()
MedianHHInc<-c()
OwnerVacRate<-c()
RenterVacRate<-c()
pEd_LThs<-c()

pctNumHHPubAssist<-c()
pctPoor<-c()
pctUnemployedCLF<-c()
pctVehicles0<-c()
pmarital1<-c()
pmarital3<-c()

for (i in 1 : N) {
  ARIMA.pred.12moavg_yr[i]<-Ascore[community_code[i]]
  MedianHHInc[i]<-Income[community_code[i]]
  OwnerVacRate[i]<-OwnerVac[community_code[i]]
  RenterVacRate[i]<-RenterVac[community_code[i]]
  pEd_LThs[i]<-edu[community_code[i]]
  pctNumHHPubAssist[i]<-PubAssist[community_code[i]]
  pctPoor[i]<-Poor[community_code[i]]
  pctUnemployedCLF[i]<-Unemployed[community_code[i]]
  pctVehicles0[i]<-Vehicles0[community_code[i]]
  pmarital1[i]<-Marital1[community_code[i]]
  pmarital3[i]<-Marital3[community_code[i]]
}

#community variables
V <- cbind(MedianHHInc,OwnerVacRate,RenterVacRate,pEd_LThs,pctPoor,pctUnemployedCLF,pctVehicles0,pmarital1,pmarital3)

#Individual level variables:
mage<-trunc(betadist(min=12,max=61,shape1=2,shape=4,N=N))
mage.sq<-mage*mage
maxt<-betadist(min=48.28,max=87.19,shape1=7,shape=7,N=N)
summary(ARIMA.pred.12moavg_yr)
hist(ARIMA.pred.12moavg_yr)

concept.season<-sample(x=seq(1:4),size=N,prob=c(0.25,0.25,0.25,0.25),replace=T)
ed<-sample(x=seq(1:3),size=N,prob=c(0.30,0.45,0.25),replace=T)
ins.pdd<-sample(x=seq(1:3),size=N,prob=c(0.49,0.49,0.02),replace=T)
race.eth<-sample(x=seq(1:8),size=N,prob=c(0.26,0.05,0.005,0.12,0.005,0.02,0.50,0.04),replace=T)
concept.yr<-sample(x=seq(1:2),size=N,prob=c(0.75,0.25),replace=T)
parity.cat<-sample(x=seq(0:2),size=N,prob=c(0.40,0.30,0.30),replace=T)

#Individual variables
W <- cbind(ed,ins.pdd,concept.season,concept.yr,parity.cat,race.eth,mage,mage.sq,maxt)

##Create the two M variables (mediators, diabetes and preeclampsia) to be a function of A, W, and V 
x.pred<-scale(cbind(ARIMA.pred.12moavg_yr,pEd_LThs,mage,mage.sq,pctPoor,pctUnemployedCLF))
ed2<-ifelse(ed==2,1,0)
ed3<-ifelse(ed==3,1,0)
x.matrix<-as.matrix(cbind(rep(1,N),x.pred,ed2,ed3))

#Diabetes prev = 0.08
model.coef.diab<-c(log((0.08)/(1-0.08)),log(1.22),log(1.15),log(1.03),log(1.15),log(1.14),log(1.18),log(0.97),log(0.95))
x.matrix.diab<- x.matrix %*% as.matrix(model.coef.diab)
pi_x.diab<- 1 / (1 + exp(-x.matrix.diab))
pc_gdiabetes_c <- rbinom(n=N, size=1, prob=pi_x.diab)
summary(pc_gdiabetes_c)

#Preeclampsia prev = 0.06
model.coef.preeclampsia<-c(log((0.06)/(1-0.06)),log(1.10),log(1.25),log(1.04),log(1.20),log(1.14),log(1.10),log(0.98),log(0.95))
x.matrix.preeclampsia<- x.matrix %*% as.matrix(model.coef.preeclampsia)
pi_x.preeclampsia<- 1 / (1 + exp(-x.matrix.preeclampsia))
pc_preeclampsia_c <- rbinom(n=N, size=1, prob=pi_x.preeclampsia)
summary(pc_preeclampsia_c)

M <- cbind(pc_gdiabetes_c,pc_preeclampsia_c)

#Create outcomes here:

#Create the matrix of predictors
#make sure the continuous ones are on the same scale
x.pred<-scale(cbind(ARIMA.pred.12moavg_yr,pEd_LThs,mage,mage.sq,pctPoor,pctUnemployedCLF))
x.matrix<-as.matrix(cbind(rep(1,N),x.pred,pc_preeclampsia_c,pc_gdiabetes_c,ed2,ed3))

#Intercept should be the overall average, so the invlogit of the coefficient should be the mean for that variable (ptb) = 0.07, or logit(0.07) = intercept
model.coef.ptb<-c(log((0.07)/(1-0.07)),log(1.22),log(1.15),log(1.03),log(1.15),log(1.20),log(1.25),log(1.30),log(1.21),log(0.97),log(0.95))
x.matrix.ptb<- x.matrix %*% as.matrix(model.coef.ptb)
pi_x.ptb<- 1 / (1 + exp(-x.matrix.ptb))
ptb <- rbinom(n=N, size=1, prob=pi_x.ptb)
summary(ptb)

Y1 <- ptb

dat<-as.data.frame(cbind(V, W, ARIMA.pred.12moavg_yr, M, Y,community_code))

use <- select(dat, 
              owner_vacancy_rate = OwnerVacRate, 
              renter_vacancy_rate = RenterVacRate,
              #public_assit_prop = PubAssist,
              poverty_line_prop = pctPoor, 
              unemployed_prop = pctUnemployedCLF,
              no_car_prop = pctVehicles0, 
              never_married_prop = pmarital1, 
              separated_prop = pmarital3, 
              maternal_age = mage, 
              edu_level = ed, 
              insurance_status = ins.pdd, 
              race_ethnicity = race.eth, 
              parity = parity.cat, 
              neighborhood_violence = ARIMA.pred.12moavg_yr,
              preterm_birth = Y1, 
              community_code)

use <- mutate(use, 
              across(all_of(c("edu_level", "insurance_status", "race_ethnicity", "parity")), as.factor))

write.csv(use, "data/crime_toy.csv", row.names = F)

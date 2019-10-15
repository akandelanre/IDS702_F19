
###########################################################################
###########################################################################
####################### The mercury in bass analysis ######################
###########################################################################
###########################################################################

###### Clear environment and load libraries
rm(list = ls())
library(ggplot2)
library(lme4) #for the lmer function
library(lattice)

###### Load the data
bass <- read.csv("data/mercurydata.txt", header = T)
bass$station <- factor(bass$station)
bass$river <- factor(bass$river)

###### View properties of the data
dim(bass)
summary(bass)
table(bass$station) #we don't have enough data in some stations

#### Exploratory data analysis
#first, mercury by itself
ggplot(bass,aes(mercury)) +
  geom_histogram(alpha=.8,fill=rainbow(10),bins=10)
#not so normal. let's try log
ggplot(bass,aes(log(mercury))) +
  geom_histogram(alpha=.8,fill=rainbow(10),bins=10)
#hmmm...skewed the other way now
#how about sqrt
ggplot(bass,aes(sqrt(mercury))) +
  geom_histogram(alpha=.8,fill=rainbow(10),bins=10)
#this looks better.
#We will continue with the untransformed variable first

#let's look at mercury vs. weight and vs. length
river_col <- c("red4","blue4")
pairs(mercury~ weight+length,data=bass,pch=3,col=river_col[as.numeric(bass$river)])
#similar relationships for mercury vs. length and mercury vs. weight
#river does not seem to make a difference

#check correlations among the predictors to look for colinearity
cor(bass$weight, bass$length)
#quite correlated! We may have to drop one of them
#Let's mean center both for now
bass$length_c <- bass$length - mean(bass$length)
bass$weight_c <- bass$weight - mean(bass$weight)

#now let's actually separate the relationships out by river
#first weight
ggplot(bass,aes(x=weight, y=mercury)) +
  geom_point(alpha = .5,colour="blue4") +
  geom_smooth(method="lm",col="red3") +
  labs(title="Mercury vs Weight") +
  facet_wrap(~river,ncol=2)
#different lines for each river. 
#Try with sqrt(mercury)
ggplot(bass,aes(x=weight, y=sqrt(mercury))) +
  geom_point(alpha = .5,colour="blue4") +
  geom_smooth(method="lm",col="red3") +
  labs(title="Mercury vs Weight") +
  facet_wrap(~river,ncol=2)
#looks like we may have a problem with linearity if we use the sqrt transformation
#we really should be controlling for the stations so let's split by stations
ggplot(bass,aes(x=weight, y=mercury)) +
  geom_point(alpha = .5,colour="blue4") +
  geom_smooth(method="lm",col="red3") +
  labs(title="Mercury vs Weight") +
  facet_wrap(~station,ncol=4)
#maybe some different slopes, but sample sizes are too small to tell

#let's look at length by river
ggplot(bass,aes(x=length, y=mercury)) +
  geom_point(alpha = .5,colour="blue4") +
  geom_smooth(method="lm",col="red3") +
  labs(title="Mercury vs Length") +
  facet_wrap(~river,ncol=2)
#less obvious here
#how about by station
ggplot(bass,aes(x=length, y=mercury)) +
  geom_point(alpha = .5,colour="blue4") +
  geom_smooth(method="lm",col="red3") +
  labs(title="Mercury vs Weight") +
  facet_wrap(~station,ncol=4)
#againmaybe some different slopes, but sample sizes are too small to tell


#Let's try a simple linear model first
mercreg <- lm(mercury ~ length_c + weight_c, data = bass)
car::vif(mercreg)
#highly correlated
#try a linear model using interactions with station
mercreg <- lm(mercury ~ length_c*station + weight_c, data = bass)
car::vif(mercreg)
#Let's just drop one of length and weight
#lets double check to see the stronger predictor
pairs(mercury~ weight+length,data=bass,pch=3,col=river_col[as.numeric(bass$river)])
#let's keep length
mercreg <- lm(mercury ~ length_c*station, data = bass)

#diagnostics
plot(y = mercreg$resid, x=bass$length, xlab = "Length", ylab = "Residual")
abline(0,0)
boxplot(mercreg$resid~bass$station, xlab = "Station", ylab = "Residual")
#pretty good fit... let's go with it. you could try logs as well.
#look at results.
summary(mercreg)

#see if interaction effects are useful
mercregnoint <- lm(mercury ~ length_c + station, data = bass)
anova(mercregnoint, mercreg)
#they seem to be useful overall, but many are based on very small sample sizes.
#it would be good to borrow strength across the stations to get better estimates
#that is what a hierarchical model does!



#### let's do a hierarchical regression (random effects regression) with lmer
#this call just uses a random intercept
mercreglmerint <- lmer(mercury ~ length_c + (1 | station), data = bass) 
summary(mercreglmerint)
confint(mercreglmerint)
#look at the intercepts (and the common slope) for each station
coef(mercreglmerint)
#these equal the fixed effects plus the random effect
fixef(mercreglmerint)
ranef(mercreglmerint)
#we can plot the different lines for each station
newdata <- data.frame(station=rep(unique(bass$station),each=20),
                      length_c=seq(min(bass$length_c),max(bass$length_c),20))
newdata$station <- as.factor(newdata$station)
newdata$pred <- predict(mercreglmerint,newdata=newdata,type="response")
ggplot(newdata, aes(x = length_c, y = pred, colour = station)) +
  geom_line(size = 0.7) +
  labs(x = "Length (Centered)", y = "Predicted Mercury") +
  theme(legend.position="none") +
  geom_text(data = newdata,aes(label = station), hjust = 0.5, vjust = 1)


#Let's add a random slope
mercreglmerintslope <- lmer(mercury ~ length_c + ( 1 + length_c  | station), data = bass) 
summary(mercreglmerintslope)
#the intercepts and slope for each station
coef(mercreglmerintslope)
fixef(mercreglmerintslope)
ranef(mercreglmerintslope)
#we can plot the different lines for each station
newdata$pred <- predict(mercreglmerintslope,newdata=newdata,type="response")
ggplot(newdata, aes(x = length_c, y = pred, colour = station)) +
  geom_line(size = 0.7) +
  labs(x = "Length (Centered)", y = "Predicted Mercury") +
  theme(legend.position="none") +
  geom_text(data = newdata,aes(label = station), hjust = 0.5, vjust = 1)

#the predicted values of mercury for each bass
preds <- fitted(mercreglmerintslope)

#residual plots
plot(mercreglmerintslope)
plot(y = residuals(mercreglmerintslope), x = bass$length, xlab= "Length", ylab = "Residuals")
#looks pretty good -- model seems to fit the data well.

#how about plots of residuals by station?
ggplot(bass,aes(x=length, y=residuals(mercreglmerintslope))) +
  geom_point(alpha = .5,colour="blue4") +
  geom_smooth(method="lm",col="red3") +
  labs(title="Mercury vs Length") +
  facet_wrap(~station,ncol=5)
xyplot(residuals(mercreglmerintslope) ~ length | station, data = bass)
#also reasonable


#Let's compare models using AIC and BIC
AIC(mercreglmerint) #old model
AIC(mercreglmerintslope) #new model
BIC(mercreglmerint) #old model
BIC(mercreglmerintslope) #new model
#new model results in a reduction of both




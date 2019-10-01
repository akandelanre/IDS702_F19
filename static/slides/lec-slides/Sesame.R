###########################################################################
###########################################################################
########################## Sesame street analysis #########################
###########################################################################
###########################################################################

###### Clear environment and load libraries
rm(list = ls())
library(arm)
library(pROC)
library(e1071)
library(caret)
library(nnet)
library(knitr)



###### Load the data
sesame = read.csv("data/sesame.txt", header = T)



###### View properties of the data
sesame$sex <- factor(sesame$sex,levels=c(1,2),labels=c("Male","Female"))
sesame$setting <- factor(sesame$setting,levels=c(1,2),labels=c("Home","School"))
#let's make a dummy variable with encouragement = 2 as the baseline, since we want to interpret the effect of encouragement
sesame$viewenc <- factor(sesame$viewenc,levels=c(2,1),labels=c("NotEncouraged","Encouraged"))
str(sesame)
head(sesame)
#viewcat is the outcome variable
#1=rarely watched the show; 2=once or twice a week; 
#3=three to five times a week; 4=watched the show on average more than 5 times a week
#viewcat is actually ordered, but we will treat it as unordered for this analysis.
dim(sesame)
table(sesame$viewcat) #Almost evenly spread out



###### Exploratory data analysis
#let's look at associations of viewcat with various predictors.
#Ignore all the post-test variables
#let's start looking at plots with continuous predictors
boxplot(prenumb ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Pretest Numbers",col=rainbow(10))
boxplot(prelet ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Pretest Letters",col=rainbow(10))
boxplot(preform ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Pretest Forms",col=rainbow(10))
boxplot(preclasf ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Pretest Classification",col=rainbow(10))
boxplot(prerelat ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Pretest Relationships",col=rainbow(10))
boxplot(prebody ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Pretest Body Parts",col=rainbow(10))
boxplot(age ~ viewcat, data = sesame, xlab = "Viewing Category", ylab = "Age",col=rainbow(10))

#now some categorical predictors
table(sesame$viewcat, sesame$viewenc)
prop.table(table(sesame$viewcat, sesame$viewenc), 2)
#prop.table(table(sesame$viewcat, sesame$viewenc), 1) gives the row percentages
chisq.test(table(sesame$viewcat, sesame$viewenc))
table(sesame$viewcat, sesame$site)
prop.table(table(sesame$viewcat, sesame$site), 2)
chisq.test(table(sesame$viewcat, sesame$site))

#center all the continuous predictors (for interpretation) for the model
sesame$cprenumb <- sesame$prenumb - mean(sesame$prenumb)
sesame$cprelet <- sesame$prelet - mean(sesame$prelet)
sesame$cpreform <- sesame$preform - mean(sesame$preform)
sesame$cpreclasf <- sesame$preclasf - mean(sesame$preclasf)
sesame$cprerelat <- sesame$prerelat - mean(sesame$prerelat)
sesame$cprebody <- sesame$prebody - mean(sesame$prebody)
sesame$cage <- sesame$age - mean(sesame$age)

#we will not explore interactions for this data
#we may not have enough to see any meaningful interactions or even fit them



###### Model fitting
#fit a multinomial regression.  note that viewcat = 1 is the reference level.
viewcatreg1 <- multinom(viewcat ~ viewenc + as.factor(site) + cprenumb + cprelet + cpreform 
                       + cpreclasf + cprerelat + cprebody + cage, data=sesame)
summary(viewcatreg1)
#exponentiate to get to interpretations in terms of multiplicative factors for odds
exp(coef(viewcatreg1))
#get confidence intervals
confint(viewcatreg1)
exp(confint(viewcatreg1))
#notice that there are no p-values but we can use the CIs for inference
output1 <- summary(viewcatreg1)
z_value <- output1$coefficients/output1$standard.errors
p_value <- (1 - pnorm(abs(z_value), 0, 1))*2 
#we are using two-tailed z test, that is, a normal approximation
full_summary1 <- lapply(c(2:4),function(x) rbind(output1$coefficients[as.character(x),],
                                                 output1$standard.errors[as.character(x),],
                                                 z_value[as.character(x),],
                                                 p_value[as.character(x),]))
kable(lapply(full_summary1,function(x) {rownames(x) <- c("Coefficient","Std. Errors","z-value","p-value"); x}))
#too many p-values to check simulataneously, let's use deviance test instead



###### Deviance test
#Let's test if site is a useful  predictor using a change in deviance test
viewcatreg1nosite <- multinom(viewcat ~ viewenc + cprenumb + cprelet + cpreform + cpreclasf 
                              + cprerelat + cprebody  + cage, data = sesame)
anova(viewcatreg1, viewcatreg1nosite, test = "Chisq")
#p-value is small -- site looks like a useful predictor

#test if encouragement is a useful  predictor using a change in deviance test
viewcatreg1noenc <- multinom(viewcat ~ as.factor(site) + cprenumb + cprelet + cpreform 
                             + cpreclasf + cprerelat + cprebody + cage, data = sesame)
anova(viewcatreg1, viewcatreg1noenc, test = "Chisq")
#p-value is small -- encouragement looks like a useful predictor


#test if setting is a useful  predictor using a change in deviance test
viewcatreg1withsetting <- multinom(viewcat ~ viewenc + as.factor(site) + cprenumb + cprelet + cpreform 
                                   + cpreclasf + cprerelat + cprebody + cage + setting, data=sesame)
anova(viewcatreg1withsetting, viewcatreg1, test = "Chisq")
#p-value is large -- looks like setting is not a useful predictor



###### Predictions
#predicted probabilities for cases in the model
predprobs <- fitted(viewcatreg1) 
#look at first five rows just to see what results
predprobs[1:5,]

##interpreting results with predictions
#let's get predicted probabilities for someone at average values of continuous predictors and in site 1
newdata <- sesame[1:2,]
newdata$site <- 1
newdata$cprenumb <- 0
newdata$cprelet <- 0
newdata$cprerelat <- 0
newdata$cpreform <- 0
newdata$cprebody <- 0
newdata$cpreclasf <- 0
newdata$cage <- 0
newdata$viewenc[1] <- "Encouraged"
newdata$viewenc[2] <- "NotEncouraged"
newdata

predict(viewcatreg1, newdata, type = "probs")
#           1         2         3         4
#1 0.02655388 0.3858047 0.3838143 0.2038272
#2 0.27441129 0.2144506 0.3088348 0.2023033
 
#now let's repeat for a kid in site 2
newdata$site <- 2
predict(viewcatreg1, newdata, type = "probs")
#           1         2         3         4
#1 0.02561089 0.2577973 0.3219774 0.3946144
#2 0.24999038 0.1353515 0.2447120 0.3699461




###### Diagnostics
####diagnostics comparing average raw residuals across bins based on predictor values
#for viewcat = 1:  create a raw residual using only the first column of the predicted probabilities
rawresid1 <- (sesame$viewcat == 1) -  predprobs[,1]

#for viewcat = 2:  create a raw residual using only the second column of the predicted probabilities
rawresid2 <- (sesame$viewcat == 2) -  predprobs[,2]

#for viewcat = 3:  create a raw residual using only the third column of the predicted probabilities
rawresid3 <- (sesame$viewcat == 3) -  predprobs[,3]

#for viewcat = 4:  create a raw residual using only the fourth column of the predicted probabilities
rawresid4 <- (sesame$viewcat == 4) -  predprobs[,4]

##can do binned plots for continuous variables
#make a 2 by 2 graphical display
par(mfcol = c(2,2))
binnedplot(sesame$cprenumb, rawresid1, xlab = "Prenumbers", ylab = "Raw residuals", main = "Binned plot: viewcat = 1")
binnedplot(sesame$cprenumb, rawresid2, xlab = "Prenumbers", ylab = "Raw residuals", main = "Binned plot: viewcat = 2")
binnedplot(sesame$cprenumb, rawresid3, xlab = "Prenumbers", ylab = "Raw residuals", main = "Binned plot: viewcat = 3")
binnedplot(sesame$cprenumb, rawresid4, xlab = "Prenumbers", ylab = "Raw residuals", main = "Binned plot: viewcat = 4")

#repeat for the other continuous predictors.....  all looks okay!

##if you want to change the reference level, e.g., using level 4 as reference, use the relevel command
sesame$viewcat2 <- relevel(as.factor(sesame$viewcat), ref = "4")

#then refit the model
viewcatreg2 = multinom(viewcat2 ~ viewenc + as.factor(site) + cprenumb + cprelet + cpreform + cpreclasf + cprerelat + cprebody + cage, data = sesame)
summary(viewcatreg2)

##tables for factor variables
tapply(rawresid1, sesame$site, mean)
tapply(rawresid2, sesame$site, mean)
tapply(rawresid3, sesame$site, mean)
tapply(rawresid4, sesame$site, mean)

tapply(rawresid1, sesame$viewenc, mean)
tapply(rawresid2, sesame$viewenc, mean)
tapply(rawresid3, sesame$viewenc, mean)
tapply(rawresid4, sesame$viewenc, mean)


## Accuracy
pred_classes <- predict(viewcatreg1)
Conf_mat <- confusionMatrix(as.factor(pred_classes),as.factor(sesame$viewcat))
Conf_mat$table
Conf_mat$overall["Accuracy"];
Conf_mat$byClass[,c("Sensitivity","Specificity")]


## Individual ROC curves for the different levels
par(mfcol = c(2,2))
roc((sesame$viewcat==1),predprobs[,1],plot=T,print.thres="best",legacy.axes=T,print.auc =T,
    col="red3",percent=T,main="Group 1")
roc((sesame$viewcat==2),predprobs[,2],plot=T,print.thres="best",legacy.axes=T,print.auc =T,
    col="gray3",percent=T,main="Group 2")
roc((sesame$viewcat==3),predprobs[,2],plot=T,print.thres="best",legacy.axes=T,print.auc =T,
    col="green3",percent=T,main="Group 3")
roc((sesame$viewcat==4),predprobs[,2],plot=T,print.thres="best",legacy.axes=T,print.auc =T,
    col="blue3",percent=T,main="Group 4")
#we can also combine them into a single plot


## Multi-class ROC curve (average of all pairwise comparisons)
par(mfcol = c(3,4))
multiclass.roc(sesame$viewcat,predprobs,plot=T,print.thres="best",legacy.axes=T,print.auc =T,col="red3",percent=T)

###########################################################################
###########################################################################
### Play around with the data some more and see if you can do better!!! ###
###########################################################################
###########################################################################




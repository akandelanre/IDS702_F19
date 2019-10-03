###########################################################################
###########################################################################
########################## Political ideology analysis #########################
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
library(MASS)



###### Load the data
political <- read.table("data/political.txt", header = T,stringsAsFactors = T)



###### View properties of the data
political$Ideology <- ordered(political$Ideology,levels=c("Very Liberal","Slightly Liberal","Moderate","Slightly Conservative","Very Conservative"))
str(political)
head(political)
dim(political)
table(political$Ideology) #we definitely have more data points at the moderate level compared to other levels
table(political)


###### Exploratory data analysis

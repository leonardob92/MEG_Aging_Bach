## ANOVA TESTING AGE, POTENTIALLY CONFOUNDING VARIABLES AND BRAIN DATA

#install packages
install.packages("lme4")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("haven")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("broom")
install.packages("readxl")
install.packages("jmv")
install.packages("effectsize")
install.packages("car")
install.packages("heplots")

#Load packages
library(lme4)
library(tidyr)         #data management package from the tidyverse suite.
library(ggplot2)       #plot package from the tidyverse suite.
library(tidyverse)     #for data manipulation and visualization
library(haven)         #
library(ggpubr)        #for creating easily publication ready plots
library(rstatix)       #for easy pipe-friendly statistical analyses
library(broom)         #for printing a summary of statistical tests as data frames
library(readxl)        #for reading excel files
library(jmv)           #for manova analyses
library(effectsize)    #for calculating effect sizes
library(lmerTest)
library(car)
library(heplots)

########

condition = 2 #1 = previously memorized musical sequences; 2 = novel T1; 3 = novel T3

########

setwd("/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG")
#Load data
tsa <- read_excel("R_Old_NewT1_NewT3.xlsx",col_names = TRUE, sheet = condition)
tsa <- tsa[tsa$Subject != 35,] #removing subject 35 which did not have proper MEG data for at least 2 out of 3 conditions
tsa <- tsa[tsa$Subject != 1,] #removing subject 1 who was 38 years old (so they did not belong to any age group)
tsa <- na.omit(tsa) #removing NaN values

#Summary statistics
young <- tsa[tsa$age==1,]
elderly <- tsa[tsa$age==2,]

summary(young)
#years of education
mean(young$`years education`,na.rm=TRUE) #omitting NA scores
sd(young$`years education`,na.rm=TRUE) #omitting NA scores
#WM
mean(young$WM,na.rm=TRUE) #omitting NA scores
sd(young$WM,na.rm=TRUE) #omitting NA scores

summary(elderly)
#years of education
mean(elderly$`years education`,na.rm=TRUE) #omitting NA scores
sd(elderly$`years education`,na.rm=TRUE) #omitting NA scores
#WM
mean(elderly$WM,na.rm=TRUE) #omitting NA scores
sd(elderly$WM,na.rm=TRUE) #omitting NA scores

###### MANCOVA

dependent <- cbind(tsa$LAC, tsa$LHIT, tsa$LIFG, tsa$MC, tsa$RAC, tsa$RHIT, tsa$RIFG, tsa$VMPFC)
age <- tsa$age
mancova_model <- manova(dependent ~ age + tsa$WM + tsa$`years education` + tsa$`years musical training` + tsa$sex) #manova + mancova
summary(mancova_model, test = "Wilks")
summary.aov(mancova_model) #univariate tests
eta_squared(mancova_model, partial = TRUE) #effect size

# Check multivariate normality assumption
qqnorm(residuals(mancova_model))
qqline(residuals(mancova_model))


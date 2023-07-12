################################################################################

#MULTIVARIATE ANALYSIS ON BEHAVIOURAL DATA OF YOUNG AND OLDER ADULTS

################################################################################

#LOAD DATA

setwd("/Users/au685269/Library/CloudStorage/OneDrive-AarhusUniversitet/AdditionalProjects/TemporalSequencesAge/Behavioral_MEG/") #working directory

library(readxl) #for reading excel files

recog_age <- read_excel("TSA_Recognition.xlsx") #original file
recog_age <- recog_age[-c(1,35),] #remove participants 1 and 35

#create data frame with response variables and covariates
resp <- subset(recog_age,select = c(Subject,Group,Old_Cor,New_T1_Cor,New_T3_Cor,
                                    Years_Education,WM,Years_Musical_Training,
                                    Sex))

#create data frame with reaction time variables and covariates
rtime <- subset(recog_age,select = c(Subject,Group,Old_Cor_RT,New_T1_Cor_RT,
                                     New_T3_Cor_RT,Years_Education,WM,
                                     Years_Musical_Training,Sex))

################################################################################

#DESCRIPTIVE STATISTICS AND DATA VISUALIZATION

summary(recog_age) #descriptive statistics for each variable

library(ggpubr) #for creating plots

#boxplot with response variables
ggboxplot(resp, x = "Group",
          y = c("Old_Cor","New_T1_Cor","New_T3_Cor"),
          merge = TRUE,
          ylab = "Number of correct responses")

#boxplot with reaction time variables
ggboxplot(rtime, x = "Group",
          y = c("Old_Cor_RT","New_T1_Cor_RT","New_T3_Cor_RT"),
          merge = TRUE,
          ylab = "Reaction time (ms)")

################################################################################

#MULTIVARIATE AND UNIVARIATE OUTLIERS
#checking that there are no multivariate or univariate outliers

library(rstatix) #for statistical analyses
options(scipen=999) #remove scientific notation

#checking for multivariate outliers
#Mahalanobis Distance measures the distance between two points in multivariate 
#space and is used to locate outliers (p-value < .001)
resp$mahalanobis<- mahalanobis(resp[,3:5], colMeans(resp[,3:5]), 
                               cov(resp[,3:5])) #calculate mahalanobis distance
resp$pvalue <- pchisq(resp$mahalanobis, df=3, lower.tail=FALSE) #report p-values
rtime$mahalanobis<- mahalanobis(rtime[,3:5], colMeans(rtime[,3:5]), 
                                cov(rtime[,3:5])) #calculate mahalanobis distance
rtime$pvalue <- pchisq(rtime$mahalanobis, df=3, lower.tail=FALSE) #report p-values

#checking for univariate outliers
identify_outliers(resp,Old_Cor)
identify_outliers(resp,New_T1_Cor)
identify_outliers(resp,New_T3_Cor)
identify_outliers(rtime,Old_Cor_RT)
identify_outliers(rtime,New_T1_Cor_RT)
identify_outliers(rtime,New_T3_Cor_RT)

#creating a new data frame without multivariate outliers
resp$Old_Cor[resp$Old_Cor == 2] <- NA #replacing participant 32's score with NA
resp$New_T1_Cor[resp$New_T1_Cor == 2] <- NA #replacing participant 32's score with NA
resp$New_T3_Cor[resp$New_T3_Cor <= 5] <- NA #replacing participant 25, 32, 77 and 78's scores with NA

################################################################################

#HOMOCEDASTICITY ASSUMPTION
#checking whether the variances are homogeneous

library(car) #use to perform tests, create visualizations and transform data

leveneTest(Old_Cor ~ Group, resp)
leveneTest(New_T1_Cor ~ Group,resp)
leveneTest(New_T3_Cor ~ Group,resp)
leveneTest(Old_Cor_RT ~ Group,rtime)
leveneTest(New_T1_Cor_RT ~ Group,rtime)
leveneTest(New_T3_Cor_RT ~ Group,rtime)

################################################################################

#MULTIVARIATE NORMALITY ASSUMPTION
#difficult to check this assumption, so we check for univariate normality instead

library(ggpubr)

#Density plots
#density plots provide a visual judgment about whether the distribution is bell shaped

#response variables
ggdensity(resp$Old_Cor)
ggdensity(resp$New_T1_Cor)
ggdensity(resp$New_T3_Cor)

#reaction time variables
ggdensity(rtime$Old_Cor_RT)
ggdensity(rtime$New_T1_Cor_RT)
ggdensity(rtime$New_T3_Cor_RT)

#Q-Q plots
#q-q plots draw the correlation between a given sample and the normal 
#distribution (a 45-degree reference line is also plotted)

#response variables
ggqqplot(resp$Old_Cor)
ggqqplot(resp$New_T1_Cor)
ggqqplot(resp$New_T3_Cor)

#reaction time variables
ggqqplot(rtime$Old_Cor_RT)
ggqqplot(rtime$New_T1_Cor_RT)
ggqqplot(rtime$New_T3_Cor_RT)

#Shapiro-Wilk test (test of normality)
shapiro.test(resp$Old_Cor)
shapiro.test(resp$New_T1_Cor)
shapiro.test(resp$New_T3_Cor)
shapiro.test(rtime$Old_Cor_RT)
shapiro.test(rtime$New_T1_Cor_RT)
shapiro.test(rtime$New_T3_Cor_RT)

################################################################################

#LINEARITY ASSUMPTION
#check whether there is a linear relationship between each pair of dependent variables for each group of the independent variable

plot(resp,pch = 20)
plot(rtime,pch = 20)

################################################################################

#HOMOGENEITY OF VARIANCE-COVARIANCE MATRICES

#Box M's test is extremely sensitive to sample size and departures from 
#normality

library(rstatix) #for statistical analyses

box_m(resp[,3:5],as.factor(resp$Group))
box_m(rtime[,3:5],as.factor(rtime$Group))

################################################################################

#MANOVA

library(jmv) #for manova analyses
library(effectsize) #for calculating effect sizes

#response
dependent <- cbind(resp$Old_Cor,resp$New_T1_Cor,resp$New_T3_Cor) #dependent variables
independent <- resp$Group #independent variable
manova_model <- manova(dependent ~ independent, data = resp) #computing the model
summary(manova_model) #summary of model
eta_squared(manova_model) #effect size
summary.aov(manova_model) #post-hoc analyses

#reaction time
dependent <- cbind(rtime$Old_Cor_RT,rtime$New_T1_Cor_RT,rtime$New_T3_Cor_RT)
independent <- rtime$Group
manova_model <- manova(dependent ~ independent, data = rtime) #computing the model
summary(manova_model) #summary of model
eta_squared(manova_model) #effect size
summary.aov(manova_model) #post-hoc analyses (ANOVA) #alpha-level = .05/3 = .017

################################################################################

#MULTIVARIATE KRUSKAL-WALLIS
#non-parametric alternative to MANOVA

install.packages("remotes")
remotes::install_github("jacobmaugoust/ULT") #reference: https://github.com/jacobmaugoust/ULT
library(ULT)

#complete dataset
multkw(as.factor(resp$Group), resp[,3:5])
#without outliers
multkw(as.factor(resp$Group), resp[,6:8])

#POST-HOC TESTS
#complete dataset
kruskal_test(resp, resp$Old_Cor ~ resp$Group)
kruskal_test(resp, resp$New_T1_Cor ~ resp$Group)
kruskal_test(resp, resp$New_T3_Cor ~ resp$Group)
#without outliers
kruskal_test(resp, resp$Old_Cor_NA ~ resp$Group)
kruskal_test(resp, resp$New_T1_Cor_NA ~ resp$Group)
kruskal_test(resp, resp$New_T3_Cor_NA ~ resp$Group)

################################################################################

#MANCOVA (response variables)

#dependent variables
dependent <- cbind(resp$Old_Cor,resp$New_T1_Cor,resp$New_T3_Cor)

#independent variable
independent <- resp$Group

#calculate MANCOVA
mancova_model <- manova(dependent ~ independent + resp$Years_Education + 
                          resp$WM + resp$Years_Musical_Training + resp$Sex) 
summary(mancova_model, test = "Wilks") #use Wilk's lambda
eta_squared(mancova_model, partial = TRUE) #effect size
summary.aov(mancova_model) #univariate tests
qqnorm(residuals(mancova_model)) #q-q plot of residuals
qqline(residuals(mancova_model)) #add q-q line

################################################################################

#MANCOVA (reaction time variables)

#dependent variables
dependent <- cbind(rtime$Old_Cor_RT,rtime$New_T1_Cor_RT,rtime$New_T3_Cor_RT)

#independent variable
independent <- rtime$Group

#calculate MANCOVA
mancova_model <- manova(dependent ~ independent + rtime$Years_Education + 
                          rtime$WM + rtime$Years_Musical_Training + rtime$Sex)
summary(mancova_model, test = "Wilks") #use Wilk's lambda
eta_squared(mancova_model, partial = TRUE) #effect size
summary.aov(mancova_model) #univariate tests
qqnorm(residuals(mancova_model)) #q-q plot of residuals
qqline(residuals(mancova_model)) #add q-q line

################################################################################

#RAINCLOUD PLOTS

#source: Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., van Langen, J.,
#& Kievit, R. A.Raincloud plots: a multi-platform tool for robust data visualization 
#[version 2; peer review: 2 approved] Wellcome Open Research 2021, 4:63. 
#https://doi.org/10.12688/wellcomeopenres.15191.2

if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('jorvlan/raincloudplots')

library(raincloudplots)
library(ggplot2)

#Response variables
#subsetting the data by group
young <- resp[resp$Group == "Y",] #young participants
older <- resp[resp$Group == "O",] #older participants

#for the function to work, all variables must be the same length
y_old <- young$Old_Cor
o_old <- older$Old_Cor
length(y_old) <- length(o_old)
y_nt1 <- young$New_T1_Cor
o_nt1 <- older$New_T1_Cor
length(y_nt1) <- length(o_nt1)
y_nt3 <- young$New_T3_Cor
o_nt3 <- older$New_T3_Cor
length(y_nt3) <- length(o_nt3)

#creating the data frame for the plot
df_2x3 <- data_2x2(
  array_1 = y_old,
  array_2 = o_old,
  array_3 = y_nt1,
  array_4 = o_nt1,
  array_5 = y_nt3,
  array_6 = o_nt3,
  labels = (c('Young','Older')),
  jit_distance = .09,
  jit_seed = 321)

#colors for the plot
colors <- rep(c("#E33000", "#217CFF"), 3)

#run the function to create the plot
raincloud_2x3_repmes(data = df_2x3, colors = colors, fills = colors, size = 1, 
                     alpha = .5, ort = "v") +
  scale_x_continuous(breaks = c(1,2,3), limits = c(0.8, 4.3), 
                    labels = c("Memorized", "Novel T1", "Novel T3")) +
  ylab("Number of correct responses") +
  xlab("") +
  theme_bw() +
  theme(axis.title.y = element_blank())

#Reaction time variables
#subsetting the data by group
young <- rtime[rtime$Group == "Y",] #young participants
older <- rtime[rtime$Group == "O",] #older participants

#for the function to work, all variables must be the same length
y_old <- young$Old_Cor_RT
o_old <- older$Old_Cor_RT
length(y_old) <- length(o_old)
y_nt1 <- young$New_T1_Cor_RT
o_nt1 <- older$New_T1_Cor_RT
length(y_nt1) <- length(o_nt1)
y_nt3 <- young$New_T3_Cor_RT
o_nt3 <- older$New_T3_Cor_RT
length(y_nt3) <- length(o_nt3)

#creating the data frame for the plot
df_2x3 <- data_2x2(
  array_1 = y_old,
  array_2 = o_old,
  array_3 = y_nt1,
  array_4 = o_nt1,
  array_5 = y_nt3,
  array_6 = o_nt3,
  labels = (c('Young','Older')),
  jit_distance = .09,
  jit_seed = 321)

#colors for the plot
colors <- rep(c("#E33000", "#217AFF"), 3)

#run the function to create the plot
raincloud_2x3_repmes(data = df_2x3, colors = colors, fills = colors, size = 1, 
                     alpha = .5, ort = "v") +
  scale_x_continuous(breaks = c(1,2,3), limits = c(0.8, 4.3), 
                     labels = c("Memorized", "Novel T1", "Novel T3")) +
  ylab("Number of correct responses") +
  xlab("") +
  theme_bw() +
  theme(axis.title.y = element_blank())

################################################################################

#CORRELATION MATRIX

library(psych)
library(corrplot)

#calculate correlation and p-values
cor_resp <- cor(resp[3:8], use = "complete.obs")
cor_resp_p = cor.mtest(cor_resp, conf.level = 0.95)
cor_rtime <- cor(rtime[3:8], use = "complete.obs")
cor_rtime_p = cor.mtest(cor_rtime, conf.level = 0.95)

#create correlation matrix
corrplot(cor_resp, method = "color", diag = FALSE, tl.pos = "n", type = "lower",
         col = colorRampPalette(c("dodgerblue3","white","firebrick"))(200),cl.lim=c(-1,1), cl.length = 5,
         p.mat = cor_resp_p$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig")

corrplot(cor_rtime, method = "color", diag = FALSE, tl.pos = "n", type = "lower",
         col = colorRampPalette(c("dodgerblue3","white","firebrick"))(200),cl.lim=c(-1,1), cl.length = 5,
         p.mat = cor_rtime_p$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig")

#correlation plot
library(ggplot2)
gg <- ggplot(resp, aes(x = Old_Cor, y = Years_Education))
gg <- gg + geom_point()
gg <- gg + geom_smooth(alpha=0.3, method="lm")
print(gg)









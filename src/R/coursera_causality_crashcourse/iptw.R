#install pakcages (if needed)
install.packages("tableone")
install.packages("ipw")
install.packages("sandwich")
install.packages("survey")

# Load packages
library(tableone)
library(ipw)
library(sandwich) # for robust variance estimation
library(survey)

# Custom functions
expit <- function(x) {1/(1+exp(-x))}
logit <- function(p) {log(p)-log(1-p)}

# read in data
load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))
# View data
View(rhc)

# Treatment variabels is swang1
#X variables that we will use
#cat1: primary disease category
#age
#sex
#meanbp1: mean blood pressure

# Crearte a data set with just these variabels, for simplicity
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
COPD<-as.numeric(rhc$cat1=='COPD')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')
female<-as.numeric(rhc$sex=='Female')
died<-as.integer(rhc$death=='Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1=='RHC')
meanbp1<-rhc$meanbp1

# new dataset
mydata<-cbind(ARF, CHF, Cirr, colcan, Coma, lungcan, MOSF, sepsis, age, female, meanbp1, treatment, died)
mydata <- data.frame(mydata)

# Covariates we will use 
xvars <- c("age","female","meanbp1","ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF", "sepsis")

# Look at a table 1
table1 <- CreateTableOne(vars = xvars, 
                         strata = "treatment",
                         data = mydata,
                         test = FALSE)
## include standardized mean difference (SMD)
print(table1, smd = TRUE)

# propensity score model with logistic function
psmodel <- glm(treatment ~ age + female + meanbp1, + ARF + CHF + Cirr + colcan + Coma + lungcan + MOSF + sepsis,
               data = mydata,
               family = binomial(link = "logit"))

# Value of propensity score for each subject
ps <- predict(psmodel, type = "response")

# Create weights
# Vectorised ifelse statement to create weights
# Treated: 1/P(A=1|X)
# Control: 1/P(A=0|X)
weight <- ifelse(treatment ==1, 1/(ps), 1/(1-ps))

# apply weights to data
weighteddata <- svydesign(ids = ~1, # formula input. 1 for no clusters. Check ?svydesign for syntax
                          data = mydata, 
                          weights = ~weight) # formula input for computed weights

# Weighted Table 1
weightedtable <- svyCreateTableOne(vars = xvars,
                                   strata = "treatment",
                                   data = weighteddata,
                                   test = FALSE)

# ## Show table with SMD
print(weightedtable, smd = TRUE)

# to get a weighted mean for a single covariate directly:
# For example: age
mean(weight[treatment ==1]* age[treatment==1])/(mean(weight[treatment==1]))

# get causal risk difference
glm.obj <- glm(died ~ treatment, weights = weight, 
               family = quasibinomial(link="identity"))

summary(glm.obj)

betaiptw <- coef(glm.obj)
SE <- sqrt(diag(vcovHC(glm.obj, type = "HC0")))

causalrd <- (betaiptw[2])
lcl <- (betaiptw[2] - 1.96*SE[2])
ucl <- (betaiptw[2] + 1.96*SE[2])
c(lcl, causalrd, ucl)

# Get Causal relative risk. Weighted GLM
glm.obj <- glm(died ~ treatment, weights = weight, 
               family = quasibinomial(link = log))

summary(glm.obj)
betaiptw<-coef(glm.obj)
SE <- sqrt(diag(vcovHC(glm.obj, type = "HC0")))

causalrr <- exp(betaiptw[2])
lcl <- exp(betaiptw[2] - 1.96*SE[2])
ucl <- exp(betaiptw[2] + 1.96*SE[2])
c(lcl, causalrr, ucl)

#######################
# truncate weights at 10
#######################
truncweight <- replace(weight, weight >10, 10)

# Get causal risk difference
glm.obj <- glm(died ~ treatment, weights = truncweight,
               family = quasibinomial(link = "identity"))
summary(glm.obj)

betaiptw <- coef(glm.obj)
SE <- sqrt(diag(vcovHC(glm.obj, type = "HC0")))

causalrd <- (betaiptw[2])
lcl <- (betaiptw[2] - 1.96*SE[2])
ucl <- (betaiptw[2] + 1.96*SE[2])
c(lcl, causalrd, ucl)


#############################
#alternative: use ipw package
#############################

# Fit propensity scoring model to get weights
weightmodel <- ipwpoint(exposure = treatment, 
                        family = "binomial",
                        link = "logit",
                        denominator = ~ age + female + meanbp1+ARF+CHF+Cirr+colcan+
                          Coma+lungcan+MOSF+sepsis, data=mydata)

# numeric summary of weights
summary(weightmodel$ipw.weights)

# plot of weights
ipwplot(weights = weightmodel$ipw.weights,
        logscale = FALSE,
        main = "weights",
        xlim = c(0,22))
mydata$wt <- weightmodel$ipw.weights

# fit propensity score model to get weights but with truncation (percentile basis)
weightmodel <- ipwpoint(exposure = treatment,
                        family = "binomial",
                        link = "logit",
                        denominator= ~ age + female + meanbp1+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis,
                        data = mydata,
                        trunc = 0.01) # this will truncate 1st and 99th percentile for weights range

# umeric summary of weights
summary(weightmodel$weights.trun)

# plot of weights
ipwplot(weights = weightmodel$weights.trun, 
        logscale = FALSE,
        main = "weights",
        xlim = c(0,22))
mydata$wt <- weightmodel$weights.trun

# fit a marginal structural model (risk difference)
msm <- (svyglm(died ~ treatment, 
               design = svydesign(ids = ~1, 
                                  weights = ~wt,
                                  data = mydata)))
coef(msm)
confint(msm)

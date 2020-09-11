library(tableone)
library(Matching)

load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))

View(rhc)

#create a data set with just these variables, for simplicity
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
died<-as.numeric(rhc$death=='Yes')
age <- rhc$age
treatment <- as.numeric(rhc$swang1 == "RHC")
meanbp1 <- rhc$meanbp1

# newdataset
mydata <-cbind(ARF,CHF,Cirr, colcan, Coma, lungcan, MOSF, sepsis, age, female, meanbp1, treatment, died)
mydata<- data.frame(mydata)

# Covariates we will use
xvars <- c("ARF","CHF", "Cirr","colcan","Coma","lungcan", "MOSF", "sepsis", "age","female", "meanbp1")

# Look at a table 1
table1 <- CreateTableOne(vars = xvars, strata = "treatment", data=mydata, test=FALSE)
print(table1, smd = FALSE)

greedymatch <- Match(Tr = treatment, M = 1, X = mydata[xvars], replace = FALSE)
matched <- mydata[unlist(greedymatch[c("index.treated", "index.control")]), ]

matchedtab1 <-CreateTableOne(vars = xvars, strata = "treatment",
                             data = matched, test = FALSE)

print(matchedtab1, smd = TRUE)

# outcome analysis
y_trt <-matched$died[matched$treatment == 1]
y_con <-matched$died[matched$treatment == 0]

# pairwise difference
diffy <- y_trt - y_con

# paired t-test
t.test(diffy)

# McNemar test
table(y_trt, y_con)

mcnemar.test(matrix(c(973, 513, 395, 303), 2, 2))

############
# propensity score matching
#########

# fit a propensity score model with logistic regression

psmodel <- glm(treatment ~ ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+age+female+meanbp1,
               family = binomial(),
               data = mydata)

# show coeffients
summary(psmodel)

# create propensity score
pscore<-psmodel$fitted.values

logit <- function(p) {log(p)-log(1-p)}
psmatch <-Match(Tr=mydata$treatment, M = 1, X = logit(pscore),
                replace = FALSE, caliper = 0.2)
matched<-mydata[unlist(psmatch[c("index.treated", "index.control")]),]

# get standardized differences
matchedtab1 <- CreateTableOne(vars = xvars, strata = "treatment", data = matched, test = FALSE)

print(matchedtab1, smd = TRUE)

# OUTCOME ANALYSIS
y_trt <-matched$died[matched$treatment==1]
y_con <- matched$died[matched$treatment==0]

# pairwise difference
diffy <- y_trt - y_con

#paired t.test
t.test(diffy)

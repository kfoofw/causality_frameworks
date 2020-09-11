library(ivpack)

# Read dataset
data(card.data)

names(card.data)
# IV is nearc4 (near 4 year college)
# outcome is lwage (log of wage)
# "treatment" is educ (number of years of education)

# summary stats
mean(card.data$nearc4)
par(mfrow = c(1,2))
hist(card.data$lwage)
hist(card.data$educ)

# Is the IV strongly associated with the treatment? strength of IV
mean(card.data$educ[card.data$nearc4 == 1])
mean(card.data$educ[card.data$nearc4 == 0])

# Making educ binary
educ12 <- card.data$educ >12

# Assume no defiers
# Estimate compliers prop: P(A = 1|Z = 1) - P(A = 1|Z = 0) = Population_prop(Always_takers & Compliers) - Population_prop(Always_takers)
propcomp <- mean(educ12[card.data$nearc4==1]) - mean(educ12[card.data$nearc4==0])
propcomp

# Intention to treat effect: P(Y|Z = 1) - P(Y|Z = 0)
itt <- mean(card.data$lwage[card.data$nearc4 == 1]) - mean(card.data$lwage[card.data$nearc4 == 0])
itt

# Complier Average Treatment Effect (CATE): itt / propcomp = [P(Y|Z = 1) - P(Y|Z = 0)] / [P(A = 1|Z = 1) - P(A = 1|Z = 0)]
itt/propcomp

########################
# Two stage least squares
########################

# Stage 1: Regress A on Z
s1 <- lm(educ12 ~ card.data$nearc4)
summary(s1)
## get predicted value of A given Z for each subject
# A_hat
predtx <- predict(s1, type = "response")
table(predtx)

# Stge 2: Regress outcome on A_hat
s2 <- lm(card.data$lwage ~ predtx)
coef(s2)[2]

# 2SLS with ivpack
ivmodel <- ivreg(lwage ~ educ12, # outcome ~ treatment 
                 ~ nearc4, # ~ encouragement
                 x = TRUE, # return X values
                 data = card.data)
# Find robust standard errors of coefficients (instead of using default summary(ivmodel))
robust.se(ivmodel)

# 2SLS with covariates
ivmodel_cov <- ivreg(lwage ~ educ12 + exper + reg661 + reg662 + reg663 + reg664 + reg665+ reg666 + reg667 + reg668,
                     ~ nearc4 + exper + reg661 + reg662 + reg663 + reg664 + reg665+ reg666 + reg667 + reg668,
                     x = TRUE,
                     data = card.data)

robust.se(ivmodel_cov)


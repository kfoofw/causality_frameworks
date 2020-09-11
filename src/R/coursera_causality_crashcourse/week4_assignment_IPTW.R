library(tableone)
library(Matching)
library(ipw)
library(survey)

library(MatchIt)
data(lalonde)

# Fit a propensity score model. 
# Use a logistic regression model, where the outcome is treatment. 
# Include the 8 confounding variables in the model as predictors, with no interaction terms or non-linear terms (such as squared terms). Obtain the propensity score for each subject. Next, obtain the inverse probability of treatment weights for each subject.

xvars <- c("age", "educ", "black", "hispan", "married", "nodegree", "re74", "re75")

names(lalonde)
ps_model <- glm(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75,
                data = lalonde,
                family = binomial())

summary(ps_model)

ps_scores <- ps_model$fitted.values

weights = ifelse(lalonde$treat == 1, 1/ps_scores, 1/(1-ps_scores))

# Question 1
# What are the minimum and maximum weights?
summary(weights)

min_weight = min(weights)
max_weight = max(weights)

ipwplot(weights)

print(paste0("Min weight is ",round(min_weight, 4) ," and Max weight is ", round(max_weight,4)))

# Question 2
# Find the standardized differences for each confounder on the weighted (pseudo) population. What is the standardized difference for nodegree?

# Create pseudo population using weighted data
weighted_lalonde <- svydesign(ids = ~1,
                              data = lalonde,
                              weights =  ~weights)

# Create tableone of pseudo population
weighted_tableone <- svyCreateTableOne(vars = xvars,
                                       strata = "treat",
                                       data = weighted_lalonde,
                                       test = FALSE)

print(weighted_tableone, smd = TRUE)

print(paste0("Standardized difference for 'nodegree' is 0.112."))

# Question 3
# Using IPTW, find the estimate and 95% confidence interval for the average causal effect. This can be obtained from svyglm

# Marginal structural model of pseudo population
msm <- svyglm(re78 ~ treat,
              design = weighted_lalonde)

coef(msm)
confint(msm)

print(paste0("Using IPTW, the estimate of the average causal effect is ", round(coef(msm)[2], 2), " and the 95 % CI is between ", round(confint(msm)[2,1],3) ," and ", round(confint(msm)[2,2],3)))

# Question 4
# Now truncate the weights at the 1st and 99th percentiles. This can be done with the trunc=0.01 option in svyglm.
# Using IPTW with the truncated weights, find the estimate and 95% confidence interval for the average causal effect

# Need to use ipwpoint function for percentile truncation of weights
weight_trunc_model <- ipwpoint(exposure = treat, 
                               family = "binomial", 
                               link = "logit", 
                               denominator = ~ age+educ+black+hispan+married+nodegree+re74+re75, 
                               data = lalonde, 
                               trunc = 0.01)

# plot out distribution of weights after truncation
ipwplot(weight_trunc_model$weights.trunc)

# Create pseudo population with trunc weights
trunc_weighted_lalonde <- svydesign(ids = ~1,
                                    data = lalonde,
                                    weights =  ~ weight_trunc_model$weights.trunc,
                                    trunc = 0.01)

# Create tableone based on trunc weight pseudo population
trunc_weighted_tableone <- svyCreateTableOne(vars = xvars,
                                             strata = "treat",
                                             data = trunc_weighted_lalonde,
                                             test = FALSE)
print(trunc_weighted_tableone, smd = TRUE)

# Marginal structural model of trunc weights pseudo population
msm_trunc <- svyglm(re78 ~ treat,
                    design = trunc_weighted_lalonde)

coef(msm_trunc)
confint(msm_trunc)

print(paste0("Using IPTW with truncated weights at 0.01,the estimate of the average causal effect is ", round(coef(msm_trunc)[2], 2), " and the 95 % CI is between ", round(confint(msm_trunc)[2,1],3) ," and ", round(confint(msm_trunc)[2,2],3)))

library(tableone)
library(Matching)
library(MatchIt)

data(lalonde)

head(lalonde)
names(lalonde)

# treatment is "treat"
# outcome is "re78"
# confounding variables are: age, educ, black, hispan, married, nodegree, re74, re75

# Question 1
# Find the standardized differences for all of the confounding variables (pre-matching). What is the standardized difference for married (to nearest hundredth)?
xvars = c("age", "educ", "black", "hispan", "married", "nodegree", "re74", "re75")
prematched_tableone <- CreateTableOne(vars = xvars, 
                                      strata = "treat", 
                                      data = lalonde,
                                      test = FALSE)

print(prematched_tableone, smd = TRUE)

print(paste0("SMD for 'married is ", 0.719))

# Question 2
# What is the raw (unadjusted) mean of real earnings in 1978 for treated subjects minus the mean of real earnings in 1978 for untreated subjects?
mean_re78_treat = mean(lalonde$re78[lalonde$treat==1])
mean_re78_untreat = mean(lalonde$re78[lalonde$treat==0])

print(paste0("Raw mean of real earnings in 1978 for treated subjects minus untreated subjects is ", mean_re78_treat - mean_re78_untreat))

# Question 3
# Fit a propensity score model. Use a logistic regression model, where the outcome is treatment. Include the 8 confounding variables in the model as predictors, with no interaction terms or non-linear terms (such as squared terms). Obtain the propensity score for each subject.

head(lalonde)

ps_model <- glm(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75,
                data = lalonde,
                family = binomial()) 

summary(ps_model)

ps_scores <- ps_model$fitted.values

summary(ps_scores)
print(paste0("Min value of propensity score is ",round(min(ps_scores), 5)," and Max value of propensity score is ", round(max(ps_scores),5)))

# Question 4
# Use options to specify pair matching, without replacement, no caliper.
# Match on the propensity score itself, not logit of the propensity score. Obtain the standardized differences for the matched data.
# What is the standardized difference for married?
set.seed(931139)

ps_match <- Match(Tr = lalonde$treat, 
                  X = ps_scores, # X is input basis for matching. Since we are matching on propensity scores, it should be ps_scores. If we are matching on covariates, we can put 'lalonde[xvars]'
                  replace = FALSE)

# create subset of original data based on matching of propensity scores
matched_lalonde <- lalonde[unlist(ps_match[c("index.treated", "index.control")]),]

matched_tableone <- CreateTableOne(vars = xvars,
                                   strata = "treat",
                                   data = matched_lalonde,
                                   test = FALSE
                                   )
print(matched_tableone, smd = TRUE)

# Question 5
# For the propensity score matched data:
# Which variable has the largest standardized difference?

print(paste0("For the propensity scored matched data, the variable with the largest standardized difference is 'black' with SMD value of 0.852."))

# Question 6
# Re-do the matching, but use a caliper this time. Set the caliper=0.1 in the options in the Match function.
# Again, before running the Match function, set the seed
# How many matched pairs are there?
  
set.seed(931139)
ps_match_caliper <- Match(Tr = lalonde$treat, 
                          X = ps_scores, # X is input basis for matching. Since we are matching on propensity scores, it should be ps_scores. If we are matching on covariates, we can put 'lalonde[xvars]'
                          replace = FALSE,
                          caliper = 0.1) # caliper is the distance which is acceptable for any match

matched_caliper_lalonde <- lalonde[unlist(c(ps_match_caliper$index.treated, ps_match_caliper$index.control)),]

matched_caliper_tableone <- CreateTableOne(vars = xvars, strata = "treat", data = matched_caliper_lalonde, test = FALSE)

print(matched_caliper_tableone, smd = TRUE)

print(paste0("Number of matched pairs are 111"))

# Question 7
# Use the matched data set (from propensity score matching with caliper=0.1) to carry out the outcome analysis.
# For the matched data, what is the mean of real earnings in 1978 for treated subjects minus the mean of real earnings in 1978 for untreated subjects?
mean_ps_re78_treat <- mean(matched_caliper_lalonde$re78[matched_caliper_lalonde$treat == 1])
mean_ps_re78_untreat <- mean(matched_caliper_lalonde$re78[matched_caliper_lalonde$treat == 0])

print(paste0("For matched data, mean of real earnings in 1978 for treated subjects minus untreated subjects is ", mean_ps_re78_treat - mean_ps_re78_untreat))

ps_re78_treat <- matched_caliper_lalonde$re78[matched_caliper_lalonde$treat == 1]
ps_re78_untreat <- matched_caliper_lalonde$re78[matched_caliper_lalonde$treat == 0]

# Obtain vector of pairwise difference
diffy <- ps_re78_treat - ps_re78_untreat
head(diffy)

t.test(diffy)

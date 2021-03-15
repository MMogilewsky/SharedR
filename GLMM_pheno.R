############## GLMM for Phenology Analysis #####################

### Preparation for setting up model######################################

### Import dataset with all phenological data, each tree represented once per monitoring visit (biweekly for transect)
library(readxl)


AllTrees_Pheno_2017 <- read_excel("C:/Users/msmog/Documents/GitHub/SharedR/AllTrees_Pheno_NoBehaviorRepeats_2017.xlsx")
View(AllTrees_Pheno_2017)
AllTrees <- drop_na(AllTrees_Pheno_2017, T_New_Lvs, T_Fruit) # data set with all trees, rows missing values in T_New_Lvs & T_Fruit columns removed
describeBy(AllTrees)
AllTrees$T_New_Lvs.t <- AllTrees$T_New_Lvs + 1 # added term to enable handling zero data
AllTrees$T_Fruit.t <- AllTrees$T_Fruit + 1 # added term to enable handling zero data
colnames(AllTrees)
AllTrees_s <- AllTrees %>% mutate_at(c(12, 16, 20, 21), funs(c(scale(.)))) ### after running several models, determined that continuous variables needed to be scaled so that eigenvalues weren't so large.
LianaTrees <- drop_na(AllTrees, L_New_Lvs, L_Fruit) # data set with only trees with lianas
LianaTrees$L_New_Lvs.t <- LianaTrees$L_New_Lvs + 1
LianaTrees$L_Fruit.t <- LianaTrees$L_Fruit + 1
LianaTrees_s <- LianaTrees %>% mutate_at(c(12, 16, 20, 21), funs(c(scale(.))))

### Create data set with just transect trees
trans <- filter(AllTrees_Pheno_2017, Species == "Transect")
View(trans)

### Import lemur feeding trees data set

library(readxl)
Feeding <- read_excel("C:/Users/msmog/Documents/GitHub/SharedR/AllTransect_OnlyLemurFeeding_Pheno_NoBehaviorRepeats_2017.xlsx", 
                      col_types = c("text", "text", "text", 
                                    "text", "numeric", "text", "numeric", 
                                    "numeric", "numeric", "date", "text", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric"))
View(Feeding)
feeding <- filter(Feeding, Species != "Transect")
View(feeding)
### Libraries used
library(readxl)
# run flattenCorMatrix.R
library(PerformanceAnalytics)
library(Hmisc)
library(lattice)
library(tidyr)
library(dplyr)
library(ggplot2)
library(MASS)
library(car)
library(fitdistrplus)
library(psych)
library(lme4)
#run overdispersion_function.R


### Determining correlation between response variables
colnames(trans)
resp <- trans[, c(23:26,30:33)]

cor(resp, use = "complete.obs", method = "spearman") ## correlation coefficients using rows with complete data
respres <- rcorr(as.matrix(resp)) ## p-values

flattenCorrMatrix(respres$r, respres$P)
chart.Correlation(resp, histogram = TRUE, method = "spearman", pch=19)

colnames(feeding)
resp2 <- feeding[, c(19:22,26:29)]

cor(resp2, use = "complete.obs", method = "spearman") ## correlation coefficients using rows with complete data
resp2res <- rcorr(as.matrix(resp2)) ## p-values

flattenCorrMatrix(resp2res$r, resp2res$P)
chart.Correlation(resp2, histogram = TRUE, method = "spearman", pch=19)
## None of the response variables have a correlation ~ 0.8 for either dataset, the value Dr. Pan recommended as a threshold for dropping.  The highest correlations are between buds and new leaves of their respective trees or lianas (which makes sense as buds were probably leaf buds that would soon become new leaves).  Results mean that I should run 4 models (ugh!), one for each of the following response variables:  tree young leaves, tree fruit, liana young leaves, liana fruit

### Using plots to select best probability distribution, since it is definitely not normal/Guassian (see TestingAssumptionsofLinearModels.R); have to do this step for each response variable that I want to model
## Transect Trees first
trans <- drop_na(trans, T_New_Lvs, T_Fruit) # Creating a dataset with no missing data in the tree response variables so that distributions that can't handle missing data will run
dim(trans)

# Tree Young Leaves

trans$T_New_Lvs.t <- trans$T_New_Lvs +1 # This is necessary because several distributions can't handle zeros
qqp(AllTrees_Pheno_2017$T_New_Lvs.t, "norm") # plot normal distribution first
qqp(AllTrees_Pheno_2017$T_New_Lvs.t, "lnorm") # log normal distribution is better but still bad, note several extreme outliers
# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter 

nbinom <- fitdistr(trans$T_New_Lvs.t, "Negative Binomial")
qqp(trans$T_New_Lvs.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(trans$T_New_Lvs.t, "Poisson")
qqp(AllTrees_Pheno_2017$T_New_Lvs.t, "pois", poisson$estimate, lambda = 0.5)
?fitdistr

gamma <- fitdistr(trans$T_New_Lvs.t, "gamma")
qqp(trans$T_New_Lvs.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

## beta distribution is recommended for % canopy cover and other % data (data bounded between 0-1) so first converting my % data to 0-1 by dividing by 100
trans$T_New_Lvs.d <- trans$T_New_Lvs/100 
beta <- fitdist(trans$T_New_Lvs.d, "beta", method = "mme") ## using library package fitdistrplus because the function fitdist will provide reasonable start list whereas I must provide a start list for function fitdistr in MASS, uses moment matching estimate
plot(beta, las = 1)



##### Function for selecting parameters for beta distribution if defaults don't work; defaults worked okay for me but here's the code just in case ##################################################
beta_mom <- function(x) {
  
  m_x <- mean(x, na.rm = TRUE)
  s_x <- sd(x, na.rm = TRUE)
  
  alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
  beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)
  
  return(list(alpha = alpha, beta = beta))
  
}

beta_mom(AllTrees_Pheno_2017$T_New_Lvs.d) ### provides a start 

#######################################################################

#### None of these distributions are really good, but lnorm and beta appears to be the absolute closest to something decent.  Beta is specifically for percent data (data bounded between 0-1) while lnorm is more like a basic log transformation.  Maybe try the lnorm first and then discuss the options with Andres?

######################## Choosing a method for estimating effect sizes
describeBy(AllTrees_Pheno_2017, group = "Species") # Need the mean for response variables because PQL (one of the methods for estimating effect sizes) gets squirrely when means are <5.  In this case they are and so I should avoid PQL.

### Laplace works well with <3 random effects and non-normally distributed data that violates PQL requirements; REML works well with non-normal distributions that can be transformed.  Not sure which to use in my case given that the lnorm and beta probability distributions were about equivalent.

### Had John look at the QQ plot of log distribution vs QQ plot of beta distribution (without him knowing which was which) and he thought that the log was the best fit, so for now, I'm going to proceed with a log normal distribution.

###################################################################################################
glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)) ## after receiving lots of errors about failure to converge, internet recommendations were to use a different optimizer - "bobyqa" and increase the iterations 100,000.

############ Tree New Leaves
lm <- lm(T_New_Lvs.t ~ Distance_to_Gap, data = AllTrees_s)
summary(lm)

Laplace <- glmer(T_New_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID), data = AllTrees_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace) # overdispersed
summary(Laplace) # AIC = 41565.5  BIC = 41591.8, errors include failure to converge and large eigenvalue

Laplace.A <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + (1 | Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.A) # overdispersed
summary(Laplace.A) # AIC = 41530.8  BIC = 41570.3, no errors

Laplace.B <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.B) # overdispersed
summary(Laplace.B) # AIC = 41300.6  BIC = 41346.6, no errors

Laplace.C <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.C) # overdispersed
summary(Laplace.C) # AIC = 38727.8  BIC = 38773.9, better than just Tree_ID in residuals

Laplace.D <- glmer(T_New_Lvs.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.D) # still overdispersed?
summary(Laplace.D) # AIC = 34622.7  BIC = 34655.6, no errors

Laplace.E <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.E) # still overdispersed?
summary(Laplace.E) # AIC = 34441.7  BIC = 34481.1, no errors, not much improvement on model Laplace.D

Laplace.F <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.F) # finally fixed the overdispersion problem!
summary(Laplace.F) # AIC = 30822.1  BIC = 30868.2, best model so far, no errors

Laplace.G <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.G) # not overdispersed
summary(Laplace.G) # AIC= 29905.0  BIC = 29957.6, only a slight improvement over Laplace.F, no errors, BEST MODEL

##### Plot to see how well model fits data
plot(fitted(Laplace.G), residuals(Laplace.G), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Laplace.G), residuals(Laplace.G))) ### plot of solid line and dashed line have some definite weirdnesses but there are large areas of overlap, which is good.

anova(Laplace.G, Laplace.F, test="Chisq") # comparing models that were not overdispersed or had errors, but I get error message, saying that models were not fitted to the same size of dataset but they were all based on the same dataset and I didn't get any weird dropped errors for any of these models?

Laplace.H <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Gap_Area + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.H) # not overdispersed
summary(Laplace.H) # AIC = 30596.3  BIC = 30648.9, no errors, worse than G but better than F

Laplace.I <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Gap_Area + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.I) # not overdispersed
summary(Laplace.I) # AIC = 29616.8  BIC = 29675.9, model failed to converge but better than G

Laplace.J <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Gap_Area + Tree_DBH + Distance_to_Gap*Gap_Area + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.J) # not overdispersed
summary(Laplace.J) # lots of strange errors and everything significant, AIC = 29633.3  BIC = 29699.0 but overall, I don't trust these results

Laplace.K <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Gap_Area + Season + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.K) # not overdispersed
summary(Laplace.K) # AIC = 30551.4  BIC = 30610.5, did not converge and not as good as G.

Laplace.L <- glmer(T_New_Lvs.t ~ Species + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.L) # not overdispersed
summary(Laplace.L) # AIC= 30050.8  BIC = 30096.9, not as good as G

Laplace.M <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Sex + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.M) # not overdispersed
summary(Laplace.M) # AIC = 29665.4  BIC = 29731.2, but got  a lot of error messages, including failure to converge, slightly lower AIC value than the other better model G

Laplace.N <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Species + Species*Sex + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.N) # not overdispersed
summary(Laplace.N) #  29573.6  BIC = 29645.9, but lots of errors, including dropped columns


####################### Tree Fruit
lm <- lm(T_Fruit.t ~ Distance_to_Gap, data = AllTrees_s)
summary(lm)

Laplace <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace) # overdispersed
summary(Laplace) # AIC = 35668.1  BIC = 35694.4, no errors

Laplace.A <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + (1 | Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.A) # overdispersed
summary(Laplace.A) # AIC =  35654.7  BIC = 35694.2, no errors, very similar to Laplace

Laplace.B <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.B) # overdispersed
summary(Laplace.B) # AIC = 35447.0  BIC = 35493.1, no errors, not very different from Laplace

Laplace.C <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.C) # overdispersed
summary(Laplace.C) # AIC = 35435.7  BIC = 35481.8, failure to converge, very large eigenvalue, not very different from Laplace

Laplace.D <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.D) # still overdispersed?
summary(Laplace.D) # AIC = 27279.7  BIC = 27312.6, no errors, much better than other models

Laplace.E <- glmer(T_Fruit.t ~ Distance_to_Gap + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.E) # still overdispersed
summary(Laplace.E) # AIC = 27051.3  BIC = 27090.7, no errors, slightly better than Laplace.D

Laplace.F <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.F) # finally fixed the overdispersion problem!
summary(Laplace.F) # AIC = 26779.2  BIC = 26825.3, no errors, only slightly better than Lacplace.E but fixed overdispersion problem, overall BEST MODEL

Laplace.G <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.G) # overdispersed
summary(Laplace.G) # AIC= 26514.9  BIC = 26567.5, no errors, but overdispersed and not much better than Laplace.F

Laplace.H <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Season + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.H) # overdispersed so summary not reliable
summary(Laplace.H)

Laplace.I <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Sex +(1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.I) # not overdispersed
summary(Laplace.I) # AIC = 25228.7  BIC = 25288.0, but had a dropped column?  so not sure how reliable the summary is

Laplace.J <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Sex + Species*Sex + (1 | Week/Tree_ID), data = AllTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.J) # not overdispersed
summary(Laplace.J) # AIC = 24996.1  BIC = 25061.9, best AIC model but 6 dropped columns?  so not sure how reliable the summary is

#################### Young Liana Leaves
lm <- lm(L_New_Lvs.t ~ Distance_to_Gap, data = LianaTrees_s)
summary(lm)

Laplace <- glmer(L_New_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID), data = LianaTrees_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace) # overdispersed
summary(Laplace) # AIC = 29524.8  BIC = 29549.1, errors include failure to converge and large eigenvalue

Laplace.A <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + (1 | Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.A) # overdispersed
summary(Laplace.A) # AIC = 29521.1  29557.6, same scores as above, errors include failure to converge and large eigenvalue

Laplace.B <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.B) # overdispersed
summary(Laplace.B) # AIC = 29232.5  BIC = 29275.1, no errors

Laplace.C <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.C) # overdispersed
summary(Laplace.C) # AIC = 25823.3  BIC = 25865.9, no errors, better than just Tree_ID in residuals

Laplace.D <- glmer(L_New_Lvs.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.D) # not overdispersed
summary(Laplace.D) # AIC = 23077.4  BIC = 23107.8, errors include failure to converge and large eigenvalue, better than Laplace.C but very problematic

Laplace.E <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Tree_DBH + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.E) # not overdispersed
summary(Laplace.E) # AIC = 22807.1  BIC = 22843.6, no errors, a little better than Laplace.D and much better than earlier models, BEST MODELS

Laplace.F <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.F) # not overdispersed
summary(Laplace.F) # AIC = 22970.2  23012.8, slightly worse than Laplace.E, no errors

Laplace.G <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.G) # not overdispersed
summary(Laplace.G) # AIC =  22740.6  BIC = 22789.2, no errors, very similar to Laplace.E

Laplace.H <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + Season + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.H) # not overdispersed
summary(Laplace.H) # AIC = 22334.2  BIC = 22388.9, no errors, BEST MODEL

Laplace.I <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + Season + Sex + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.I) # not overdispersed
summary(Laplace.I) # AIC = 22302.5  BIC = 22369.4, failed to converge and dropped a column, not really better than H

Laplace.J <- glmer(L_New_Lvs.t ~ Distance_to_Gap + Species + Tree_DBH + Sex + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.J) # not overdispersed
summary(Laplace.J) # AIC = 22532.8  BIC = 22593.5, failure to converge, dropped a column, and not as good as H

################### Liana fruit
lm <- lm(L_Fruit.t ~ Distance_to_Gap, data = LianaTrees_s)
summary(lm)

Laplace <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Tree_ID), data = LianaTrees_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace) # overdispersed
summary(Laplace) # AIC = 16867.6  BIC = 16891.9, no errors

Laplace.A <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + (1 | Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.A) # overdispersed
summary(Laplace.A) # AIC = 16452.3  BIC = 16488.8, no errors, slightly better than Laplace.A

Laplace.B <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.B) # overdispersed
summary(Laplace.B) # AIC = 16234.9  BIC = 16277.5, no errors, slightly better than B

Laplace.C <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.C) # overdispersed
summary(Laplace.C) # AIC = 17935.5 BIC = 17978.1, no errors, worse than A & B

Laplace.D <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.D) # not overdispersed
summary(Laplace.D) # AIC = 10048.8  BIC = 10079.2, no errors, much better than all other models

Laplace.E <- glmer(L_Fruit.t ~ Distance_to_Gap + Tree_DBH + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.E) # not overdispersed
summary(Laplace.E) # AIC = 9918.5   BIC = 9955.0, no errors, slightly better than D

Laplace.F <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.F) # not overdispersed
summary(Laplace.F) # AIC = 8369.0   BIC = 8411.6, no errors, much better than D & E

Laplace.G <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.G) # not overdispersed
summary(Laplace.G) # AIC =  8295.5   BIC = 8344.1, no errors, not really much better than F

Laplace.H <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + Season + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.H) # not overdispersed
summary(Laplace.H) # AIC = 7681.6   7736.3, no errors, BEST MODEL

Laplace.I <- glmer(L_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + Sex + (1 | Week/Tree_ID), data = LianaTrees_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.I) # not overdispersed
summary(Laplace.I) # AIC = 3185.7   BIC = 3246.5, failure to converge, one dropped column

#############################################################################################

######## Lemur Feeding Trees Only##########################################################################
colnames(feeding)
feeding_s <- feeding %>% mutate_at(c(12, 15, 16, 17), funs(c(scale(.)))) ### after running several models, determined that continuous variables needed to be scaled so that eigenvalues weren't so large.
str(feeding_s)
feeding_s$Species <- as.factor(feeding_s$Species)
feeding_s$Sex <- as.factor(feeding_s$Sex)
feeding_s <- drop_na(feeding_s, T_Fruit, T_New_Lvs)
summary(feeding_s)
describeBy(feeding, group = "Species")


################################## Tree fruit
feeding$T_Fruit.t <- feeding$T_Fruit + 1
feeding_s$T_Fruit.t <- feeding_s$T_Fruit + 1 # This is necessary because several distributions can't handle zeros
qqp(feeding$T_Fruit.t, "norm") # plot normal distribution first
qqp(feeding$T_Fruit.t, "lnorm") # log normal probability distribution isn't that bad, actually

## beta distribution is recommended for % canopy cover and other % data (data bounded between 0-1) so first converting my % data to 0-1 by dividing by 100
feeding$T_Fruit.d <- feeding$T_Fruit/100 
beta <- fitdist(feeding$T_Fruit.d, "beta", method = "mme") ## using library package fitdistrplus because the function fitdist will provide reasonable start list whereas I must provide a start list for function fitdistr in MASS, uses moment matching estimate
plot(beta, las = 1)
### log normal and beta distribution are similar

#####################



############## Tree Fruit

Laplace <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Tree_ID), data = feeding_s, family = gaussian(link = "log")) # simplest model
overdisp_fun(Laplace) # overdispersed
summary(Laplace) # AIC = 6355.9   BIC = 6373.8, no errors

Laplace.A <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + (1 | Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.A) # overdispersed
summary(Laplace.A) # AIC = 6350.7   6373.1 , no errors, very similar to Laplace

Laplace.B <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.B) # overdispersed
summary(Laplace.B) # AIC = 6293.2   BIC = 6320.0, no errors, slightly better than Laplace

Laplace.C <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.C) # overdispersed
summary(Laplace.C) # AIC = 5200.1   BIC = 5226.9, failure to converge, very large eigenvalue

Laplace.D <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.D) # still overdispersed
summary(Laplace.D) # AIC = 5919.4  BIC = 5941.7, no errors, worse than C

Laplace.E <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Session), data = feeding_s, family = gaussian(link = "log")) # trying to addres the overdispersion problem
overdisp_fun(Laplace.E) # still overdispersed
summary(Laplace.E) # AIC = 5607.7   BIC = 5625.6, large eigenvalue, failure to converge, better than C, worse than D

View(feeding_s)
feeding_s$Session <- as.factor(feeding_s$Session)
Laplace.F <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) # still trying to correct the overdispersion problem
overdisp_fun(Laplace.F) # finally fixed the overdispersion problem!
summary(Laplace.F) # AIC = 4008.6   BIC = 4039.9, no errors, best model yet

Laplace.G <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.G) # not overdispersed
summary(Laplace.G) # AIC = 3966.9  BIC = 4002.6, no errors, only slightly better than F


Laplace.H <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + Gap_Area + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.H) # not overdispersed
summary(Laplace.H) # AIC = 3787.8   BIC = 3827.5, no errors

Laplace.I <- glmer(T_Fruit.t ~ Distance_to_Gap + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.I) # not overdispersed
summary(Laplace.I) # AIC = 3965.5   BIC = 3996.7, no errors, worse than H

Laplace.J <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.J) # not overdispersed
summary(Laplace.J) # AIC = 3785.8  BIC = 3821.1, no errors, same as H but slightly simpler, BEST MODEL

Laplace.K <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Distance_to_Gap*Gap_Area + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.K) # not overdispersed
summary(Laplace.K) # AIC = 3787.4   BIC = 3827.1, no errors, same as J but more complex

Laplace.L <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Distance_to_Gap*Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.L) # not overdispersed
summary(Laplace.L) # AIC = 3785.0   BIC = 3824.8, no errors, same as J but more complex

Laplace.M <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Season + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.M) # not overdispersed
summary(Laplace.M) # AIC = 3784.7   3824.4, similar to J but more complex

Laplace.N <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Sex + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.N) # not overdispersed
summary(Laplace.N) # AIC = 3789.7   BIC = 3833.8, similar to J but more complex

###################### Young feeding tree leaves
feeding$T_New_lvs.t <- feeding$T_New_Lvs + 1
feeding_s$T_New_lvs.t <- feeding_s$T_New_Lvs + 1
Laplace <- glmer(T_New_lvs.t ~ Distance_to_Gap + (1 | Tree_ID), data = feeding_s, family = gaussian(link = "log")) # simplest model
overdisp_fun(Laplace) # overdispersed
summary(Laplace) # AIC = 5142.1   BIC = 5160.0, no errors

Laplace.A <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + (1 | Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.A) # overdispersed
summary(Laplace.A) # AIC = 5142.0   5164.3 , no errors, very similar to Laplace

Laplace.B <- glmer(T_Fruit.t ~ Distance_to_Gap + Tree_DBH + (1 | Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.B) # overdispersed
summary(Laplace.B) # AIC = 6299.4   BCI = 6321.7, no errors, much worse than Laplace

Laplace.C <- glmer(T_New_lvs.t ~ Distance_to_Gap + (1 | Week), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.C) # overdispersed
summary(Laplace.C) # AIC = 4792.3   4810.2, no errors, better than all previous models

Laplace.D <- glmer(T_New_lvs.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = feeding_s, family = gaussian(link = "log")) # trying to address overdispersion
overdisp_fun(Laplace.D) # not overdispersed!
summary(Laplace.D) # AIC = 3234.1   3256.5, no errors, much, much better model

Laplace.E <- glmer(T_New_lvs.t ~ Distance_to_Gap + (1 | Session), data = feeding_s, family = gaussian(link = "log")) # trying to addres the overdispersion problem
overdisp_fun(Laplace.E) # overdispersed
summary(Laplace.E) # AIC = 5068.7   BIC = 5086.6, much worse model

Laplace.F <- glmer(T_New_lvs.t ~ Distance_to_Gap + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.F) # not overdispersed!
summary(Laplace.F) # AIC = 1375.4 BIC = 1402.2, errors abound! but lowest AIC score

Laplace.G <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.G) # not overdispersed!
summary(Laplace.G) # AIC = 1376.6   BIC = 1407.9, similar score to F but seemed to address some of the errors, BEST MODEL

Laplace.H <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + Gap_Area + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.H) # not overdispersed!
summary(Laplace.H) # AIC = 1358.0   BIC = 1393.4, similar score to G but more complex

Laplace.I <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.I) # not overdispersed!
summary(Laplace.I) # AIC = 1374.9   BIC = 1410.7, similar score to G but more complex

Laplace.J <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + Distance_to_Gap*Species + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.J) # not overdispersed!
summary(Laplace.J) # AIC = 1377.7   BIC = 1413.5, similar to G but more complex

Laplace.K <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + Sex + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.K) # not overdispersed!
summary(Laplace.K) # AIC = 1380.5   BIC = 1420.8, slightly worse than model G

Laplace.L <- glmer(T_New_lvs.t ~ Distance_to_Gap + Species + Season + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.L) # not overdispersed!
summary(Laplace.L) # AIC = 1352.2   BIC = 1388.0, same AIC as G but more complex

#####################################

View(feeding_s)
feeding_s$Session <- as.factor(feeding_s$Session)
Laplace.F <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log")) # still trying to correct the overdispersion problem
overdisp_fun(Laplace.F) # finally fixed the overdispersion problem!
summary(Laplace.F) # AIC = 4008.6   BIC = 4039.9, no errors, best model yet

Laplace.G <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.G) # not overdispersed
summary(Laplace.G) # AIC = 3966.9  BIC = 4002.6, no errors, only slightly better than F


Laplace.H <- glmer(T_Fruit.t ~ Distance_to_Gap + Species + Tree_DBH + Gap_Area + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.H) # not overdispersed
summary(Laplace.H) # AIC = 3787.8   BIC = 3827.5, no errors

Laplace.I <- glmer(T_Fruit.t ~ Distance_to_Gap + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.I) # not overdispersed
summary(Laplace.I) # AIC = 3965.5   BIC = 3996.7, no errors, worse than H

Laplace.J <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.J) # not overdispersed
summary(Laplace.J) # AIC = 3785.8  BIC = 3821.1, no errors, same as H but slightly simpler, BEST MODEL

Laplace.K <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Distance_to_Gap*Gap_Area + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.K) # not overdispersed
summary(Laplace.K) # AIC = 3787.4   BIC = 3827.1, no errors, same as J but more complex

Laplace.L <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Distance_to_Gap*Tree_DBH + (1 | Week/Session/Tree_ID), data = feeding_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.L) # not overdispersed
summary(Laplace.L) # AIC = 3785.0   BIC = 3824.8, no errors, same as J but more complex


############## Young liana leaves
feedliana_s <- drop_na(feeding_s, L_Fruit, L_New_Lvs)
feedliana_s$L_Fruit.t <- feedliana_s$L_Fruit + 1
feedliana_s$L_Lvs.t <- feedliana_s$L_New_Lvs + 1
feedliana_s$Session <-as.factor(feedliana_s$Session)

lia.1 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(lia.1) # overdispersed

lia.2 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Session), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(lia.2) # overdispersed

lia.3 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Week), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(lia.3) # overdispersed, failure to converge, large eigenvalues despite scaling

lia.4 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Session), data = feedliana_s, family = gaussian(link = "log")) 
# model failed to run
glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000))

lia.5 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Week), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(lia.5) # overdispersed 
summary(lia.5) # AII = 3421.6, no errors but overdispersed

lia.6 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Season), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(lia.6) # still overdispersed

lia.7 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Week/Session), data = feedliana_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.7) # model failed to run


### Not able to correct overdispersion problem.  Maybe just stick with the data set that includes transect trees and all lemur trees (rather than just feeding trees) for now

#### Liana fruit
Laplace.1 <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Tree_ID), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace) # overdispersed
summary(Laplace.1)

Laplace.2 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Session), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.2) # overdispersed

Laplace.3 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Week), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.3) # overdispersed, failure to converge, large eigenvalues despite scaling

Laplace.4 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Session), data = feedliana_s, family = gaussian(link = "log")) 
# model failed to run

Laplace.5 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Week), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.5) # overdispersed 

Laplace.6 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Season), data = feedliana_s, family = gaussian(link = "log")) 
overdisp_fun(Laplace.6) # still overdispersed, not sure what to do

Laplace.7 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Season/Session), data = feedliana_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.7) # model failed to run

Laplace.8 <- glmer(L_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Week/Session), data = feedliana_s, family = gaussian(link = "log"))
overdisp_fun(Laplace.8) # model failed to run
#####################################################################################

#### These are just lemur trees - no transect trees unless visited by a lemur
### inclusion of feeding vs non-feeding category
library(readxl)
FeedingvsNonFeeding <- read_excel("C:/Users/msmog/Dropbox/Dissertation/Data/FeedingvsNonFeeding.xlsx")
View(FeedingvsNonFeeding)
str(FeedingvsNonFeeding)
FeedingvsNonFeeding$Activity <- as.factor(FeedingvsNonFeeding$Activity)
FeedingvsNonFeeding$Species <- as.factor(FeedingvsNonFeeding$Species)
FeedingvsNonFeeding$Season <- as.factor(FeedingvsNonFeeding$Season)
FeedingvsNonFeeding$Sex <- as.factor(FeedingvsNonFeeding$Sex)
FvN <- drop_na(FeedingvsNonFeeding, T_New_Lvs, T_Fruit)
FvN <- filter(FvN, Species != "Transect")
FvN$T_New_Lvs.t <- FvN$T_New_Lvs + 1
FvN$T_Fruit.t <- FvN$T_Fruit + 1
colnames(FvN)
FvN_s <- FvN %>% mutate_at(c(14, 18, 22, 23), funs(c(scale(.)))) # scaling numerical data (distance to gap, gap area, tree height, tree DBH) 
glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)) # increased number of itirations because having trouble converging

FvN_s$T_New_Lvs.t <- FvN_s$T_New_Lvs + 1
FvN_s$T_Fruit.t <- FvN_s$T_Fruit + 1

### Finding best distribution
qqp(FvN$T_New_Lvs.t, "norm") # plot normal distribution first
qqp(FvN$T_New_Lvs.t, "lnorm") # log normal distribution is better but still bad, note several extreme outliers
# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter 

FvN$T_New_Lvs.d <- FvN$T_New_Lvs/100
FvN$T_Fruit.d <- FvN$T_Fruit/100
poisson <- fitdistr(FvN$T_Fruit.d, "Poisson")
qqp(FvN$T_Fruit.d, "pois", poisson$estimate, lambda = 0.5)
?fitdistr

gamma <- fitdistr(FvN$T_New_Lvs.t, "gamma")
qqp(FvN$T_New_Lvs.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
### All distributions poor fit but I guess that log is the best

##### Young Leaves model for all lemur trees
fmod <- glmer(T_New_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod) # overdispersed, no other errors

fmod.1 <- fmod <- glmer(T_New_Lvs.t ~ Distance_to_Gap + (1 | Tree_ID/Week), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.1) # overdispersed
summary(fmod) # overdispersed, no other errors with this one

FvN_s$Session <- as.factor(FvN_s$Session)
fmod.2 <- fmod <- glmer(T_New_Lvs.t ~ Distance_to_Gap + (1 | Week/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.2) # overdispersed
summary(fmod.2) 

fmod.3 <- glmer(T_New_Lvs.t ~ Distance_to_Gap +  (1 | Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.3) # overdispersed
summary(fmod.3)

fmod.3 <- glmer(T_New_Lvs.t ~ Distance_to_Gap +  (1 | Week), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.3) # overdispersed

fmod.4 <- glmer(T_New_Lvs.t ~ Distance_to_Gap +  (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.4) # not overdispersed
summary(fmod.4) # AIC = 10422.5

fmod.5 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Gap_Area + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.5) # not overdispersed
summary(fmod.5)  # 10152.5, slightly better model

fmod.6 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.6) # not overdispersed
summary(fmod.6) # 4909.6, much, much better model

fmod.7 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Species + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.7) # not overdispersed
summary(fmod.7) #4865.6, BEST MODEL

fmod.8 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Species + Sex + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.8) # not overdispersed
#model wouldn't run'

fmod.9 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Species + Distance_to_Gap*Gap_Area + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.9) # not overdispersed
# failed to run

fmod.10 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Tree_DBH + Species + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.10) # not overdispersed
summary(fmod.10) # AIC = 4965.8, only slightly worse model

fmod.11 <- glmer(T_New_Lvs.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Species + Activity + (1 | Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.11) # not overdispersed
summary(fmod.11) # AIC = 4867.6, but scaling problems

### Tree fruit
fmod.A <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Tree_ID), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.A) # overdispersed

fmod.B <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Week), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.B) # overdispersed

fmod.C <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Week/Tree_ID), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.C) # overdispersed, scaling and convergence errors

fmod.D <- glmer(T_Fruit.t ~ Distance_to_Gap + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.D) # not overdispersed
summary(fmod.D) # AIC = 7388.4, no errors

fmod.E <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.E) # not overdispersed
summary(fmod.E) # AIC = 7150.6, slightly improved model, no errors

fmod.F <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.F) # not overdispersed
summary(fmod.F) # AIC = 7079.6, BEST MODEL

fmod.G <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Distance_to_Gap*Gap_Area + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.G) # not overdispersed
summary(fmod.G) # AIC = 7082.0, scaling problems, failed to converge & not better than F.

fmod.H <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Distance_to_Gap*Tree_DBH + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.H) # not overdispersed
summary(fmod.H) # not an improvement, failed to converge, scaling problems

fmod.I <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Season + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.I) # not overdispersed
summary(fmod.I) #  worse model, failed to converge, scaling problems

fmod.J <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Species + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.J) # not overdispersed
summary(fmod.J) # AIC = 7081.7, slightly worse model, scaling problem, failure to converge error

fmod.K <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Sex + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.K) # not overdispersed
summary(fmod.K) # AIC = 7087.7, not a better model, failed to converge

fmod.L <- glmer(T_Fruit.t ~ Distance_to_Gap + Gap_Area + Tree_DBH + Sex*Species + (1 | Week/Tree_ID/Session), data = FvN_s, family = gaussian(link = "log"))
overdisp_fun(fmod.L) # not overdispersed
summary(fmod.L) # 6 columns dropped, failed to converge, not a better model
#### Lemur Lianas
FvNliana_s <- drop_na(FvN_s, L_Fruit, L_New_Lvs)
FvNliana_s$L_Fruit.t <- FvNliana_s$L_Fruit + 1
FvNliana_s$L_Lvs.t <- FvNliana_s$L_New_Lvs + 1
FvNliana_s$Session <-as.factor(FvNliana_s$Session)
dim(FvNliana_s)
### Young Liana Leaves
lmod.A <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Tree_ID), data = FvNliana_s, family = gaussian(link = "log"))
overdisp_fun(lmod.A) # overdispersed

lmod.B <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Week), data = FvNliana_s, family = gaussian(link = "log"))
overdisp_fun(lmod.B) # overdispersed

lmod.C <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Session), data = FvNliana_s, family = gaussian(link = "log"))
overdisp_fun(lmod.C) # overdispersed

lmod.D <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Week/Session), data = FvNliana_s, family = gaussian(link = "log"))
overdisp_fun(lmod.D) # overdispersed

lmod.E <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Week/Session/Tree_ID), data = FvNliana_s, family = gaussian(link = "log"))
## could not get this model to run

lmod.F <- glmer(L_Fruit.t ~ Distance_to_Gap + (1 | Session/Tree_ID), data = FvNliana_s, family = gaussian(link = "log"))
overdisp_fun(lmod.F) # no longer overdispersed, might want to try an alternate optimizer to bobyqa, failed to converge
summary(lmod.F) # AIC = -41827.5 (didn't know you could have negative AIC, ?)

lmod.G <- glmer(L_Fruit.t ~ Distance_to_Gap + Tree_DBH + (1 | Session/Tree_ID), data = FvNliana_s, family = gaussian(link = "log"))
# Would not run

lmod.H <- glmer(L_Fruit.t ~ Distance_to_Gap + Gap_Area + (1 | Session/Tree_ID), data = FvNliana_s, family = gaussian(link = "log"))
overdisp_fun(lmod.H) # not overdispersed but failed to converge, has scaling problems, something seems up with the optimizer
summary(lmod.H) # probably not an accurate model or AIC = -41079.7





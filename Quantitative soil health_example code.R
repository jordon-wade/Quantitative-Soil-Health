library(psych)
library(nFactors)
library(lavaan)
library(MVN)
library(GPArotation)

set.seed(12345)
setwd("~/Dropbox")

#### LOAD IN DATA ####
## North Carolina dataset ##
NC.matrix <- as.matrix(read.csv("Active Projects/Factor analysis and soil health/Data for GitHub/NC data.csv", row.names=1))

## New York-1 dataset ##
NY1.raw <- as.matrix(read.csv("Active Projects/Factor analysis and soil health/Data for GitHub/NY data - Stafford lfs.csv"))
NY1 <- NY1.raw[,1]
NY1.tmp <- NY1.raw[,-1]
NY1.mat <- apply(NY1.tmp, 2, as.numeric)
NY1.mat <- as.matrix(NY1.mat)
rownames(NY1.mat) <- NY1
NY1.mat <- round(NY1.mat, 2)

## New York-2 dataset ##
NY2.raw <- as.matrix(read.csv("Active Projects/Factor analysis and soil health/Data for GitHub/NY data - HL sil.csv"))
NY2 <- NY2.raw[,1]
NY2.tmp <- NY2.raw[,-1]
NY2.mat <- apply(NY2.tmp, 2, as.numeric)
NY2.mat <- as.matrix(NY2.mat)
rownames(NY2.mat) <- NY2
NY2.mat <- round(NY2.mat, 2)


#### EFA: NORTH CAROLINA (n = 68 observations) ####
## Choose number of factors to retain ##
ev.NC <- eigen(NC.matrix, symmetric=TRUE)
ap.NC <- parallel(subject=68,var=8, rep=100,cent=.05)
nS.NC <- nScree(x=ev.NC$values, aparallel=ap.NC$eigen$qevpea)
plotnScree(nS.NC)
#...all signs say 2 factors is preferred

# What do those factors look like? #
EFA.NC <- fa(NC.matrix, 2, rotate="oblimin", fm="ml")
EFA.NC
print(EFA.NC, digits=3, cut=.5)


#### CFA: NY1, no yield ####
CFA.NY.1 <- '
F1 =~ WSA + Soil.protein + Min.C + POXC
'
CFA.NY.1.fit <- cfa(CFA.NY.1, sample.cov = NY1.mat, sample.nobs = 16, std.lv = TRUE)
summary(CFA.NY.1.fit, fit.measures = TRUE, standardized = TRUE)
resid(CFA.NY.1.fit)

#### SEM: NY2, w/ yield ####
SEM.NY.2 <- '
F1 =~ WSA + Soil.protein + Min.C + POXC
Yield ~ F1
'
SEM.NY.2.fit <- sem(SEM.NY.2, sample.cov = NY2.mat, sample.nobs = 16, std.lv = TRUE)
summary(SEM.NY.2.fit, fit.measures = TRUE, standardized = TRUE)
resid(SEM.NY.2.fit)

## Checking CFA fit in NY2 ##
CFA.NY.2 <- '
F1 =~ WSA + Soil.protein + Min.C + POXC
'
CFA.NY.2.fit <- cfa(CFA.NY.2, sample.cov = NY2.mat, sample.nobs = 16, std.lv = TRUE)
summary(CFA.NY.2.fit, fit.measures = TRUE, standardized = TRUE)
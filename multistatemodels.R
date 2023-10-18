
##step1. simulation data######################
##the purpose is to simulate a dataset with co-occruance of diseases
#square-matrix to show transition
install.packages("mstate")
library("mstate")
tmat <- mstate::transMat(x = list(c(2, 3, 5, 6),
                                  c(4, 5, 6),
                                  c(4, 5, 6),
                                  c(5, 6),
                                  c(),
                                  c()),
                         names = c("Tx", "Rec", "AE", "Rec+AE",
                                   "Rel", "Death"))
print(tmat)

##         to
## from     Tx Rec AE Rec+AE Rel Death
##   Tx     NA   1  2     NA   3     4
##   Rec    NA  NA NA      5   6     7
##   AE     NA  NA NA      8   9    10
##   Rec+AE NA  NA NA     NA  11    12
##   Rel    NA  NA NA     NA  NA    NA
##   Death  NA  NA NA     NA  NA    NA

library("mstate")
data("ebmt4")
msebmt <- msprep(data = ebmt4, trans = tmat,
                 time = c(NA, "rec", "ae","recae", "rel", "srv"),
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
                 keep = c("match", "proph", "year", "agecl"))
msebmt[msebmt$id == 2, ]

## An object of class 'msdata'
##
## Data:
##    id from to trans Tstart Tstop time status              match proph
## 8   2    1  2     1      0    12   12      0 no gender mismatch    no
## 9   2    1  3     2      0    12   12      1 no gender mismatch    no
## 10  2    1  5     3      0    12   12      0 no gender mismatch    no
## 11  2    1  6     4      0    12   12      0 no gender mismatch    no
## 12  2    3  4     8     12    29   17      1 no gender mismatch    no
## 13  2    3  5     9     12    29   17      0 no gender mismatch    no
## 14  2    3  6    10     12    29   17      0 no gender mismatch    no
## 15  2    4  5    11     29   422  393      1 no gender mismatch    no
## 16  2    4  6    12     29   422  393      0 no gender mismatch    no
##         year agecl
## 8  1995-1998 20-40
## 9  1995-1998 20-40
## 10 1995-1998 20-40
## 11 1995-1998 20-40
## 12 1995-1998 20-40
## 13 1995-1998 20-40
## 14 1995-1998 20-40
## 15 1995-1998 20-40
## 16 1995-1998 20-40
install.packages("flexsurv")
library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE)
fits_wei <- vector(mode = "list", length = n_trans)
msebmt$years <- msebmt$time/365.25
for (i in 1:n_trans){
    fits_wei[[i]] <- flexsurvreg(Surv(years, status) ~ match + proph + year + agecl ,
                                 data = subset(msebmt, trans == i),
                                 dist = "weibull")
}


pat_2 <- data.frame(msebmt[msebmt$id == 2,
                           c("match", "proph", "year", "agecl")][1, ])
head(pat_2)

##                match proph      year agecl
## 8 no gender mismatch    no 1995-1998 20-40

yr_grid <- seq(0, 10, .1)

cumhaz_grid <- seq(0, max(msebmt$years), .01)
cumhaz_pat_2 <- msfit.flexsurvreg(fits_wei, trans = tmat,
                                  t = cumhaz_grid,
                                  newdata = pat_2,
                                  variance = FALSE)
head(cumhaz_pat_2$Haz)

##   time       Haz trans
## 1 0.00 0.0000000     1
## 2 0.01 0.1136587     1
## 3 0.02 0.1535433     1
## 4 0.03 0.1830811     1
## 5 0.04 0.2074241     1
## 6 0.05 0.2285137     1

tail(cumhaz_pat_2$Haz)

##        time       Haz trans
## 20431 16.97 0.3752974    12
## 20432 16.98 0.3754027    12
## 20433 16.99 0.3755079    12
## 20434 17.00 0.3756132    12
## 20435 17.01 0.3757183    12
## 20436 17.02 0.3758235    12

sim_stprobs_mstate_2 <- function(n_pats){
    mstate::mssample(Haz = cumhaz_pat_2$Haz,
                     trans = tmat,
                     tvec = yr_grid,
                     clock = "reset",
                     M = n_pats)
}

install.packages("msm")
library(msm)
statetable.msm(state, PTNUM, data=cav)
head(cav)


######################################
###step2-multistate models###########

##install packages
install.packages("ggplot2")
install.packages("tidymodels")
library(tidyverse)
library(tidymodels)
library(msm)
library(tidyr)
library(dplyr)
library(ggplot2)
##check data
head(cav)
set.seed(1234)
df <- tibble(cav) %>%
    mutate(o_state = state)

df

statetable.msm(state = state, subject = PTNUM, data = df)
##assume backward transition are misclassified and alter these variables
##b_age is the baseline age,


df1 <- df %>% group_by(PTNUM) %>%
    mutate(b_age = min(age),
           state = cummax(state)
    )
statetable.msm(state = state, subject = PTNUM, data = df1)

##Setting Up and Running the Model

# Intensity matrix Q:
q <- 0.01
Q <- rbind(c(0,q,0,q), c(0,0,q,q),c(0,0,0,q),c(0,0,0,0))
qnames <- c("q12","q14","q23","q24","q34")


covariates = list("1-2" = ~ years + b_age + dage ,
                  "1-4" = ~ years + b_age + dage ,
                  "2-3" = ~ dage,
                  "2-4" = ~ dage,
                  "3-4" = ~ dage)

obstype <- 1
#indicates that observations have been taken at arbitrary time points.
center <- FALSE
#means that covariates will not be centered at their means during the maximum likelihood estimation process.
deathexact <- TRUE
#indicates that the final absorbing state is exactly observed.
method <- "BFGS"
control <- list(trace = 0, REPORT = 1)

model_1 <- msm(state~years, subject = PTNUM, data = df1, center= center,
               qmatrix=Q, obstype = obstype, deathexact = deathexact, method = method,
               covariates = covariates, control = control)


#Model Status
##First check to see if the model has converged.
conv <- model_1$opt$convergence; cat("Convergence code =", conv,"\n")
## Convergence code = 0

##Next, look at a measure of how well the model fits the data proposed by using a visual test

prev <- prevalence.msm(model_1)

# reshape observed prevalence
do1 <-as_tibble(row.names(prev$Observed)) %>% rename(time = value) %>%
    mutate(time = as.numeric(time))
do2 <-as_tibble(prev$Observed) %>% mutate(type = "observed")
do <- cbind(do1,do2) %>% select(-Total)
do_l <- do %>% gather(state, number, -time, -type)
# reshape expected prevalence
de1 <-as_tibble(row.names(prev$Expected)) %>% rename(time = value) %>%
    mutate(time = as.numeric(time))
de2 <-as_tibble(prev$Expected) %>% mutate(type = "expected")
de <- cbind(de1,de2) %>% select(-Total)
de_l <- de %>% gather(state, number, -time, -type)

# bind into a single data frame
prev_l <-rbind(do_l,de_l) %>% mutate(type = factor(type),
                                     state = factor(state),
                                     time = round(time,3))


# plot for comparison
prev_gp <- prev_l %>% group_by(state)
pp <- prev_l %>% ggplot() +
    geom_line(aes(time, number, color = type)) +
    xlab("time") +
    ylab("") +
    ggtitle("")
pp + facet_wrap(~state)

# plot_prep was obtained from plot.msm()
res <- plot_prep(model_1)
time <- res[[1]]
Health <- res[[2]]
Mild_CAV <- res[[3]]
Severe_CAV <- res[[4]]
df_w <- tibble(time,Health, Mild_CAV, Severe_CAV)
df_l <- df_w %>% gather("state", "survival", -time)
p <- df_l %>% ggplot(aes(time, 1 - survival, group = state)) +
    geom_line(aes(color = state)) +
    xlab("Years") +
    ylab("Probability") +
    ggtitle("Fitted Survival Probabilities")
p

##########################################################
##################data simulation#########################
##########################################################
##https://rdrr.io/rforge/simfms/man/simfms.html
##simulation of clustered multistate data
install.packages("simfms", repos="http://R-Forge.R-project.org")
library(simfms)

simfms(nsim = NULL, tmat = NULL, clock = "forward",
       frailty = list(dist="gamma", par= .5, type="shared"),
       nclus = NULL, csize = NULL,
       covs = NULL, beta = NULL,
       marg = list(dist = "weibull", lambda = 1, rho = 1),
       cens = list(dist = "weibull", lambda = 1, rho = 1, admin = 72),
       copula = list(name = "clayton", par = 1),
       format = "long")


################################################################################
### - Cancer standard multi-state structure - ##################################
################################################################################
trans.cancer()
simfms(nsim  = NULL,
       tmat  = trans.cancer(),
       clock = "forward",
       frailty = list(dist="gamma", par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02, 8), treat=rep(2, 8)),
       marg  = list(dist="weibull", lambda=1, rho=1),
       cens  = list(dist="weibull", lambda=1, rho=1,  admin= 72),
       copula= list(name="clayton", par= 1))


##addPeriods created the longitudinal data,
dtTime <- addPeriods(dtTrial, nPeriods = 3, idvars = "id", timevars = c("Y0", "Y1",
                                                                        "Y2"), timevarName = "Y")
dtTime


def <- defData(varname = "xbase", dist = "normal", formula = 20, variance = 3)
def <- defData(def, varname = "nCount", dist = "noZeroPoisson", formula = 6)
def <- defData(def, varname = "mInterval", dist = "gamma", formula = 30, variance = 0.01)
def <- defData(def, varname = "vInterval", dist = "nonrandom", formula = 0.07)

dt <- genData(200, def)
dt[id %in% c(8, 121)]  # View individuals 8 and 121

dtPeriod <- addPeriods(dt)
dtPeriod[id %in% c(8, 121)]  # View individuals 8 and 121 only

def2 <- defDataAdd(varname = "Y", dist = "normal", formula = "15 + .1 * time", variance = 5)
dtPeriod <- addColumns(def2, dtPeriod)



#download the package
install.packages('msm')
library(msm)
#get the data
data(cav)  ##no need to download, can be used directly
cav[1:21,]
head(cav)
##1 CAV-free, 2 mild CAV, 3 moderate/sever CAV, 4 death
table(cav$statemax)
#statetable.msm provides a frequency table of pairs of consecutive states
statetable.msm(state, PTNUM, data=cav)


#to
#from    1    2    3    4
#1     1367  204   44   148
#2     46    134   54   48
#3     4     13    107  55


#specify a model

#O means no probability,
Q <- rbind ( c(0, 0.25, 0, 0.25),
             c(0.166, 0, 0.166, 0.166),
             c(0, 0.25, 0, 0.25),
             c(0, 0, 0, 0) )
Q

##initial model
Q.crude <- crudeinits.msm(state ~ years, PTNUM, data=cav,
                          qmatrix=Q)

cav.msm <- msm(state ~ years, subject=PTNUM, data=cav,
               qmatrix=Q, deathexact=4)

cav.msm

cavsex.msm <- msm(state ~ years, subject=PTNUM, data=cav,
                  qmatrix=Q, deathexact=4, covariates=~sex)


cavsex.msm

qmatrix.msm(cavsex.msm, covariates=list(sex=0)) # Male
qmatrix.msm(cavsex.msm, covariates=list(sex=1)) # Female

##trasition-specific covariates
cavsex.msm <- msm(state ~ years, subject=PTNUM, data=cav,
                  qmatrix=Q, deathexact=4,
                  covariates = list("1-2" = ~ sex, "1-4" = ~sex))

####????what is constraint??
#This constrains the effect of sex to be equal for the progression rates q12; q23, equal for the
#death rates q14; q24; q34, and equal for the recovery rates q21; q32.
cav3.msm <- msm(state ~ years, subject=PTNUM, data=cav,
                qmatrix=Q, deathexact=4,
                covariates=~sex,
                constraint = list(sex=c(1,2,3,1,2,3,2)))

##10 year transition probabilities, get probability
pmatrix.msm(cav.msm, t=10)

#           State 1    State 2    State 3   State 4
#State 1 0.30940656 0.09750021 0.08787255 0.5052207
#State 2 0.17165172 0.06552639 0.07794394 0.6848780
#State 3 0.05898093 0.02971653 0.04665485 0.8646477
#State 4 0.00000000 0.00000000 0.00000000 1.0000000
##a person in state1, who is disease-free, has a probability of 0.3 being still disease free, probabilities of 0.1 of being
#mild disease, and 0.5 of being dead in ten years.

sojourn.msm(cav.msm)
##total length of stay
totlos.msm(cav.msm)
qratio.msm(cav.msm, ind1=c(2,1), ind2=c(1,2))
hazard.msm(cavsex.msm)
#the estimated hazard ratios corresponding to each covariate effect on the transition intensities
#$sex, compard to men, the HR is 2.4 for women from state3 to state4
#                       HR            L            U
#State 1 - State 2 0.5632779042 3.333382e-01 9.518320e-01
#State 1 - State 4 1.1289701413 6.261976e-01 2.035418e+00
#State 2 - State 1 1.2905853501 4.916004e-01 3.388139e+00
#State 2 - State 3 1.0765518296 5.193868e-01 2.231408e+00
#State 2 - State 4 0.0003804824 7.241465e-65 1.999137e+57
#State 3 - State 2 1.0965531163 1.345395e-01 8.937364e+00
#State 3 - State 4 2.4135379727 1.176293e+00 4.952139e+00
#survival plot
plot(cav.msm, legend.pos=c(8, 1))


prevalence.msm(cav.msm, times=seq(0,20,2))
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)

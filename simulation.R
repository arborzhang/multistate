
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



#Fitting parametric survival models in R

# Libraries
library(flexsurvPlus)
library(tibble)
library(dplyr)
library(boot)
library(survival)
library(survminer)
library(tidyr)
library(ggplot2)

# Non-scientific notation
options(scipen=999) 

# Colours for plots
blue = rgb(69, 185, 209, max=255)
pink = rgb(211,78,147,max=255)
Dyellow = rgb(214, 200, 16, max=255)
orange<-rgb(247,139,21,max=255)

adtte <- sim_adtte(seed = 2020, rho = 0.6)
head(adtte)
#>   USUBJID ARMCD             ARM PARAMCD                     PARAM AVAL AVALU
#> 1       1     A Reference Arm A     PFS Progression Free Survival  141  DAYS
#> 2       2     A Reference Arm A     PFS Progression Free Survival  173  DAYS
#> 3       3     A Reference Arm A     PFS Progression Free Survival  197  DAYS
#> 4       4     A Reference Arm A     PFS Progression Free Survival  133  DAYS
#> 5       5     A Reference Arm A     PFS Progression Free Survival  100  DAYS
#> 6       6     A Reference Arm A     PFS Progression Free Survival  525  DAYS
#>   CNSR
#> 1    0
#> 2    0
#> 3    0
#> 4    0
#> 5    0
#> 6    0

# subset PFS data and rename
PFS_data <- adtte %>%
  filter(PARAMCD=="PFS") %>%
  transmute(USUBJID,
            ARMCD,
            PFS_days = AVAL,
            PFS_event = 1- CNSR
  )
head(PFS_data)
#>   USUBJID ARMCD PFS_days PFS_event
#> 1       1     A      141         1
#> 2       2     A      173         1
#> 3       3     A      197         1
#> 4       4     A      133         1
#> 5       5     A      100         1
#> 6       6     A      525         1

# Create survfit object
km.est.PFS <- survfit(Surv(PFS_days, PFS_event) ~ ARMCD, 
                      data = PFS_data, 
                      conf.type = 'plain')

# Plot Kaplan-Meier
plot(km.est.PFS, 
     col = c(blue, pink), # plot colours
     lty = c(1:2), # line type
     xlab = "Time (Days)", # x-axis label
     ylab = "Progression-free survival", # y-axis label
     xlim = c(0, 800)) 

legend(x = 500, 
       y = .9, 
       legend = c("Arm A", "Arm B"), 
       lty = c(1:2), 
       col = c(blue, pink))
       
       
 psm_PFS_all <- runPSM(data=PFS_data,
                     time_var="PFS_days",
                     event_var="PFS_event",
                     model.type= c("Common shape", 
                                   "Independent shape", 
                                   "Separate"),
                     distr = c('exp',
                               'weibull',
                               'gompertz',
                               'lnorm',
                               'llogis',
                               'gengamma',
                               'gamma',
                               'genf'),
                     strata_var = "ARMCD",
                     int_name="A",
                     ref_name = "B")
psm_PFS_all

## ----setup, message=FALSE-----------------------------------------------------
library(pspline.inference)

library(mgcv)
library(magrittr)
library(dplyr)

library(kableExtra)

# Housekeeping for publishing on CRAN
source("cran.R")

# Housekeeping for parallel computation
source("parallel.R")

# Customizations for figures
source("figures.R")

## ----data---------------------------------------------------------------------
obs <- read.csv(system.file("extdata", "seasonal.csv", package="pspline.inference"))
obs$cases.cum <- cumsum(obs$cases)
obs$cases.cumrel <- obs$cases.cum / max(obs$cases.cum)

## ----model--------------------------------------------------------------------
model <- gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=obs)

## ----model time---------------------------------------------------------------
sampleFreq <- 20
modelTime <- seq(min(obs$time) - 0.5, max(obs$time) + 0.5 - 1 / sampleFreq, 1 / sampleFreq)
predictors = data.frame(time=modelTime)

## ----model samples------------------------------------------------------------
n <- 20

## ----estimated cases----------------------------------------------------------
casesEst <- model %>%
  pspline.estimate.timeseries(predictors, pspline.outbreak.cases, samples=n, level=0.95)

## ----estimated cases plot-----------------------------------------------------
ggplot(casesEst, c("casesObs", "casesEst")) +
  geom_ribbon(aes(x=time, ymin=cases.lower, ymax=cases.upper, fill="casesEst")) +
  geom_point(data=obs, aes(x=time, y=cases, color="casesObs", size="casesObs")) +
  labs(x=NULL, y="Cases")

## ----estimate cum cases, eval=vignette.eval$full------------------------------
cumCasesRelEst <- model %>%
  pspline.estimate.timeseries(predictors, pspline.outbreak.cumcases.relative, samples=n, level=0.95)

## ----estimate cum cases plot, eval=vignette.eval$full-------------------------
ggplot(cumCasesRelEst, c("casesObs", "casesEst")) +
  geom_ribbon(aes(x=time, ymin=cases.cumrel.lower, ymax=cases.cumrel.upper, fill="casesEst")) +
  geom_point(data=obs, aes(x=time, y=cases.cumrel, color="casesObs", size="casesObs")) +
  labs(x=NULL, y="Relative cumulative incidence")

## ----plot cases estimates, eval=vignette.eval$full----------------------------
casesSamples <- model %>%
  pspline.sample.timeseries(predictors, pspline.outbreak.cases, samples=15)
  
ggplot(casesSamples, c("casesObs", "casesEstSamples")) +
  geom_line(aes(x=time, y=cases, group=pspline.sample, color="casesEstSamples")) +
  geom_point(data=obs, aes(x=time, y=cases, color="casesObs", size="casesObs")) +
  labs(y="Cases")

## ----custom time series calc--------------------------------------------------
seriousDiseaseCases <- function(model, data) {
  # Get predicted incidence and cumulative incidence for the given param values
  cumcases = pspline.outbreak.calc.cumcases(data$time, data$cases)
  
  # Calculate when we hit 50 cases
  treatmentAvail <- cumcases < 50

  # Calculate serious disease case counts -- 80% of (10% before we run out of treatment + 100% after). 
  seriouscases = data$cases * 0.8
  seriouscases[treatmentAvail] = seriouscases[treatmentAvail] * 0.1
  
  data %>%
    mutate(seriouscases.cum=pspline.outbreak.calc.cumcases(time, seriouscases)) %>%
    select(-cases)
}

## ----custom time series est---------------------------------------------------
seriousEst <- model %>%
  pspline.estimate.timeseries(predictors, seriousDiseaseCases, samples=n, level=0.95)

## ----custom cum cases est-----------------------------------------------------
cumCasesEst <- model %>%
  pspline.estimate.timeseries(predictors, pspline.outbreak.cumcases, samples=n, level=0.95)

## ----custom time series est plot----------------------------------------------
ggplot(seriousEst, c("casesObs", "casesEstMedian", "seriousEstMedian")) +
  geom_line(aes(x=time, y=seriouscases.cum.median, color="seriousEstMedian")) +
  geom_line(data=cumCasesEst, aes(x=time, y=cases.cum.median, color="casesEstMedian")) +
  geom_point(data=obs, aes(x=time, y=cases.cum, color="casesObs", size="casesObs")) +
  labs(y="Cases")

## ----custom calc parameterized, eval=vignette.eval$full-----------------------
seriousDiseaseCases <- function(progressionRate, treatmentSuccessRate, treatmentMax) {
  function(model, data) {
    # Calculate cumulative incidence
    cumcases = pspline.outbreak.calc.cumcases(data$time, data$cases)
    
    # Calculate when we hit 100 cases
    treatmentAvail <- cumcases < treatmentMax
  
    # Calculate serious disease case counts -- 80% of (10% before we run out of treatment + 100% after). 
    seriouscases = data$cases * progressionRate
    seriouscases[treatmentAvail] = seriouscases[treatmentAvail] * (1 - treatmentSuccessRate)
    
    data %>%
      mutate(seriouscases.cum=pspline.outbreak.calc.cumcases(time, seriouscases)) %>%
      select(-cases)
  }
}

## ----custom est parameterized, eval=vignette.eval$full------------------------
seriousEst <- model %>%
  pspline.estimate.timeseries(
    predictors, 
    seriousDiseaseCases(progressionRate=0.8, treatmentSuccessRate=0.9, treatmentMax=50), 
    samples=n, 
    level=0.95
  )

## ----estimate supply duration, message=FALSE----------------------------------
medicationSupplyEnd <- function(medicationQuantity) {
  function(model, data) {
    # Get predicted cumulative incidence for the given param values
    data %<>% mutate(
      cumcases=pspline.outbreak.calc.cumcases(time, cases)
    )
    
    # Find where cumulative incidence exceeds the amount of medication we have
    idx = sum(data$cumcases < medicationQuantity)
    
    # Use the point halfway between latest time below the threshold and the earliest time above the threshold
    data.frame(supply.duration=mean(c(data$time[idx], data$time[idx + 1])))
  }
}

supplyDurationEst <- model %>% 
  pspline.estimate.scalars(predictors, medicationSupplyEnd(medicationQuantity=50))

## ----estimate supply duration table, echo=TRUE--------------------------------
kable_styling(
  kable(
    supplyDurationEst, 
    caption="Estimated time before medication supply is exhausted, in weeks since July 1st", 
    col.names=c("Lower CL", "Median", "Upper CL")
  )
)

## ----plot supply duration, warning=FALSE, eval=vignette.eval$full-------------
supplyDurationSamples <- model %>% 
  pspline.sample.scalars(predictors, medicationSupplyEnd(medicationQuantity=50), samples=n)

ggplot(seriousEst, c("casesObs", "casesEst", "supplyEnd")) +
  geom_ribbon(data=cumCasesEst, aes(x=time, ymin=cases.cum.lower, ymax=cases.cum.upper, fill="casesEst", color="casesEst")) +
  geom_point(data=obs, aes(x=time, y=cases.cum, color="casesObs", size="casesObs")) +
  geom_density(
    data=supplyDurationSamples,
    aes(x=supply.duration, y=..density..*50, fill="supplyEnd", color="supplyEnd"),
    trim=TRUE
  ) +
  labs(x=NULL, y="Cumulative incidence")

## ----generate truth, eval=vignette.eval$full----------------------------------

generateTruths <- function() {
  timeMin = 1
  timeMax = 52
  deltaT = 0.05

  # Times
  time = seq(timeMin - 0.5, timeMax + 0.5 - deltaT, deltaT)
  
  # Rough shape of a peak using a single cos period
  onset = runif(1, 10, 40)
  duration = runif(1, 5, 15)
  max = runif(1, 25, 75)
  cases0 = 0.5 * log(max) * (1 - cos((time - onset) * 2 * pi / duration)) * (time > onset) * (time < onset + duration)
  cases0 = round(exp(cases0) - 1)
  data0 = data.frame(time=time, cases=cases0)

  # Make it polynomial by running it through a simple spline
  data.frame(
    time=time,
    cases=gam(cases ~ s(time, k=4, bs="cp", m=3), family=poisson, data=data0) %>%
      predict(type="response")
  )
}

## ----generate observations, eval=vignette.eval$full---------------------------
generateObservations <- function(truth) {
  # Initialize with observation sampling times
  observed = data.frame(time = seq(min(truth$time) + 0.5, max(truth$time), 1))
  
  # Merge with true case values and run through rpoiss
  truth %>% 
    inner_join(observed, by="time") %>%
    mutate(cases=rpois(nrow(.), cases))
}

## ----make model, eval=vignette.eval$full--------------------------------------
makeModel <- function(data) {
  gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=data)
}

## ----validate supply duration, message=FALSE, warning=FALSE, include=FALSE, eval=vignette.eval$full----
set.seed(0)
validationResults = pspline.validate.scalars(
  fun.truth=generateTruths,
  n.truths=10,
  fun.observations=generateObservations,
  n.observations=10,
  fun.model=makeModel,
  fun.outcome=medicationSupplyEnd(50),
  n.samples=20,
  level=0.95
)

## ----validate supply duration table, eval=vignette.eval$full------------------
kable_styling(
  kable(
    sprintf("%.2f%%", validationResults$summary$supply.duration.good * 100), 
    caption="Coverage: frequency with which true value is included in the 95% CI", 
    col.names=c("Supply duration")
  )
)


library(dplyr)
library(doParallel)
library(foreach)
library(beepr)
library(here)


source(here("Vfunctions", "vStateUpdate_withHomogeneity.R"))
source(here("Vfunctions", "V_DS_T2.R"))
source(here("Vfunctions", "vSI_withHomogeneity.R"))
source(here("Vfunctions", "effcChecker.R"))
source(here("Vfunctions", "effcIdGenT2.R"))
source(here("Vfunctions", "reseedingT2.R"))

# Set Up for Parallel Computing

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


## Model Definition ##
runs <- 1
generations <- 100 # the number of generations (or timesteps) to run the model over
output <- "counts" # defines the output as "matrix" (full disease state and ID of units) or "counts" (number of units in S,I,R states or "I_counts" (returns only count of infected units over time)

## PARAMETERS ##

## Population Parameters ##

N <- 1000 # population size (may exxclude in favour of calculating in state_update function)
I0 <- 10
state_tracker <- c(rep("I", I0), rep("S", N - I0)) # initial disease state of units in population
population_sample <- "large" # large or small population (determines sampling strategy)

## Epidemiological Parameters ##

R_period <- Inf # the refractory period: 1 = no immunity, Inf = lifelong immunity
R0_estimate <- 3 # average number of effective contacts per unit
k_transmission <- Inf # overdispersion parameter, level of heterogeneity in effective contact potential: 0.1 = high heterogneiety, Inf = homogeneous

heterogeneity <- "fixed" # "fixed" or "variable", is unit transmission potential constant through time?
contact_list1 <- replicate(N, c(1:N), simplify = FALSE) 
contact_list <- lapply(1:length(contact_list1), 
                       function(x) contact_list1[[x]][-x]) 

## Reseeding Parameters ##
reseed <- FALSE # reseed 5 infections periodically during the first 10% of generations
reseed_threshold <- 0# threshold generation for reseeding events
reseed_frequency <- 0
introductions <- 0 # number of new cases introduced when reseeding occurs

## Vaccination Parameters ##
V_prop <- 0 # proportion to be vaccinated in each round
V_start <- 0 # start generation for vaccination
V_period <- 0 # duration of immunity from vaccination (by number of timesteps)
V_rounds <- 0 # number of rounds (-1 for starting round) of vaccination in programme
V_schedule <- 0 # how frequent is vacciantion (ie every 2 "years")
V_mode <- "random" #random, target or non_target, NB: target/non-target only possible if heterogneiety is fixed


## SIMULATION ##




Runs_tracker <-
  state_update(
    
    ## Population Parameters ##
    state_tracker, # a vector of the the initial state of the population
    N, # population size (may exxclude in favour of calculating in state_update function)
    population_sample, # large or small population (determines sampling strategy)
    
    ## Model Parameters ##
    generations, # the number of generations (or timesteps) to run the model over
    output, # defines the output as "matrix" (full disease state and ID of units) or "counts" (number of units in S,I,R states or "I_counts" (returns only count of infected units over time)
    
    ## Epidemiological Parameters ##
    R_period, # the refractory period (number of timesteps that units are immune)
    R0_estimate, # average number of effective contacts per unit
    k_transmission, # overdispersion parameter, level of heterogeneity in effective contact potential
    heterogeneity, # is heterogeneity in effective contacts "fixed" or "variable"
    contact_list, # potential contacts of each unit 
    
    ## Reseeding Parameters ##
    reseed, # reseed 5 infections periodically during the first 10% of generations
    reseed_threshold, # threshold generation for reseeding events
    reseed_frequency,
    introductions, # number of new cases introduced when reseeding occurs
    
    ## Vaccination Parameters ##
    V_prop, # proportion to be vaccinated in each round
    V_start, # start generation for vaccination
    V_period, # duration of immunity from vaccination (by number of timesteps)
    V_rounds, # number of rounds (-1 for starting round) of vaccination in programme
    V_schedule, # how frequent is vacciantion (ie every 2 "years")
    V_mode #random, target or non_target
  )
















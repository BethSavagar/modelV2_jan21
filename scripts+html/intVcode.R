
## Model Definition ##
runs <- 500
generations <- 2000 # the number of generations (or timesteps) to run the model over
output <- "I_counts" # defines the output as "matrix" (full disease state and ID of units) or "counts" (number of units in S,I,R states or "I_counts" (returns only count of infected units over time)

## Population Parameters ##

N <- 10000 # population size (may exxclude in favour of calculating in state_update function)
I0 <- 10
state_tracker <- c(rep("I", I0), rep("S", N - I0)) # initial disease state of units in population
population_sample <- "large" # large or small population (determines sampling strategy)

## Epidemiological Parameters ##

R_period <- 1 # the refractory period: 1 = no immunity, Inf = lifelong immunity
R0_estimate <- 3 # average number of effective contacts per unit
k_transmission <- Inf # overdispersion parameter, level of heterogeneity in effective contact potential: 0.1 = high heterogneiety, Inf = homogeneous

heterogeneity <- "fixed" # "fixed" or "variable", is unit transmission potential constant through time?

if(population_sample == "small"){
  contact_list1 <- replicate(N, c(1:N), simplify = FALSE)
  contact_list <- lapply(1:length(contact_list1), function(x) contact_list1[[x]][-x])
}else{
  contact_list <- vector()
}


## Reseeding Parameters ##
reseed <- TRUE # reseed 5 infections periodically during the first 10% of generations
reseed_threshold <- 200 # threshold generation for reseeding events
reseed_frequency <- 20
introductions <- 10 # number of new cases introduced when reseeding occurs

## Vaccination Parameters ##
V_prop <- 0 # proportion to be vaccinated in each round
V_start <- 1500 # start generation for vaccination
V_period <- Inf # duration of immunity from vaccination (by number of timesteps)
V_rounds <- 1 # number of rounds (-1 for starting round) of vaccination in programme
V_schedule <- 0 # how frequent is vacciantion (ie every 2 "years")
V_mode <- "random" #random, target or non_target, NB: target/non-target only possible if heterogneiety is fixed


# With 2000 generations ##
het <- rep(c("fixed", "fixed", "fixed", "homogeneous"), 11)
k <- rep(c(0.1,0.3,Inf,Inf), 11)
pVac <- rep(seq(0, 1, 0.1), each = 4)

pars <- data.frame("k" = k, 
                   "het" = het, 
                   "pVac" = pVac)

R0_estimate <- 3
R_period <- 1
generations <- 2000
runs <- 500
V_start <- 1500

## Model ##

intV_list <- list()

system.time(
  
  for (i in 1:nrow(pars)) {
    
    k_transmission <- pars[i, 1]
    heterogeneity <- pars[i, 2]
    V_prop <- pars[i, 3]
    
    print(i)
    
    Runs_tracker <-
      foreach (z = 1:runs, .combine = "rbind") %dopar% {
        state_update(
          state_tracker,# a vector of the the initial state of the population
          N, # population size (may exxclude in favour of calculating in state_update function)
          population_sample,# large or small population (determines sampling strategy)
          
          generations,# the number of generations (or timesteps) to run the model over
          output,# defines the output as "matrix" (full disease state and ID of units) or "counts" (number of units in S,I,R states or "I_counts" (returns only count of infected units over time)
          
          R_period, # the refractory period (number of timesteps that units are immune)
          R0_estimate, # average number of effective contacts per unit
          k_transmission, # overdispersion parameter, level of heterogeneity in effective contact potential
          heterogeneity, # is heterogeneity in effective contacts "fixed" or "variable"
          contact_list, # potential contacts of each unit
          
          reseed, # reseed 5 infections periodically during the first 10% of generations
          reseed_threshold, # threshold generation for reseeding events
          reseed_frequency,
          introductions, # number of new cases introduced when reseeding occurs
          
          V_prop, # proportion to be vaccinated in each round
          V_start, # start generation for vaccination
          V_period, # duration of immunity from vaccination (by number of timesteps)
          V_rounds, # number of rounds (-1 for starting round) of vaccination in programme
          V_schedule, # how frequent is vacciantion (ie every 2 "years")
          V_mode #random, target or non_target
        )
        
      }
    
    
    intV_list[[i]] <- Runs_tracker
    
  }
)

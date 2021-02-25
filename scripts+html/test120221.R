
pars <- cbind("k_transmission" = rep(c(0.1,0.3,Inf,Inf), 11), 
              "heterogeneity" = rep(c("fixed", "fixed", "fixed", "homogeneous"), 11), 
              "V_prop" = rep(seq(0.6, 0.7, 0.01), each = 4))


R0_estimate <- 3
R_period <- 1

VList_main <- list()


system.time(
  
  for (i in 1:nrow(pars)) {
    
    k_transmission <- pars[i, 1]
    heterogeneity <- pars[i, 2]
    V_prop <- pars[i, 3]
    
    print(k_transmission)
    print(heterogeneity)
    print(V_prop)
    
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
    
    
    VList_main[[i]] <- Runs_tracker
    
  }
)




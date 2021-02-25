# Updated with vaccination and working as per 11.01.21

# given the initial disease state of the population, the state_tracker function will return the change in the disease state of individuals in the population over time

state_update <- function(state_tracker, # vector of initial population states
                         generations, # number of timesteps/generations to simulate
                         R_period, # the immune period (number of timesteps that a unit spends in the recovered state)
                         N, # population size
                         heterogeneity, # is heterogeneity in effective contacts "fixed" or "variable"
                         k_estimate, # overdispersion parameter, level of heterogeneity in effective contact potential
                         R_estimate, # average number of effective contacts per unit
                         output, # defines the output as "matrix" or "counts" or "I_counts" (number of infections over time)
                         reseed, # defines whether to re-seed new infections in early generations (taking arguments "TRUE" or "FALSE")
                         threshold, # the threshold generation for re-seeding events, if reseed = "TRUE"
                         frequency, # the frequency (number of generations) of reseeding events, if reseed = "TRUE
                         introductions, # number of introductions to seed if re-seeding = "TRUE"
                         contact_list, # list of potential contacts for each unit in the popualtion
                         population, # "large" or "small", defines how to sample cases from the population (if large sample from N. if small sample from N - index case)
                         p_vac, # proportion to be vaccinated in each round
                         V_start, # start generation for vaccination
                         V_period, # duration of immunity from vaccination (by number of timesteps)
                         V_rounds, # number of rounds of vaccination in programme
                         V_schedule, # how frequent is vacciantion (ie every 2 "years")
                         V_mode  # mode of vacc (random, target, non-target)
) {
  
  
  # -------------
  ## SET UP ##
  # -------------
  
  # Create a vector to store the number of timesteps spent in R state
  # this will ensure that units spend the number of generations defined by R_period in the "R" state before reverting to susceptibility
  RTime <- c(rep(0, N))
  
  VTime <- c(rep(0, N))
  
  V_timing <- seq(from = V_start, to = V_start+((V_rounds-1)*V_schedule), by = V_schedule)
  
  # If heterogeneity == "fixed"
  # Define the effective contact potential of units, when heterogeneity == "fixed"
  # If the number of effective contacts per unit is fixed over time then define a vector "fixed_effc" to contain the number of contacts per unit in the population.
  fixed_effc <- rnbinom(N, # population size
                        size = k_estimate, # dispersion parameter (level of heterogeneity)
                        mu = R_estimate) # R0
  
  # Use effcChecker function to validate that the number of effective contacts per unit does not exceed the total population size
  fixed_effc <- effcChecker(fixed_effc, N, k_estimate, R_estimate)
  
  
  # Set the initial state of the population, this will be updated with each iteration of the generations FOR loop below.
  
  disease_state <- state_tracker
  
  
  # Create vectors to store the number of units in each state per generation, only if the output is "counts" / "I_counts"
  
  if(output == "counts" | 
     output == "I_counts"){
    S_counts <- sum(disease_state == "S")
    I_counts <- sum(disease_state == "I")
    R_counts <- sum(disease_state == "R")
    V_counts <- sum(disease_state == "V")
  }
  
  # --------------
  ## SIMULATION ##
  # --------------
  
  for (gen in 2:generations) { # run the simulation for the number of timsteps specified in 'generations' variable
    
    # new disease state takes a vector with the current state of the population and simulates a state change for each unit, producing a new disease state vector as output. 
    new_disease_state <- newDiseaseState(disease_state,
                                         R_period,
                                         N,
                                         heterogeneity,
                                         k_estimate,
                                         R_estimate,
                                         output,
                                         reseed, # reseed infection periodically during the first 10% of generations
                                         RTime, # tracks the number of timesteps spent in "R" state for each individual unit
                                         fixed_effc,
                                         gen, # tracks which generation is being simulated
                                         threshold, 
                                         frequency, 
                                         introductions, 
                                         contact_list,
                                         population,
                                         p_vac, # proportion to be vacciantion in each round
                                         V_period, # duration of immunity from vaccination (by number of timesteps)
                                         V_timing, # vector defining timesteps for vaccination
                                         V_mode,  # mode of vacc (random, target, non-target)
                                         VTime # vector tracking hte number of generations each unit spends in V state
    )
    
    # ---------------------------------
    ## UPDATE RTime, VTime##
    # ---------------------------------
    
    
    # Update RTime vector according to new state of units
    RTime[new_disease_state == "R"] <- RTime[new_disease_state == "R"] +1
    RTime[new_disease_state != "R"] <- 0
    
    # Update VTime vector according to new state of units
    #VTime[newV_index] <- 0
    
    VTime[new_disease_state == "NV"] <- 0
    new_disease_state[new_disease_state == "NV"] <- "V"
    
    VTime[new_disease_state == "V"] <- VTime[new_disease_state == "V"] +1
    VTime[new_disease_state != "V"] <- 0
    
    
    # ---------------------------------
    ## UPDATE DISEASE STATE (matrix/counts) ##
    # ---------------------------------
    
    # Update tracker of units ID and disease states/numbers of units in each state across generations
    
    
    if (output == "matrix") {
      state_tracker <-
        cbind(state_tracker, new_disease_state) # update column of state_matrix with new population state
      
    } else if (output == "counts") {
      # update vectors containing counts of units in each disease state over time
      
      S_counts <- c(S_counts, sum(new_disease_state == "S"))
      I_counts <- c(I_counts, sum(new_disease_state == "I"))
      R_counts <- c(R_counts, sum(new_disease_state == "R"))
      V_counts <- c(V_counts, sum(new_disease_state == "V"))
      
    } else if (output == "I_counts") {
      I_counts <- c(I_counts, sum(new_disease_state == "I"))
    }
    
    # update the disease state of the population for the next loop
    disease_state <- new_disease_state
  } 
  ## END OF FOR LOOP
  
  # -----------------------
  ## DEFINE OUTPUT ##
  # -----------------------
  
  # Define the formate of the function output, as specified by the 'output' argument
  # If output is set to 'matrix' return the full matrix containing the disease state of each individual unit over time
  # If output is set to 'counts' return a dataframe tracking only the number (not the identity) of units in each disease state over time
  
  
  if (output == "matrix") {
    
    #rename state_tracker columns as generation number
    colnames(state_tracker) <- seq(1:generations)
    
    return(state_tracker) # return the full matrix tracking the state of each unit in the population across time
    
  } else if (output == "I_counts") {
    
    return(I_counts)  # returns just the number of units in the infected state over time (tracks only epidemic trajectory)
    
  } else if (output == "counts") { # return the number (but not ID) of units in each state across time
    
    # bind state_count vectors into a data_frame called state_counts, with each disease state occupying a different row (S, I, R)
    state_counts <- as.data.frame(rbind(S_counts, I_counts, R_counts, V_counts))
    # rename columns and rows of state_counts 
    colnames(state_counts) <- seq(1:generations)
    rownames(state_counts) <- c("S", "I", "R", "V")
    
    return(state_counts) # return the number (but not ID) of units in each state across time
  }
  
}

# Updated with vaccination and working as per 11.01.21

# given the initial disease state of the population, the state_tracker function will return the change in the disease state of individuals in the population over time

state_update <- function(
  
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
) {
  
  
  # -------------
  ## SET UP ##
  # -------------
  
  # Create a vector to store the number of timesteps spent in R state
  # this will ensure that units spend the number of generations defined by R_period in the "R" state before reverting to susceptibility
  RTime <- c(rep(0, N))
  
  VTime <- c(rep(0, N))
  
  V_timing <- seq(from = V_start, # vaccination start generation
                  to = V_start+((V_rounds-1)*V_schedule), # start generation + number of rounds* frequency of rounds 
                  by = V_schedule) # frequency of vaccination
  
  # Define Effective Contact potentials for Fixed Heterogeneity or Homogeneous population
  
  # If the number of effective contacts per unit is fixed over time, define vector "fixed_effc" to contain the number of contacts per unit in the population.
  
  if(heterogeneity == "fixed"){
    fixed_effc <- rnbinom(N, # population size
                          size = k_transmission, # dispersion parameter (level of heterogeneity)
                          mu = R0_estimate) # R0
    
    # effcChecker validates that the max number of effective contacts per unit does not exceed the total population size
    fixed_effc <- effcChecker(fixed_effc, N, k_transmission, R0_estimate)
    
  } else if(heterogeneity == "homogeneous"){ # if homogeneous then set contact rate for each unit in the population to R0
    
    fixed_effc <- rep(R0_estimate, N)
  }else if(heterogeneity == "variable"){
    fixed_effc <- vector()
  }
  
  
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
                                         N,
                                         population_sample,
                                         
                                         gen, # tracks which generation is being simulated
                                         R_period,
                                         RTime, # tracks the number of timesteps spent in "R" state for each individual unit
                                         R0_estimate,
                                         k_transmission,
                                         heterogeneity,
                                         contact_list,
                                         fixed_effc,
                                         
                                         reseed, # reseed infection periodically during the first 10% of generations
                                         reseed_threshold, 
                                         reseed_frequency, 
                                         introductions, 
                                         
                                         V_prop, # proportion to be vacciantion in each round
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
    state_counts <- as.data.frame(rbind(S_counts, 
                                        I_counts, 
                                        R_counts, 
                                        V_counts))
    
    # rename columns and rows of state_counts 
    colnames(state_counts) <- seq(1:generations)
    rownames(state_counts) <- c("S", "I", "R", "V")
    
    return(state_counts) # return the number (but not ID) of units in each state across time
  }
  
}

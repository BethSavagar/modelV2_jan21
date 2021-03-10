# Updated with homogeneity and working as per 19/02/21
# Investigsting addition of exposure variaiton (19/-02/21)
# adding new exposure correlation code from 09 Mar 2021 - see newExposurecode_010321.Rmd

# given the initial disease state of the population, the state_tracker function will return the change in the disease state of individuals in the population over time

state_update <- function(
  
  ## Population Parameters ##
  state_tracker, # a vector of the the initial state of the population
  N, # population size (may exxclude in favour of calculating in state_update function)
  population_sample, # large or small population (determines sampling strategy)
  contact_list, # potential contacts for each unit (when population_sample == small)
  
  ## Model Parameters ##
  generations, # the number of generations (or timesteps) to run the model over
  output, # defines the output as "matrix" (full disease state and ID of units) or "counts" (number of units in S,I,R states or "I_counts" (returns only count of infected units over time)
  
  ## Epidemiological Parameters ##
  R_period, # the refractory period (number of timesteps that units are immune)
  R0_estimate, # average number of effective contacts per unit
  k_transmission, # defines overdispersion of contacts in population (transmission and exposure)
  heterogeneity, # population variation is “fixed” (constant over time), “variable” (changes each timestep) or “homogeneous” (all units have same potential)
  fixed_exp, # TRUE or FALSE, if TRUE then exposure potential is equal for all units, if FALSE then exposure potential is drawn from NBD as transmission
 
  
  # Exposure and Transmission Correlation
  C_lower, # lower bound for correlation range
  C_upper, # upper bound for correlation range
  sd0, # starting standard deviation, usually 0
  sd_increment, # increment for increasing standard deviation in while loop used to generate EXPdist
  
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
  
 
  ############################################################################
  ############# EXPOSURE AND TRANSMISSION POTENTIAL ##################

  if(heterogeneity == "fixed"){
    effc_dist <- rnbinom(N, # population size
                          size = k_transmission, # dispersion parameter (level of heterogeneity)
                          mu = R0_estimate) # R0
    
    # effcChecker validates that the max number of effective contacts per unit does not exceed the total population size
    effc_dist <- effcChecker(effc_dist, N, k_transmission, R0_estimate)
    
    ## Exposure Potential ##
    
    correlation <- 1
    
    if(C_lower > 0){ # if correlation is positive
      while(correlation < C_lower | correlation > C_upper){
        
        R0rnd <- effc_dist + rnorm(N,
                                   mean=R0_estimate,
                                   sd=sd0)
        
        exp_dist <- c()
        exp_dist[order(R0rnd)] <- sort(effc_dist) # sort both from small to large (positive correlation)
        
        correlation <- cor(effc_dist, exp_dist)
        sd0 <- sd0+sd_increment # increase sd to decrease correlation
        
      }
    }else if(C_upper < 0){ # if correlation is negative
      
      R0rnd <- effc_dist + rnorm(N,
                                 mean=R0_estimate,
                                 sd=sd0)
      
      exp_dist <- c()
      exp_dist[order(R0rnd)] <- rev(sort(effc_dist))
      
      correlation <- cor(effc_dist, exp_dist)
      sd0 <- sd0+sd_increment # increase sd to decrease correlation
    } 
    
    
  } else if(heterogeneity == "homogeneous"){ # if homogeneous then set contact rate for each unit in the population to R0
    effc_dist <- rep(R0_estimate, N) # all units transmission potential == R0
    exp_dist <- rep(1, N) # all units exposure potential == 1
    
  }else if(heterogeneity == "variable"){
    effc_dist <- vector() # empty vector will be populated every 'generation'
    exp_dist <- vector() # empty vector will be populated every 'generation'
  }
  
  if(fixed_exp == TRUE){
    exp_dist <- rep(1,N)
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
                                         contact_list,
                                         
                                         
                                         gen, # tracks which generation is being simulated
                                         R_period,
                                         RTime, # tracks the number of timesteps spent in "R" state for each individual unit
                                         R0_estimate,
                                         k_transmission,
                                        
                                         heterogeneity,
                                         fixed_exp,

                                         effc_dist,
                                         exp_dist,
                                         
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

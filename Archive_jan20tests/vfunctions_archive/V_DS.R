
# Updated with vaccination and working as per 11.01.21


#newDiseaseState function takes as input the current state of the population, produces as output a vector containing the new state of the population


newDiseaseState <-
  function(disease_state,
           R_period,
           N,
           heterogeneity,
           k_estimate,
           R_estimate,
           output,
           reseed, 
           RTime,
           fixed_effc,
           gen, 
           threshold, 
           frequency, 
           introductions, 
           contact_list,
           population, # the effective contact potential of each unit when heterogeneity is fixed
           p_vac, # proportion to be vacciantion in each round
           # start_date, # start generation for vaccination
           V_period, # duration of immunity from vaccination (by number of timesteps)
           # N_rounds, # number of rounds of vaccination in programme
           V_timing, # vector defining timesteps for vaccination
           V_mode,  # mode of vacc (random, target, non-target)
           VTime # vector tracking hte number of generations each unit spends in V state
  ){
  
  S_index <- which(disease_state == "S") # ids of Susceptible individuals
  I_index <- which(disease_state == "I") # ids of Infected individuals
  R_index <- which(disease_state == "R") # ids of Recovered individuals
  V_index <- which(disease_state == "V") # ids of Vacc individuals
  
  # Vector to hold new disease state
  new_disease_state <-
    vector(length = length(disease_state)) # vector to hold disease state of units in population at next timestep
  
  
  # ------------
  ## RECOVERY ##
  # ------------
  
  # All infecteds become recovered in a given timestep
  
  if (length(I_index) > 0) {
    new_disease_state[I_index] <- "R"
  }
  
  
  
  
  
  # -------------------------------
  ## LOSS OF IMMUNITY ##
  # -------------------------------
  
  ## Recovereds become Susceptible again after being in R_state for >= R_period timesteps
  
  if (length(R_index) > 0) {
   # RTime[R_index] <- RTime[R_index] + 1 # update R units in RTime vector to reflect +1 timestep in R state
    newS_index <- which(RTime == R_period) # identify units which have been in R state for >= ~R_period, store in newS_index
    
    # Update RTime vector for units reverting to susceptibility
  #  RTime[newS_index] <- 0 # becomes 0 since units are no longer immune
    
    # Update new_disease_state vector
    new_disease_state[R_index] <- "R"
    new_disease_state[newS_index] <- "S" # units identified in newS_index lose immunity
  }
  
  
  if (length(V_index) > 0) {
    newS_index <- which(VTime == V_period) 
   
    # Update new_disease_state vector
    new_disease_state[V_index] <- "V"
    new_disease_state[newS_index] <- "S" # units identified in newS_index lose immunity
  }
  
  
  
  # -----------
  ## INFECTION
  # -----------
  
  ## Infected units produce secondary cases according to negative binomial distribution (NBD)
  # See the S_to_I function description above for a detailed explanation
  # S_to_I takes a vector of susceptible units and a vector of infected units and generated new infected units (cases) based on anegative binomial distribution. 
  # S_to_I outputs a vector (cases_index) storing the ID of new infected units (cases)
  
  if (length(S_index > 0)) {
    cases_index <- S_to_I(I_index, # S_to_I generates a vector of new effective contacts i.e. new cases
                          S_index,
                          N,
                          k_estimate,
                          R_estimate,
                          heterogeneity, # heterogeneity variable determines whether unit effective contacts are fixed or variable over time.
                          fixed_effc, 
                          contact_list,
                          population # stores the effective contact potential of units, used if heterogeneity == "fixed"
    ) 
    
    # Reseeding Infection
    
    if(reseed == TRUE){
      # reseeding generates IDs of new infections in the population, up to 5 units
      reseed_I <- reseeding(gen, threshold, frequency, introductions, S_index, I_index, cases_index)
      
      newI_index <- c(cases_index, reseed_I) # update case index vector with new cases from re-seeding infection
    }else if(reseed == FALSE){
      newI_index <- cases_index
    }
    
    new_disease_state[S_index] <- "S" # all S units remain in S state
    new_disease_state[newI_index] <- "I" # update ids of units in cases_index to "I"
  }
  
  
  # ------------
  ## VACCINATION ##
  # ------------
  
  if(length(V_timing) > 0 & p_vac > 0){
    
    if(any(gen == V_timing)){
      
      if(V_mode == "random"){
        
        newV_index <- sample(N, round(N*p_vac), replace = FALSE)
        
        new_disease_state[newV_index] <- "NV"
        
      }
      
      if(V_mode == "target" & heterogeneity == "fixed"){
        
        newV_index <- (1:N)[rev(order(fixed_effc))][1:round(p_vac*N)]
        
        new_disease_state[newV_index] <- "NV"
      }
      
      
      if(V_mode == "non_target" & heterogeneity == "fixed"){
        
        newV_index <- (1:N)[order(fixed_effc)][1:round(p_vac*N)]
        
        new_disease_state[newV_index] <- "NV"
        
      }
      
      #VTime[newV_index] <- 0
      
    }
  }
  
  return(new_disease_state)
}

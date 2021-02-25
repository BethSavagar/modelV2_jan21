
# Updated with vaccination and working as per 11.01.21


#newDiseaseState function takes as input the current state of the population, produces as output a vector containing the new state of the population


newDiseaseState <-
  function(disease_state, # the current disease state of units in the population
           N, # population size
           population_sample, # large or small, determines infection sampling strategy
           gen, # tracks which generation is being simulated
           
           R_period, # refractory period (number of timesteps spent immune)
           RTime, # tracks the number of timesteps spent in "R" state for each individual unit
           R0_estimate, # R0
           k_transmission, # dispersion parameter for transmission 
           k_exposure, # dispersion parameter for exposure, required to generate probability distribution when independent of tranmsission
           heterogeneity, # "fixed" or "variable" depending on whether unit transmission potential changes with time
           correlation, # relationship between transmission & exposure correlated postively (pos), negatively (neg) or unrelated (none)
           contact_list, # potential contacts for all units in the populatio 
           fixed_effc, # effective contact rate for units, required if heterogneeity is fixed
           fixed_exp,
           
           reseed, # reseed infection periodically during the first 10% of generations
           reseed_threshold, # threshold generation for reseeding events
           reseed_frequency, # frequency of reseeding events
           introductions, # number of infections introduced during reseeding
           
           V_prop, # proportion to be vacciantion in each round
           V_period, # duration of immunity from vaccination (by number of timesteps)
           V_timing, # vector defining timesteps for vaccination
           V_mode,  # mode of vacc (random, target, non-target)
           VTime # used to avoid problem of unused argument in VTime : a vector tracking hte number of generations each unit spends in V state
  ){
  
  S_index <- which(disease_state == "S") # ids of Susceptible individuals
  I_index <- which(disease_state == "I") # ids of Infected individuals
  R_index <- which(disease_state == "R") # ids of Recovered individuals
  V_index <- which(disease_state == "V") # ids of Vacc individuals
  
  # Vector to hold new disease state
  #new_disease_state <-
    #vector(length = length(disease_state)) # vector to hold disease state of units in population at next timestep
  new_disease_state <- disease_state # 28/01/21 create new vector to store new state of units in population
  
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
    
    newS_index <- which(RTime == R_period) # identify units which have been in R state for >= ~R_period, store in newS_index
    
    # Update new_disease_state vector
    # new_disease_state[R_index] <- "R" - commented out 28/01/21 due to change in new_DS vector <- disease_state (line 42)
    new_disease_state[newS_index] <- "S" # units identified in newS_index lose immunity
  }
  
  # Vaccinated units become susceptible again after V_period 
  
  if (length(V_index) > 0) {
    newS_index <- which(VTime == V_period) 
   
    # Update new_disease_state vector
    # new_disease_state[V_index] <- "V"  28/01/21 commented out due to line 42
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
    cases_index <- S_to_I(I_index, # ID of infected units
                          S_index, # ID of susceptible units
                          N, # population size
                          population_sample, # large or small, determines sampling of infections
                          R0_estimate,
                          k_transmission,
                          k_exposure,
                          heterogeneity, # fixed or variable, depending on whether unit transmission potential is constant
                          contact_list,
                          fixed_effc, # effective contact rate of units if fixed heterogeneity
                          fixed_exp
    ) 
    
    # Reseeding Infection
    
    if(reseed == TRUE){
      # reseeding generates IDs of new infections in the population, up to 5 units
      reseed_I <- reseeding(gen, # current generation
                            reseed_threshold, 
                            reseed_frequency, 
                            introductions, 
                            S_index, 
                            I_index, 
                            cases_index
                            )
      
      newI_index <- c(cases_index, reseed_I) # update case index vector with new cases from re-seeding infection
    
      }else if(reseed == FALSE){
        
      newI_index <- cases_index
    }
    
    #new_disease_state[S_index] <- "S" # commented out per line 42 on 28/01/21
    new_disease_state[newI_index] <- "I" # update ids of units in cases_index to "I"
  }
  
  
  # ------------
  ## VACCINATION ##
  # ------------
  
  if(length(V_timing) > 0 & V_prop > 0){
    
    if(any(gen == V_timing)){ # if the current generation is in the V_timing vector (vac schedule)
      
      # if above is true then conduct vaccination according to V_mode
      
      if(V_mode == "random"){ # random vaccination
        
        newV_index <- sample(N, #population
                             round(N*V_prop), # number of units to vaccinate
                             replace = FALSE # cannot vaccinate same unit twice in 1 round
                             )
        
        new_disease_state[newV_index] <- "NV" # NV state is used so that VTime can be updated in stateUpdate function
        
      }
      
      # targetted/selective vaccination is only possible if units have fixed transmission potential
      
      if(V_mode == "target" & heterogeneity == "fixed"){
        
        newV_index <- (1:N)[rev(order(fixed_effc))][1:round(V_prop*N)] # subset population for vac by reverse order of NBD, subset proportion corresponding to V_prop
        
        new_disease_state[newV_index] <- "NV"
      }
      
      
      if(V_mode == "non_target" & heterogeneity == "fixed"){
        
        newV_index <- (1:N)[order(fixed_effc)][1:round(V_prop*N)]
        
        new_disease_state[newV_index] <- "NV"
        
      }
    }
  }
  
  return(new_disease_state)
}

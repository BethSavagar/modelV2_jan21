

# The S_to_I function takes as arguments: the ids of infected units (I_index), the ids of susceptible units (S_index), the population size (N), an estimate of the overdispersion of effective contacts (k_estimate), the average number of effective contacts per infected unit (R_estimate), an argument to define whether the effective contact potential of units is fixed or variable across timesteps (heterogeneity) and a vector of effective contact potentials for each unit (fixed_effc, for use when heterogeneity == "fixed")


S_to_I <- function(I_index,
                   S_index,
                   N,
                   population_sample,
                   R0_estimate,
                   k_transmission,
                   heterogeneity, # defines whether the effective contact potential of units is fixed or variable across timesteps
                   contact_list,# list of possible contacts for each unit (with self removed)
                   fixed_effc # stores the effective contact potential of units, required when heterogeneity is fixed across timesteps
                   ) {
  # In 2 steps: (i) Simulate the number of effective contacts per infected unit, (ii) Generate secondary cases (from effective contacts with S units)
  
  # ------------------------------------------------------------------------------------------
  ## SIMULATE NUMBER OF EFFECTIVE CONTACTS PER INFECTED UNIT
  # ------------------------------------------------------------------------------------------
  # Calculate the number of effective contacts made by each infected unit in a given timestep
  # Stored initially in effc_num vector
  
  ## Heterogeneity == FIXED ##
  # When heterogeneity == "fixed" the effective contact potential of each unit is constant across timesteps and is stored in the vector fixed_effc
  
  if (heterogeneity == "fixed" |
      heterogeneity == "homogeneous" ) { 
    
    # define effc_num vector to store the effc value (stored in fixed_effc) for infected individuals 
    effc_num <- fixed_effc[I_index] # subsetting effc_fixed vector with I_index
    
    ## Heterogeneity == VARIABLE ##   
    # When heterogeneity == "variable" the effective contact potential of each unit can vary across timesteps (i.e. stochastic), redrawn from the Negative Binomial Distribution each timestep.  
    
  } else if (heterogeneity == "variable") { 
    
    I0 <- length(I_index) # the number of currently infected cases
    
    # define vector effc_raw to store the number of effective contacts for each infected unit (drawn from a Negative Binomial Distribution)
    effc_raw <- rnbinom(I0,   # I0 is the number of infected units in the population (defines number of elements in effc_raw vector)
                        size = k_transmission, # dispersion parameter of effective contacts
                        mu = R0_estimate) # average number of effective contacts for a unit in the population
    
    ## Effective Contact Verification ##
    # verify that no unit produces >= N effective contacts (an impossibility) using the effcChecker function (see above)
    effc_num <- effcChecker(effc_raw, 
                            N, 
                            k_transmission, 
                            R0_estimate) 
    
  }
  
  # ------------------------------------------------------------------------------------------
  ## GENERATE IDs of SECONDARY CASES
  # ------------------------------------------------------------------------------------------
  # 2 step process: (i) Contact IDs (all units), (ii) Case IDs (only S units)
  
  # If the number of effective contacts (effc) is non-zero simulate contacts between units in the population by random sampling
  # An 'effective contact' can occur between an infected unit and any other unit in the population, excluding itself (ie units in S, I or R state can be contacted). See effcIdGenerator   
  # A contacted unit only becomes infected if it is in the S state (see line 287)
  
  ## Clean the effc vector ##
  
  # subset the 'effc_num' vector to contain only non-zero effective contacts
  effc <- effc_num[effc_num > 0] # this will be used to simulate contacts between units in the population using random sampling.
  
  ## Store Infected Unit (with non-zero contacts) IDs: ##
  
  # A vector 'effc_id' is defined which store the ids of infected units which generated more than 0 effective contacts (effc_num > 0)
  # this will be used to ensure that no self-self contacts are simulated 
  effc_id <- I_index[effc_num > 0]
  
  # generate a list of potential contacts for the I_index cases - to remove self-self contacts
  if(population_sample == "small"){
    possible_contacts <- contact_list[effc_id]
  }else{
    possible_contacts <- vector()
  }
  
  
  ## Generate Contact IDs ##
  
  # If >0 effective contacts are made (i.e. length(effc)>0) use effcIdGenerator to produce list of effectively contacted unit ids 
  
  if (length(effc) > 0) {
    
    # effcIdGenerator produces a list of effectively contacted unit ids (stored in contact_ids vector) using random sampling of population IDs (see effcIdGenerator code)
    contact_ids <- effcIdGen2(possible_contacts, # potential contacts of infected units (ie removing infected unit ID)
                              effc, # vector containing number of effective contacts per infected unit
                              N, # population size
                              population_sample # determines the sample strategy depending on large/small population
                              )  
    
    ## Generate Case IDs ##
    
    # Define the cases_index vector to store the ids of newly infected units (effectively contacted S units)
    # Cases are generated from effective contacts between an infected unit and a susceptible unit
    # Within the vector of contacted unit ids (contact_ids), identify those ids which are in the Susceptible state (S_index)
    
    cases_index <- contact_ids[contact_ids %in% S_index]
    
  } else{
    
    #if no secondary cases are produced then case_index vector <- 0
    cases_index <- numeric(0)
  }
  
  return(cases_index)
}

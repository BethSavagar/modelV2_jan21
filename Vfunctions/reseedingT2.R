



reseeding <-
  function(gen,
           reseed_threshold,
           reseed_frequency,
           introductions,
           S_index,
           I_index, 
           cases_index) {
    
    
    if (gen %% reseed_frequency == 0 & # every X generations
        gen < reseed_threshold) { # if current generation is below threhsold
      
      # Define sample population for introductions
      
      if (length(cases_index) > 0) {
        # remove new infected cases from S_index for sampling
        S_pop <- S_index[-cases_index] 
        }else{
          S_pop <- S_index # if no new infected cases then all of S_index used for sampling.
        }
      
      # Sample population (ids in vector S_pop)
      if (length(S_pop) > introductions) {
        # if the number of units in S_index is > than the number of introductions specified, then sample from S_index
        reseed_I <- sample(S_pop, # susceptible population to sample from 
                           introductions, # number of introductions
                           replace = FALSE) # randomly sample from susceptible population (excluding ID of new I cases)
      } else {
        # if the number of units in S_index are < the number of introductions, make all remaining S into I
        reseed_I <- S_pop
      }
      
    } else{
      
      reseed_I <- numeric(0)
    }
    
    return(reseed_I)
    
  }
  
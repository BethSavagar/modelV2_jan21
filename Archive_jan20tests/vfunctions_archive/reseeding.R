



reseeding <-
  function(gen,
           threshold,
           frequency,
           introductions,
           S_index,
           I_index, 
           cases_index) {

    
    if (gen %% frequency == 0 & # every X generations
        gen < threshold) { # if current generation is below threhsold
      
      # Define sample population for new reseeding of infections
      if (length(cases_index) > 0) {
        Z <- S_index[-cases_index] #remove new infected cases from S_index for sampling
      } else{
        Z <- S_index # if no new infected cases then all of S_index used or sampling.
      }
      
      
      # Sample population (ids in vector Z)
      if (length(Z) > introductions) {
        # if the number of units in S_index is > than the number of introductions specified, then sample from S_index
        reseed_I <- sample(Z,
                           introductions,
                           replace = FALSE) # randomly sample from susceptible population (excluding ID of new I cases)
      } else{
        # if the number of units in S_index are < the number of introductions, make all remaining S into I
        reseed_I <- Z
      }
    } else{
      reseed_I <- numeric(0)
    }
    
    return(reseed_I)
    
  }
  
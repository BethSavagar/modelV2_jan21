effcChecker <- function(effc_num, N, k_estimate, R_estimate) {
  # First identify whether any infected units produced >= N effective contacts
  # The vector effc_N_id stores the position of infected units which caused >= N effective contacts in the effc_num vector
  effc_N_id <- which(effc_num >= N)
  
  # Check whether there are any units which have >= N effective contacts ( length(effc_N_id) > 0)
  # Run the loop described above
  
  if (length(effc_N_id) > 0) {
    for (i in effc_N_id) {
      # iterate over the vector efffc_N_id which contains positions of units generating >= N effective contacts in the effc_num vector
      
      effc_i <- effc_num[i] # assign variable effc_i to the number of effective contacts for the indexed unit (position i in effc_num vector)
      
      # Use a while loop to re-run the negative binomial simulation of effective contacts for the indexed unit (effc_i) until the number of effective contacts is < N.
      while (effc_i >= N) {
        new_i <- rnbinom(1, # simulate for single infected unit
                         size = k_estimate, # dispersion parameter, as above
                         mu = R_estimate) # average number of effective contacts, as above
        effc_i <- new_i # update effc_i within while loop
      }
      
      effc_num[i] <- effc_i # Update the effc_num vector with the new number of effective contacts for the infected unit.
    }
  }
  return(effc_num)
}


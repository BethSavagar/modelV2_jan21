effcIdGen2 <- function(contact_list, # possible contacts
                       effc, # vector of effective contacts per infected unit
                       N, # population size
                       population_sample                       
                       ){
  
  
  
  if(population_sample == "small"){ # if population is relatively small, use correct sampling procdure, with ID of infected unit removed
    
    C_list <- lapply(1:length(contact_list), 
                     function(x) sample(contact_list[[x]], # possible contacts for indexed case
                                        effc[x], # number of contacts for indexed case
                                        replace = FALSE) # infected individual cannot contact the same unit twice
                     )
    
    contact_ids <- unique(unlist(C_list)) # contact_ids contains the list of all effectively contacted units
    
  } else if (population_sample == "large"){ # if population is relatively large, the effect of keeping index case ID in sampling population is negligible
    
    effcSum <- sum(effc)
    contact_ids <- unique(sample(N, # population to sample from
                                 effcSum, # total effective contacts by all infected uunits
                                 replace=TRUE # same unit can be contacted by multiple infected units
                                 ) 
                          )
    
  }
  
  return(contact_ids)
  
}




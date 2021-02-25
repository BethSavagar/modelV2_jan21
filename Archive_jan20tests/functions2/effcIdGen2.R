effcIdGen2 <- function(contact_list2, Vsize, N, population){
  vSum <- sum(Vsize)
  
  if(population == "small"){ # if populaiton is relatively small, use correct sampling procdure, with ID of infected unit removed
    
    L <- lapply(1:length(contact_list2), 
                function(x) 
                  sample(contact_list2[[x]], Vsize[x], replace = FALSE))
    
    contact_ids <- unique(unlist(L))
    
  }else if(population == "large"){ # if population is relatively large, the effect of keeping index case ID in sampling population is negligible
    
  contact_ids <- unique(sample(N,
                               vSum,
                               replace=TRUE))
  
}
  return(contact_ids)
  
}




library(ggplot2)
library(dplyr)
# with non endmeic sims excluded


endemicity <- function(Tmat){
  
  P_endemic <- sum(Tmat[,ncol(Tmat)]!=0) / nrow(Tmat) 

  return(P_endemic)
}


kvec <- c(0.1,
          0.3,
          0.5,
          0.7,
          1, 
          Inf)
#kvec <- c(0.1, 0.3, 0.5, 0.7, 1, Inf)
MasterList <- list(k0.1List,
                   k0.3List,
                   k0.5List,
                   k0.7List,
                   k1List,
                   kInfList)
#MasterList <- list(k0.1List,k0.3List,k0.5List,k0.7List,k1List, kInfList)

endemicity_df <- vector()

for(k in 1:length(kvec)){
  
  kval <- kvec[k]
  setList <- MasterList[[k]]
  
  pEnd <- sapply(setList, endemicity)
  
  kEnd <- cbind(kval,pars,pEnd)
  
  endemicity_df <- rbind(endemicity_df, kEnd)
}



ggplot(data = endemicity_df, aes(x = R0, y = pEnd))+
  geom_line(aes(group = factor(kval), colour = factor(kval)), position=position_jitter(w=0.02, h=0))+
  facet_wrap(~Rperiod)+
  theme_bw()+
  expand_limits(x = 0, y = 0)


ggplot(data = endemicity_df, aes(x = factor(kval), y = pEnd))+
  geom_line(aes(group = factor(R0), colour = factor(R0)))+
  facet_wrap(~Rperiod)+
  theme_bw()

endemicity_df %>% filter(pEnd != 0) %>%
  ggplot(aes(x = R0, y = pEnd))+
  geom_point(aes(group = factor(kval), colour = factor(kval)), position=position_jitter(w=0.02, h=0))+
  facet_wrap(~Rperiod)+
  theme_bw()+
  geom_line(aes(group = factor(kval), colour = factor(kval)), position=position_jitter(w=0.02, h=0))
  


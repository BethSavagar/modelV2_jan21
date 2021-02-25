library(ggplot2)


kvec <- c(0.1, 0.3, 0.5, 0.7, 1, Inf)
#kvec <- c(0.1, 0.3, 0.5, 0.7, 1, Inf)
MasterList <- list(k0.1List,
                   k0.3List,
                   k0.5List,
                   k0.7List,
                   k1List,
                   kInfList)
#MasterList <- list(k0.1List,k0.3List,k0.5List,k0.7List,k1List, kInfList)

medInf <- vector()

for(k in 1:length(kvec)){
  
  kval <- kvec[k]
  setList <- MasterList[[k]]
  
  EQ_medInf <- sapply(setList, function(x)
    equilibrium <- median(x[, 900:1000])
  )
  
  k_medInf <- cbind(kval, pars, EQ_medInf)
  
  medInf <- rbind(medInf, k_medInf)
}


ggplot(medInf, aes(x = as.factor(Rperiod), y = R0, fill = EQ_medInf))+
  geom_tile()+
  scale_fill_distiller(palette = "Blues", name = "Median Incidence")+
  xlab("Refractory Period")+
  scale_y_continuous(breaks = c(seq(0,5,0.5)),
                     labels =  c(as.character(seq(0,5,0.5))))+
  facet_wrap(~kval, labeller = label_both)+
  theme_bw()

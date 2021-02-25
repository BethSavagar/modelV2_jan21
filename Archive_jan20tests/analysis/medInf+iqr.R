library(ggplot2)
library(dplyr)
# with non endmeic sims excluded

kvec <- c(0.1, 0.5, 1)
#kvec <- c(0.1, 0.3, 0.5, 0.7, 1, Inf)
MasterList <- list(k0.1List,k0.5List,k1List)
#MasterList <- list(k0.1List,k0.3List,k0.5List,k0.7List,k1List, kInfList)


EpiStat <- function(mat) {
  if (median(mat) == 0) {
    M <- median(mat)
    iqr <- IQR(mat)
    q1 <- quantile(mat)[2]
    q3 <- quantile(mat)[4]
  } else{
    vec <- as.vector(mat[(mat != 0)])
    M <- median(vec)
    iqr <- IQR(vec)
    q1 <- quantile(vec)[2]
    q3 <- quantile(vec)[4]
  }
  sumStat <- c(M, iqr, q1, q3)
  return(sumStat)
}

MedIQR <- vector()

for(k in 1:length(kvec)){
  
  kval <- kvec[k]
  setList <- MasterList[[k]]
  
  Eq <- lapply(setList, function(x)
    eq <- x[, 900:1000]
  )
  
  sumStat <- t(sapply(Eq, EpiStat))
  colnames(sumStat) <- c("med", "iqr", "q1", "q3")
  parStat <- cbind(kval, pars, sumStat)
  
  MedIQR <- rbind(MedIQR, parStat)
}

MedIQR <- as.data.frame(MedIQR)


ggplot(data = MedIQR, 
       mapping = aes(x = as.factor(R0), 
                     y = med, 
                     ymin = q1, 
                     ymax = q3)) +
  geom_pointrange(aes(group = factor(kval), colour = factor(kval)))+
  facet_wrap(~Rperiod)+
  theme_bw()


ggplot(data = MedIQR, 
       mapping = aes(x = as.factor(R0), 
                     y = med 
                      )) +
  geom_point(aes(group = factor(kval), colour = factor(kval)))+
  geom_errorbar(aes(ymin = q1, 
                ymax = q3))+
  facet_grid(kval~Rperiod, labeller = label_both)+
  xlab("R0")+
  ylab("median incidence G900-1000")+
  theme_bw()


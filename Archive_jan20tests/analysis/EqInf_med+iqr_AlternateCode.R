# median incidence at equilibrium:


# first filter for last 100 generations


endStats <- function(mat){
  eq <- mat[, (ncol(mat)-100) :ncol(mat)]
  id <- which(eq[,ncol(eq)] != 0) # which simulations became endemic 
  if(length(id) == 0){
    M <- 0 
    iqr <- 0
    lq <- 0
    uq <- 0
  }else if(length(id)>0){
    mEnd <- as.vector(eq[id,]) # simulations which become endemic, subset for last 100 generations
    M <- median(mEnd)
    iqr <- IQR(mEnd)
    lq <- quantile(mEnd)[2]
    uq <- quantile(mEnd)[4]
  }
  sumStat <- c(Med = M, iqr = iqr, lq = lq, uq = uq)
  return(sumStat)  
  }

EqStats <- t(sapply(kInfList, endStats))
pars_eq <- cbind(pars, EqStats)


ggplot(data = pars_eq, aes(x = R0, y = Med))+
  geom_line(aes(group = factor(Rperiod), colour = factor(Rperiod)))+
  xlab("R0")+
  ylab("Median incidence at Equilibrium")+
  scale_colour_discrete("Refractory Period")+
  theme_bw()

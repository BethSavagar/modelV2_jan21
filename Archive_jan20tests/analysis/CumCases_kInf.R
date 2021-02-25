
library(dplyr)
library(ggplot2)

example <- kInfList[[13]]

# cumulative cases

med_CC <- sapply(kInfList, function(x)
  median(apply(x,1,sum))) # median cumulative cases over 100 simulations for each parameter combination

iqr_CC <- sapply(kInfList, function(x)
  IQR(apply(x,1,sum)))

pars_CC <- cbind(pars, med_CC, iqr_CC)

# filter just simulations where Rperiod = 11

pars_CC %>% filter(Rperiod == 1) %>%
  ggplot(aes(x = R0, y = med_CC))+
  geom_point()+
  xlab("R0")+
  ylab("Cumulative Cases (100G)")+
  theme_bw()

  
  ggplot(data = pars_CC, aes(x = R0, y = med_CC))+
    geom_line(aes(group = factor(Rperiod), colour = factor(Rperiod)))+
    xlab("R0")+
    ylab("Cumulative Cases (100G)")+
    scale_colour_discrete("Refractory Period")+
    theme_bw()
  